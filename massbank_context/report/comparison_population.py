"""候補の曖昧さを公開した統計用母集団への射影。"""
from dataclasses import dataclass, asdict, is_dataclass
from massbank_context.experiment.bundle import INCHIKEY_PATTERN
from massbank_context.stats.activity import pathway_activity_rows


@dataclass(frozen=True)
class PopulationProjection:
    compounds: list[dict]
    exclusions: dict
    representatives: dict
    summary: dict


def build_population(bundle, candidates, match_status):
    candidates = [asdict(c) if is_dataclass(c) else c for c in candidates]
    grouped, exclusions = {}, {}
    linked = {s.feature_id for s in bundle.spectrum_links}
    for feature in bundle.features:
        fid = feature.feature_id
        active = [c for c in candidates if c['feature_id'] == fid and c['status'] == 'active']
        status = match_status.get(fid, {})
        reason = None
        mapping = bundle.manifest.get('provenance', {}).get('msms_link_statuses', {}).get(fid, {})
        if mapping.get('status') in ('precursor_mismatch', 'metadata_missing', 'out_of_range'):
            reason = 'msms_mapping_failed'
        elif fid in linked and status.get('status') != 'matched':
            reason = 'matching_incomplete'
        elif status.get('truncated'):
            reason = 'candidate_truncated'
        elif any(getattr(feature, k) is None for k in ('log2fc', 'p_value', 'q_value', 'mean_a', 'mean_b')):
            reason = 'missing_comparison'
        elif not active:
            reason = 'no_candidate'
        elif any(not c['inchikey'] or not INCHIKEY_PATTERN.fullmatch(c['inchikey']) for c in active):
            reason = 'unsupported_identifier'
        elif len({c['inchikey'] for c in active}) != 1:
            reason = 'candidate_conflict'
        if reason:
            exclusions[fid] = dict(reason=reason, candidate_ids=[c['candidate_id'] for c in active])
            continue
        ik = active[0]['inchikey']
        grouped.setdefault(ik, []).append(feature)
    compounds, representatives = [], {}
    for ik, features in sorted(grouped.items()):
        signs = {(f.log2fc > 0) - (f.log2fc < 0) for f in features} - {0}
        if len(signs) > 1:
            for f in features:
                exclusions[f.feature_id] = dict(reason='direction_conflict', inchikey=ik)
            continue
        chosen = sorted(features, key=lambda f: (-int(f.significant), -abs(f.log2fc), f.feature_id))[0]
        item = dict(inchikey=ik, log2fc=chosen.log2fc, significant=chosen.significant,
                    feature_id=chosen.feature_id, feature_ids=sorted(f.feature_id for f in features))
        compounds.append(item)
        representatives[ik] = item
    summary = dict(population_complete=bundle.population_complete,
                   n_features_total=bundle.manifest['population']['n_features_total'],
                   n_features_available=len(bundle.features), n_excluded=len(exclusions),
                   n_compounds=len(compounds),
                   n_pending_matching=sum(x['reason'] == 'matching_incomplete' for x in exclusions.values()),
                   rule='activeのfull InChIKeyが一意、打切りなし。逆方向重複を除外し、有意優先・絶対log2fc最大で代表化。')
    return PopulationProjection(compounds, exclusions, representatives, summary)


def compute_pathways(projection, lookups, organism: str, min_members: int = 3):
    if type(min_members) is not int or min_members < 1:
        raise ValueError('min_membersは1以上です')
    base = dict(projection.summary, exclusions=projection.exclusions,
                representatives=projection.representatives, organism=organism)
    if not projection.summary['population_complete']:
        return dict(base, tested=False, reason='population_incomplete')
    if projection.summary['n_pending_matching']:
        return dict(base, tested=False, reason='matching_incomplete')
    if not organism:
        return dict(base, tested=False, reason='organism_required')
    memberships, pending = {}, []
    for c in projection.compounds:
        relevant = [l for l in lookups if l.get('inchikey') == c['inchikey']]
        if not relevant or any(l['status'] not in ('resolved', 'empty') for l in relevant):
            pending.append(c['inchikey'])
            continue
        seen, hits = set(), []
        for lookup in relevant:
            for hit in lookup['hits']:
                if hit.get('taxid') == organism and hit['pathway_id'] not in seen:
                    seen.add(hit['pathway_id'])
                    hits.append(hit)
        if hits:
            memberships[c['inchikey']] = hits
    if pending:
        return dict(base, tested=False, reason='lookup_incomplete', pending_inchikeys=pending)
    background = [c for c in projection.compounds if c['inchikey'] in memberships]
    rows, summary = pathway_activity_rows(background, memberships, min_members)
    for row in rows:
        row['feature_ids'] = sorted({fid for c in background
                                     if any(h['pathway_id'] == row['pathway_id'] for h in memberships[c['inchikey']])
                                     for fid in c['feature_ids']})
    return dict(base, tested=any(r['tested'] for r in rows), rows=rows,
                statistics=summary, n_lookup_complete=len(projection.compounds),
                n_with_species_pathways=len(background), n_without_species_pathways=len(projection.compounds)-len(background),
                reason=None if rows else 'no_species_pathways',
                note='濃縮・増減は経路活性やフラックスを直接証明しません。')
