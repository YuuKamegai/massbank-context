"""全体概要と、版を固定した特徴量別の根拠を提示する。"""
import base64
import json
from collections import Counter
from dataclasses import asdict
from massbank_context.core.evidence_ids import stable_id, canonical_json
from massbank_context.evidence.selection import candidate_dicts
from massbank_context.evidence.matching import select_ids
from massbank_context.report.comparison_population import build_population


def get_evidence(store, analysis_id: str, *, feature_ids=None, cursor=None, limit=50):
    if type(limit) is not int or not 1 <= limit <= 200:
        raise ValueError('limitは1〜200です')
    bundle, revision, state = store.read_view(analysis_id)
    matches = state.get('matches', {})
    lookup_state = state.get('lookups', {})
    context = state.get('context', {})
    candidates = candidate_dicts(bundle, matches)
    ids = select_ids([f.feature_id for f in bundle.features], feature_ids)
    filter_hash = stable_id('filter', [canonical_json(ids)])
    offset, evidence_offset, detail_page = 0, 0, False
    cur = {}
    if cursor:
        try:
            cur = json.loads(base64.urlsafe_b64decode(cursor.encode()))
        except (ValueError, TypeError, UnicodeError) as exc:
            raise ValueError('invalid_cursor') from exc
        if not isinstance(cur, dict):
            raise ValueError('invalid_cursor')
        if (cur.get('analysis_id') != analysis_id or cur.get('revision') != revision
                or cur.get('filter') != filter_hash or cur.get('limit') != limit):
            raise ValueError('stale_cursor: 解析・条件または根拠版が変わりました')
        offset = cur.get('offset')
        if type(offset) is not int or offset < 0 or offset > len(ids):
            raise ValueError('invalid_cursor')
        detail_page = cur.get('mode') == 'evidence'
        evidence_offset = cur.get('evidence_offset', 0)
        if type(evidence_offset) is not int or evidence_offset < 0:
            raise ValueError('invalid_cursor')
    shown = set(ids[offset:offset+limit])
    selected = [c for c in candidates if c['feature_id'] in shown]
    candidate_ids = {c['candidate_id'] for c in selected}
    lookups_by_id = {l['candidate_id']: l for l in lookup_state.get('lookups', [])}
    lookups = []
    for cid in sorted(candidate_ids):
        raw = lookups_by_id.get(cid, dict(candidate_id=cid, status='not_queried'))
        lookups.append(dict({k:v for k,v in raw.items() if k != 'hits'}, n_hits=len(raw.get('hits', []))))
    evidence = [e for e in lookup_state.get('evidence', []) if e['candidate_id'] in candidate_ids]
    for c in selected:
        for record in lookup_state.get('catalog', []):
            if c['inchikey'] and record['inchikey'] == c['inchikey'] and (
                    record['accession'] == (c['source_id'] if c['source'] == 'massbank' else None)):
                evidence.append(dict(record, evidence_id=stable_id('evidence', [c['candidate_id'], canonical_json(record)]),
                                     candidate_id=c['candidate_id'], feature_id=c['feature_id'],
                                     source='massbank_rdf_index', provenance=lookup_state.get('catalog_provenance', {})))
    counts = Counter(l['status'] for l in lookup_state.get('lookups', [])
                     if l['candidate_id'] in {c['candidate_id'] for c in candidates})
    linked = {l.feature_id for l in bundle.spectrum_links}
    annotated = {c['feature_id'] for c in candidates if c['status'] == 'active'}
    total = len(bundle.features)
    summary = dict(n_features_total=bundle.manifest['population']['n_features_total'],
                   n_features_available=total, population_complete=bundle.population_complete,
                   n_significant=sum(f.significant for f in bundle.features),
                   n_up=sum(f.log2fc is not None and f.log2fc > 0 for f in bundle.features),
                   n_down=sum(f.log2fc is not None and f.log2fc < 0 for f in bundle.features),
                   n_missing_log2fc=sum(f.log2fc is None for f in bundle.features),
                   n_with_candidates=len(annotated), n_with_msms=len(linked),
                   n_candidates=len(candidates), lookup_statuses=dict(counts),
                   n_not_queried=len(candidates)-sum(counts.values())+counts.get('not_queried',0))
    projection = build_population(bundle,candidates,matches.get('features',{}))
    eligibility = dict(projection.exclusions)
    for ik, representative in projection.representatives.items():
        for fid in representative['feature_ids']:
            eligibility[fid] = dict(reason='eligible',inchikey=ik,
                                    representative_feature_id=representative['feature_id'])
    def encode_cursor(**extra):
        return base64.urlsafe_b64encode(canonical_json(dict(analysis_id=analysis_id, revision=revision,
            filter=filter_hash, limit=limit, **extra)).encode()).decode()
    next_cursor = encode_cursor(offset=offset+limit) if not detail_page and offset+limit < len(ids) else None
    evidence.sort(key=lambda e:e['evidence_id'])
    n_evidence = len(evidence)
    if evidence_offset > n_evidence:
        raise ValueError('invalid_cursor')
    next_evidence = encode_cursor(offset=offset, mode='evidence', evidence_offset=evidence_offset+64) if evidence_offset+64 < n_evidence else None
    return dict(analysis_id=analysis_id, evidence_revision=revision,
                context_revision=context.get('context_revision'), context=context,
                comparison=bundle.manifest['comparison'], summary=summary,
                population_derivation=dict(n_features_before=bundle.manifest['population'].get('n_features_before'),
                    exclusions=bundle.manifest['population'].get('exclusions', []),
                    caveats=(bundle.manifest['population'].get('preprocessing_report') or {}).get('caveats', []),
                    derivation_status=bundle.manifest['population'].get('derivation_status','not_available')),
                features=[dict(asdict(f),statistical_eligibility=eligibility.get(f.feature_id))
                          for f in bundle.features if f.feature_id in shown],
                candidates=selected, evidence=evidence[evidence_offset:evidence_offset+64], lookups=lookups,
                matching={fid: matches.get('features', {}).get(fid, dict(status='not_matched' if fid in linked else 'no_msms'))
                          for fid in sorted(shown)},
                scope=lookup_state.get('scope', ['massbank_metadata', 'pubchem_pathways']),
                truncated=next_cursor is not None, next_cursor=next_cursor,
                detail_page=detail_page, n_evidence=n_evidence, evidence_truncated=next_evidence is not None,
                next_evidence_cursor=next_evidence,
                n_selected=len(ids))
