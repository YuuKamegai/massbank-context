"""出典別の候補を消さずに関連付ける。"""
from dataclasses import asdict

# 上流でMS/MSの対応付け自体に失敗した理由。「MS/MSなし」とは区別する。
MAPPING_FAILURES = ('precursor_mismatch', 'metadata_missing', 'out_of_range')


def select_ids(all_ids, requested):
    ids = sorted(all_ids if requested is None else set(requested))
    if set(ids) - set(all_ids):
        raise ValueError('未知のIDが指定されています')
    return ids


def mapping_failure(bundle, feature_id: str) -> str | None:
    """上流が記録したMS/MS対応付けの失敗理由。失敗でなければNone。"""
    status = bundle.manifest.get('provenance', {}).get('msms_link_statuses', {}).get(feature_id, {})
    return status.get('status') if status.get('status') in MAPPING_FAILURES else None


def links_by_feature(bundle) -> dict:
    grouped = {}
    for link in bundle.spectrum_links:
        grouped.setdefault(link.feature_id, []).append(link)
    return grouped


def candidate_dicts(bundle, matches):
    candidates = [asdict(c) for c in bundle.annotations] + list(matches.get('candidates', []))
    groups = {}
    for c in candidates:
        if c['inchikey']:
            groups.setdefault((c['feature_id'], c['inchikey']), []).append(c['candidate_id'])
    for c in candidates:
        c['same_structure_as'] = sorted(cid for cid in groups.get((c['feature_id'], c['inchikey']), [])
                                         if cid != c['candidate_id'])
    return candidates
