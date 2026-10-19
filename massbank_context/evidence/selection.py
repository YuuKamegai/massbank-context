"""出典別の候補を消さずに関連付ける。"""
from dataclasses import asdict


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
