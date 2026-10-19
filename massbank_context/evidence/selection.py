"""出典別の候補を消さずに関連付ける。"""
from dataclasses import asdict


def candidate_dicts(bundle, matches):
    candidates = [asdict(c) for c in bundle.annotations] + list(matches.get('candidates', []))
    for c in candidates:
        c['same_structure_as'] = sorted(other['candidate_id'] for other in candidates
                                         if other['candidate_id'] != c['candidate_id']
                                         and c['inchikey'] and c['inchikey'] == other['inchikey'])
    return candidates
