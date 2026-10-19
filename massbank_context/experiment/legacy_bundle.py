"""旧TSVの欠落行を復元せずに取り込む。"""
import hashlib
from pathlib import Path
from massbank_context.experiment.contract import parse_differential
from massbank_context.core.evidence_ids import stable_id
from massbank_context.core.evidence_models import EvidenceBundle, Feature, Candidate


def import_legacy(path: Path) -> EvidenceBundle:
    path = Path(path).resolve()
    raw = path.read_bytes()
    differential = parse_differential(raw.decode('utf-8-sig'))
    aid = stable_id('analysis', ['legacy', hashlib.sha256(raw).hexdigest()])
    manifest = dict(contract_name='legacy-differential', contract_version=1,
                    comparison=dict(differential.meta), population=dict(
                        n_features_total=differential.n_features_total,
                        population_complete=False), files=[])
    features, annotations = [], []
    ids = set()
    for row in differential.rows:
        fid = row['spot_id']
        if not fid or fid in ids:
            raise ValueError('旧TSVのspot_idが空または重複しています')
        ids.add(fid)
        features.append(Feature(feature_id=fid, **{key: value for key, value in row.items()
                                                   if key in Feature.__dataclass_fields__ and key != 'feature_id'}))
        annotations.append(Candidate(stable_id('candidate', [aid, fid]), fid,
                                     'legacy_annotation', fid, name=row['name'],
                                     inchikey=row['inchikey'], provenance={'legacy': True}))
    return EvidenceBundle(aid, manifest, tuple(features), tuple(annotations), (), {}, False, path)
