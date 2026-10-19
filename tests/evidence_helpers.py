"""研究データを含まない根拠バンドルの合成入力。"""
import csv
import hashlib
import json
from pathlib import Path

KEY_A = 'AAAAAAAAAAAAAA-BBBBBBBBBB-C'
KEY_B = 'CCCCCCCCCCCCCC-DDDDDDDDDD-E'
KEY_CONFLICT = 'EEEEEEEEEEEEEE-FFFFFFFFFF-G'


def make_bundle(root: Path, *, case='mixed') -> Path:
    root.mkdir(parents=True, exist_ok=True)
    files = []
    features = []
    annotations = []
    for i in range(8):
        features.append(dict(feature_id=f'f{i}', spot_id=str(i), mz=200+i,
                             rt=1, log2fc='' if i == 7 else (-1 if i == 6 else 1),
                             p_value=.01, q_value=.02, mean_a=10, mean_b=20,
                             significant='true' if i % 2 == 0 else 'false'))
        if i != 4:
            annotations.append(dict(annotation_id=f'a{i}', feature_id=f'f{i}',
                                    source='msdial', name='名前\t候補',
                                    inchikey='' if i == 3 else KEY_CONFLICT if i in (5, 6) else KEY_A,
                                    rank=1, provenance_json='{}', curation_json='{}',
                                    status='active'))
    annotations.append(dict(annotations[0], annotation_id='other', inchikey=KEY_B))
    for role, rows, columns in (
        ('features', features, list(features[0])),
        ('annotations', annotations, list(annotations[0])),
        ('spectrum_links', [], ['feature_id', 'spectrum_file', 'spectrum_id',
                               'source_analysis_id', 'mapping_evidence']),
    ):
        path = root / f'{role}.tsv'
        with path.open('w', encoding='utf-8-sig', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        files.append(dict(role=role, path=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    manifest = dict(contract_name='massbank-context-evidence', contract_version=2,
                    source_analysis_id='synthetic-run', exporter={'name':'fixture','version':'1'},
                    exported_at='2026-10-19T00:00:00Z',
                    comparison=dict(group_a='A',group_b='B',n_a=3,n_b=3,
                                    log2fc_sign='positive means group_b is higher',
                                    method='Welch', correction='BH', q_threshold=.05,
                                    log2fc_threshold=1, preprocess={}),
                    population=dict(n_features_total=8, population_complete=True,
                                    definition='全比較対象', exclusions=[]), files=files)
    path = root / 'manifest.json'
    path.write_text(json.dumps(manifest, ensure_ascii=False), encoding='utf-8')
    return path


def change_manifest(path, **updates):
    data = json.loads(path.read_text(encoding='utf-8'))
    data.update(updates)
    path.write_text(json.dumps(data), encoding='utf-8')


def rewrite_table(manifest_path, role, change):
    data = json.loads(manifest_path.read_text(encoding='utf-8'))
    file = next(f for f in data['files'] if f['role'] == role)
    path = manifest_path.parent / file['path']
    with path.open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        columns, rows = reader.fieldnames, list(reader)
    change(rows)
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    file['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(data), encoding='utf-8')
