"""外来の分子の走査のための、研究データを含まない v3 根拠バンドルの合成入力。"""
import csv
import hashlib
import json
from pathlib import Path

TEST_SAMPLES = [f'T{i}' for i in range(1, 7)]
REFERENCE_SAMPLES = [f'R{i}' for i in range(1, 7)]

# feature_id, mz, MS/MS のピーク（None なら MS/MS 無し）, 検出数（test, reference）, 値（test）, 値（reference の検出時）
_FEATURES = [
    ('k_parent', 238.0990, [125.0153, 207.0571, 220.0888, 238.0990], (6, 0), [100, 200, 300, 400, 500, 600]),
    ('k_hnk', 240.0786, [125.0153, 89.0386, 222.0680], (6, 0), [10, 20, 30, 40, 50, 60]),
    ('k_frag', 220.0885, None, (6, 0), [5] * 6),
    ('x_msms', 500.3000, [125.0153, 207.0571], (4, 1), [7] * 6),
    ('mass_only', 254.0942, None, (6, 6), [9] * 6),
    ('endo', 300.2000, None, (6, 6), [9] * 6),
    ('no_mz', None, None, (6, 6), [9] * 6),
]


def _write_tsv(path: Path, columns, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def make_v3_bundle(root: Path, *, version: int = 3, polarity: str | None = 'positive',
                   detection: str = 'known', with_samples: bool = True) -> Path:
    """v3（または version=2 の v2）の根拠バンドルを書き、manifest のパスを返す。"""
    root.mkdir(parents=True, exist_ok=True)
    features, links, msp = [], [], []
    for fid, mz, peaks, _det, values in _FEATURES:
        features.append(dict(feature_id=fid, spot_id=fid, mz='' if mz is None else mz, rt=1,
                             log2fc=1, p_value=.01, q_value=.02, mean_a=1, mean_b=sum(values) / 6,
                             significant='true', adduct=''))
        if peaks:
            msp.append(f'NAME: {fid}\nPRECURSORMZ: {mz}\nIONMODE: Positive\n'
                       + ''.join(f'{p} 100\n' for p in peaks))
            links.append(dict(feature_id=fid, spectrum_file='spectra.msp',
                              spectrum_id=f'spectra:{len(msp)}',
                              source_analysis_id='synthetic-run', mapping_evidence='fixture'))
    (root / 'spectra.msp').write_text('\n'.join(msp), encoding='utf-8')
    tables = [
        ('features', list(features[0]), features),
        ('annotations', ['annotation_id', 'feature_id', 'source', 'name', 'inchikey', 'rank',
                         'provenance_json', 'curation_json', 'status'], []),
        ('spectrum_links', list(links[0]), links),
    ]
    if version >= 3 and with_samples:
        samples = ([dict(sample_id=s, group='b', role='sample', injection_order=i + 1)
                    for i, s in enumerate(TEST_SAMPLES)]
                   + [dict(sample_id=s, group='a', role='sample', injection_order=i + 7)
                      for i, s in enumerate(REFERENCE_SAMPLES)]
                   + [dict(sample_id='QC1', group='', role='qc', injection_order=13)])
        intensities = []
        for fid, _mz, _peaks, (n_test, n_ref), values in _FEATURES:
            for i, sid in enumerate(TEST_SAMPLES):
                detected = i < n_test
                intensities.append(dict(feature_id=fid, sample_id=sid, value=values[i],
                                        detected='' if detection == 'unknown' else str(detected).lower()))
            for i, sid in enumerate(REFERENCE_SAMPLES):
                detected = i < n_ref
                intensities.append(dict(feature_id=fid, sample_id=sid, value=9 if detected else 1,
                                        detected='' if detection == 'unknown' else str(detected).lower()))
        tables += [('samples', ['sample_id', 'group', 'role', 'injection_order'], samples),
                   ('intensities', ['feature_id', 'sample_id', 'value', 'detected'], intensities)]
    files = [dict(role='spectra', path='spectra.msp',
                  sha256=hashlib.sha256((root / 'spectra.msp').read_bytes()).hexdigest())]
    for role, columns, rows in tables:
        path = root / f'{role}.tsv'
        _write_tsv(path, columns, rows)
        files.append(dict(role=role, path=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    comparison = dict(group_a='Vehicle', group_b='Drug', n_a=6, n_b=6,
                      log2fc_sign='positive means group_b is higher', method='Welch', correction='BH',
                      q_threshold=.05, log2fc_threshold=1, preprocess={})
    population = dict(n_features_total=len(features), population_complete=True,
                      definition='全比較対象', exclusions=[])
    manifest = dict(contract_name='massbank-context-evidence', contract_version=version,
                    source_analysis_id='synthetic-run', exporter={'name': 'fixture', 'version': '1'},
                    exported_at='2026-11-04T00:00:00Z', comparison=comparison, population=population,
                    files=files)
    if version >= 3:
        if polarity is not None:
            comparison['polarity'] = polarity
        if with_samples:
            manifest['intensities_stage'] = 'pre_imputation'
        population['detection_unavailable_reason'] = 'sample_axis_mismatch' if detection == 'unknown' else None
    path = root / 'manifest.json'
    path.write_text(json.dumps(manifest, ensure_ascii=False), encoding='utf-8')
    return path
