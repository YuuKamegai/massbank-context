"""外来の分子の走査の較正に使う v3 バンドルを、再解析の Console の出力から作り直す（spec §8）。

読むだけ: `datasets/reanalysis/<accession>/pipeline/`、`reanalysis-study/runs/`、ms-data-parser のコード。
書くのは `datasets/reanalysis/<accession>/xenobiotic_calibration/<pos|neg>/` の下だけ。

- `bundle_export/` は `export_dataset_evidence` の出力そのもの（spectrum_links は空）。
- `bundle/` は、その写しに同じアラインメントの `.dcl` と spectrum_links を足したもの。
  ms-data-parser の dataset 経路は `spectrum_links=[]` で書き、MS/MS の筋が使えない。結び方は
  ARF 経路の `alignment_links` と同じで、ID（SMF_ID == AlignmentID == DCL の番号）で結び、
  m/z と RT は検証だけに使う（m/z の近さでは結ばない）。

主の比較の結果が `legacy_result_without_dataset_hash` で戻せない run は、run の TSV の頭のメタと
結果の JSON（群の試料）の条件で計算し直し、元の結果と全特徴で一致することを確かめる。

使い方: `C:/Python314/python.exe scripts/build_xenobiotic_calibration_bundles.py [--force]`
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path

sys.dont_write_bytecode = True  # ms-data-parser の下に .pyc を書かない

MS_DATA_PARSER = Path(r'C:\Users\yuu18\Metabolomix_with_LLM')
sys.path.insert(0, str(MS_DATA_PARSER))

from metabolomix.analysis.dataset_service import compare_dataset, preprocess_dataset  # noqa: E402
from metabolomix.analysis.evidence_export import export_dataset_evidence  # noqa: E402
from metabolomix.arf2.reader import load_catalog  # noqa: E402
from metabolomix.dcl.reader import deserialize_dcl  # noqa: E402
from metabolomix.msdial.alignment_spot import precursor_agrees  # noqa: E402
from metabolomix.pipeline import session_handoff  # noqa: E402

DATASETS = Path(r'C:\Users\yuu18\datasets\reanalysis')

# 主の比較。legacy の run は元の結果の result_id（run の TSV の頭のメタ）を書く。
RUNS = (
    {'accession': 'ST003841', 'pol': 'pos', 'polarity': 'positive',
     'legacy_result': 'res_fdc3a93bb259417c995dc406e5081df9'},
    {'accession': 'ST002901', 'pol': 'pos', 'polarity': 'positive',
     'legacy_result': 'res_dc2bb24c91fb463fad7f4ebaf7929f6c'},
    {'accession': 'ST002901', 'pol': 'neg', 'polarity': 'negative',
     'legacy_result': 'res_46cd94e741d5490489999881a046bcc2'},
    {'accession': 'MTBLS9365', 'pol': 'pos', 'polarity': 'positive',
     'restored_result': 'res_995c9ad76e674040954dfe8c3efdf025'},
    {'accession': 'MTBLS9365', 'pol': 'neg', 'polarity': 'negative',
     'restored_result': 'res_238bd362e6624124a498632c41bc22f0'},
)

# 両 run の TSV の頭のメタ `# preprocess` と同じ（conservative-v1 の既定）
LEGACY_RECIPE = {'normalize': 'none', 'blank_min_fold': None, 'drift_correct': False,
                 'max_qc_rsd': None, 'impute': 'half_min', 'min_detection_rate': 0.0}


def _close(left, right, tolerance=1e-9):
    """保存された結果は NaN を None で持つので、両者を同じ欠測として比べる。"""
    left = None if left is None or math.isnan(left) else left
    right = None if right is None or math.isnan(right) else right
    if left is None or right is None:
        return left is None and right is None
    return left == right or abs(left - right) <= tolerance * max(1.0, abs(left), abs(right))


def _reproduce(ds, root: Path, legacy_id: str) -> tuple[dict, dict]:
    """元の結果の条件で計算し直し、(結果, 一致の記録) を返す。一致しなければ止める。"""
    stored = json.loads((root / 'results' / f'{legacy_id}.json').read_text(encoding='utf-8'))['data']
    preprocess_dataset(ds, dict(LEGACY_RECIPE))
    result = compare_dataset(ds, stored['samples_a'], stored['samples_b'],
                             q_threshold=stored['q_threshold'],
                             log2fc_threshold=stored['log2fc_threshold'],
                             log_transform=stored['log_transform'],
                             group_a_label=stored['a'], group_b_label=stored['b'])
    before = {str(row['feature']): row for row in stored['results']}
    after = {str(row['feature']): row for row in result['results']}
    if set(before) != set(after):
        raise SystemExit(f'{legacy_id}: 特徴の集合が元の結果と違う')
    differing = [fid for fid, row in after.items()
                 if not all(_close(row.get(k), before[fid].get(k))
                            for k in ('log2fc', 'p', 'q', 'mean_a', 'mean_b'))]
    if differing:
        raise SystemExit(f'{legacy_id}: {len(differing)} 特徴の統計が元の結果と違う: {differing[:5]}')
    return result, {'reproduced_from': legacy_id, 'n_features_compared': len(after),
                    'n_differing': 0, 'columns': ['log2fc', 'p', 'q', 'mean_a', 'mean_b'],
                    'group_a': stored['a'], 'group_b': stored['b']}


def _spectrum_links(result: dict, ds, root: Path, source_id: str) -> tuple[Path, list[dict], dict]:
    """ARF 経路の `alignment_links` と同じ規則の結び付け。m/z・RT は検証だけに使う。"""
    arf2 = sorted((root / 'input').glob('AlignResult-*.arf2'), key=lambda p: p.name, reverse=True)[0]
    dcl = arf2.with_suffix('.dcl')
    catalog = {spot['AlignmentID']: spot for spot in load_catalog(arf2)}
    entries = deserialize_dcl(str(dcl))
    links, statuses = [], {'ok': 0, 'no_msms': 0, 'precursor_mismatch': 0,
                           'metadata_missing': 0, 'out_of_range': 0}
    id_space_mismatch = 0
    for row in result['results']:
        fid = str(row['feature'])
        spot = int(fid)
        raw = catalog.get(spot, {})
        meta = ds.feature_metadata.get(row['feature'], {})
        if not 0 <= spot < len(entries):
            statuses['out_of_range'] += 1
            continue
        # SMF_ID が AlignmentID と同じ空間であることの確かめ（mzTab の値と .arf2 の値）
        if raw.get('MassCenter') is not None and meta.get('mz') is not None \
                and not precursor_agrees(meta['mz'], raw['MassCenter'], meta.get('rt'), raw.get('RT')):
            id_space_mismatch += 1
        entry = entries[spot]
        if raw.get('MassCenter') is None or raw.get('RT') is None:
            status = 'metadata_missing'
        elif not precursor_agrees(entry['precursor_mz'], raw['MassCenter'], entry['rt'], raw['RT']):
            status = 'precursor_mismatch'
        elif not entry['msms_spectrum']:
            status = 'no_msms'
        else:
            status = 'ok'
            links.append({'feature_id': fid, 'spectrum_file': dcl.name,
                          'spectrum_id': f'{dcl.stem}:{spot}', 'source_analysis_id': source_id,
                          'mapping_evidence': 'same-alignment-stem; SMF_ID == AlignmentID == DCL index; '
                                              'precursor/RT verified'})
        statuses[status] += 1
    if id_space_mismatch:
        raise SystemExit(f'{root.name}: SMF_ID と AlignmentID が {id_space_mismatch} 件食い違う')
    return dcl, links, statuses


def _add_links(export_dir: Path, bundle_dir: Path, dcl: Path, links: list[dict], statuses: dict):
    """`bundle_export` の写しに `.dcl` と spectrum_links を足す。他の表は変えない。"""
    shutil.copytree(export_dir, bundle_dir)
    table = bundle_dir / 'spectrum_links.tsv'
    columns = table.read_text(encoding='utf-8').splitlines()[0].split('\t')
    lines = ['\t'.join(columns)] + ['\t'.join(link[c] for c in columns) for link in links]
    table.write_bytes(('\n'.join(lines) + '\n').encode('utf-8'))
    shutil.copyfile(dcl, bundle_dir / dcl.name)

    manifest_path = bundle_dir / 'manifest.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    for item in manifest['files']:
        if item['role'] == 'spectrum_links':
            item['sha256'] = hashlib.sha256(table.read_bytes()).hexdigest()
    manifest['files'].append({'role': 'spectra', 'path': dcl.name,
                              'sha256': hashlib.sha256((bundle_dir / dcl.name).read_bytes()).hexdigest()})
    manifest['provenance']['spectrum_links_status'] = 'verified_alignment'
    manifest['provenance']['msms_link_status_counts'] = statuses
    manifest['provenance']['spectrum_links_added_by'] = (
        'massbank-context scripts/build_xenobiotic_calibration_bundles.py')
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, sort_keys=True,
                                        separators=(',', ':')) + '\n', encoding='utf-8')


def build(run: dict, force: bool) -> dict:
    accession, pol = run['accession'], run['pol']
    root = next((DATASETS / accession / 'pipeline' / pol).glob('pipeline_*'))
    out = DATASETS / accession / 'xenobiotic_calibration' / pol
    export_dir, bundle_dir = out / 'bundle_export', out / 'bundle'
    if force:
        for directory in (export_dir, bundle_dir):
            shutil.rmtree(directory, ignore_errors=True)
    elif export_dir.exists() or bundle_dir.exists():
        raise SystemExit(f'{out} に既にバンドルがある（作り直すなら --force）')
    out.mkdir(parents=True, exist_ok=True)

    ds = session_handoff.load_run_dataset(root)
    record = {'accession': accession, 'pol': pol, 'pipeline': root.name,
              'source_verification': ds.source_verification}
    if 'restored_result' in run:
        result = ds.results[run['restored_result']]
        record['result_from'] = {'restored': run['restored_result'],
                                 'group_a': result['a'], 'group_b': result['b']}
    else:
        result, record['result_from'] = _reproduce(ds, root, run['legacy_result'])

    summary = export_dataset_evidence(ds, result, export_dir, polarity=run['polarity'])
    manifest = json.loads((export_dir / 'manifest.json').read_text(encoding='utf-8'))
    dcl, links, statuses = _spectrum_links(result, ds, root, manifest['source_analysis_id'])
    _add_links(export_dir, bundle_dir, dcl, links, statuses)
    record.update(n_features=summary['n_features'], n_samples=summary['n_samples'],
                  n_intensity_rows=summary['n_intensity_rows'],
                  intensities_stage=summary['intensities_stage'],
                  detection_unavailable_reason=manifest['population']['detection_unavailable_reason'],
                  n_links=len(links), link_statuses=statuses, dcl=dcl.name)
    (out / 'build_record.json').write_text(
        json.dumps(record, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--force', action='store_true', help='既存のバンドルを作り直す')
    args = parser.parse_args()
    for run in RUNS:
        print(json.dumps(build(run, args.force), ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
