"""外来の分子の走査の較正 r3 のバンドルを、再解析の Console の出力から作り直す（spec §8、§10.9）。

読むだけ: `datasets/reanalysis/<accession>/pipeline/`、ms-data-parser のコード（main）、r2 の `bundle/`。
書くのは `datasets/reanalysis/<accession>/xenobiotic_calibration/<pos|neg>/` の `bundle_r3/` と
`build_record_r3.json` だけ。r2 の `bundle/`・`bundle_export/` は消さない。

`bundle_r3/` は `dataset_export_evidence_bundle` と同じ呼び方の出力そのもので、手で足さない。
ツールは動作中のセッションのデータセットを読むので、同じ関数を run のデータセットで呼ぶ。
キュレーションと `.arf2` / `.dcl` の決め方もツールと同じにする。

主の比較の結果が `legacy_result_without_dataset_hash` で戻せない run は、run の TSV の頭のメタと
結果の JSON（群の試料）の条件で計算し直し、元の結果と全特徴で一致することを確かめる。

書いたあと、r2 の手で足した links と行ごとに比べ、他の表も比べて記録する。

使い方: `C:/Python314/python.exe scripts/build_xenobiotic_calibration_bundles.py [--force]`
"""
from __future__ import annotations

import argparse
import csv
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
from metabolomix.curation import apply as curation_apply  # noqa: E402
from metabolomix.pipeline import session_handoff  # noqa: E402
from metabolomix.tools.handle_tools import _alignment_files  # noqa: E402

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


def _curation(ds) -> dict | None:
    """`dataset_export_evidence_bundle` の既定（apply_curation=True）と同じ判定。"""
    arf2 = next(filter(None, (curation_apply.arf2_for_mztab(p) for p in ds.source_files)), None)
    if arf2 is None:
        return None
    mapped = (getattr(ds, 'feature_qc', None) or {}).get('source') == 'arf'
    return {'state': 'applied' if mapped else 'unmapped', 'flag_set': curation_apply.flags_for_arf2(arf2)}


def _rows(path: Path) -> list[dict]:
    with open(path, encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def _digest(rows: list[dict], skip: tuple[str, ...] = ()) -> str:
    text = json.dumps([{k: v for k, v in row.items() if k not in skip} for row in rows],
                      sort_keys=True, ensure_ascii=False)
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


def compare_with_r2(new: Path, old: Path) -> dict:
    """r2 の手で足したバンドルとの違い。links は source_analysis_id と根拠の文を除いて比べる。"""
    record = {}
    skip = ('source_analysis_id', 'mapping_evidence')
    links_new, links_old = _rows(new / 'spectrum_links.tsv'), _rows(old / 'spectrum_links.tsv')
    key = lambda row: (row['feature_id'], row['spectrum_file'], row['spectrum_id'])  # noqa: E731
    set_new, set_old = {key(r) for r in links_new}, {key(r) for r in links_old}
    record['links'] = {'n_r3': len(links_new), 'n_r2': len(links_old),
                       'only_r3': len(set_new - set_old), 'only_r2': len(set_old - set_new),
                       'same': set_new == set_old and _digest(links_new, skip) == _digest(links_old, skip),
                       'mapping_evidence_r3': sorted({r['mapping_evidence'] for r in links_new})[:3]}
    for name in ('features', 'samples', 'intensities', 'annotations'):
        left, right = _rows(new / f'{name}.tsv'), _rows(old / f'{name}.tsv')
        record[name] = {'n_r3': len(left), 'n_r2': len(right), 'same': _digest(left) == _digest(right)}
    dcl_new, dcl_old = sorted(new.glob('*.dcl')), sorted(old.glob('*.dcl'))
    record['dcl_same'] = (len(dcl_new) == len(dcl_old) == 1 and dcl_new[0].name == dcl_old[0].name
                          and dcl_new[0].read_bytes() == dcl_old[0].read_bytes())
    return record


def build(run: dict, force: bool) -> dict:
    accession, pol = run['accession'], run['pol']
    root = next((DATASETS / accession / 'pipeline' / pol).glob('pipeline_*'))
    out = DATASETS / accession / 'xenobiotic_calibration' / pol
    bundle_dir = out / 'bundle_r3'
    if force:
        shutil.rmtree(bundle_dir, ignore_errors=True)
    elif bundle_dir.exists():
        raise SystemExit(f'{bundle_dir} が既にある（作り直すなら --force）')

    ds = session_handoff.load_run_dataset(root)
    record = {'accession': accession, 'pol': pol, 'pipeline': root.name,
              'source_verification': ds.source_verification}
    if 'restored_result' in run:
        result = ds.results[run['restored_result']]
        record['result_from'] = {'restored': run['restored_result'],
                                 'group_a': result['a'], 'group_b': result['b']}
    else:
        result, record['result_from'] = _reproduce(ds, root, run['legacy_result'])

    curation = _curation(ds)
    arf2, dcl = _alignment_files(ds, None, 'pipeline_path', str(root))  # dataset_handle_open と同じ
    summary = export_dataset_evidence(ds, result, bundle_dir, curation=curation,
                                      polarity=run['polarity'], arf2_path=arf2, dcl_path=dcl)
    manifest = json.loads((bundle_dir / 'manifest.json').read_text(encoding='utf-8'))
    provenance = manifest['provenance']
    record.update(n_features=summary['n_features'], n_samples=summary['n_samples'],
                  n_intensity_rows=summary['n_intensity_rows'],
                  intensities_stage=summary['intensities_stage'],
                  detection_unavailable_reason=manifest['population']['detection_unavailable_reason'],
                  curation_state=(curation or {}).get('state'),
                  arf2=str(arf2), dcl=str(dcl),
                  spectrum_links_status=provenance.get('spectrum_links_status'),
                  spectrum_links_reason=provenance.get('spectrum_links_reason'),
                  msms_link_status_counts=provenance.get('msms_link_status_counts'),
                  id_space_check={k: v for k, v in (provenance.get('id_space_check') or {}).items()
                                  if k != 'examples'},
                  compared_with_r2=compare_with_r2(bundle_dir, out / 'bundle'))
    (out / 'build_record_r3.json').write_text(
        json.dumps(record, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--force', action='store_true', help='既存の bundle_r3 を作り直す')
    args = parser.parse_args()
    csv.field_size_limit(10 ** 8)
    for run in RUNS:
        print(json.dumps(build(run, args.force), ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
