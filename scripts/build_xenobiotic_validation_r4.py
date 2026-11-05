"""外来の分子の走査の独立の検証 r4（設計書 §10.11）。

1. MTBLS3233・MTBLS6039 の陽陰のバンドルを、ms-data-parser の main の `export_dataset_evidence` の
   出力そのもので `bundle_r4/` に作る（r3 と同じやり方。`build_xenobiotic_calibration_bundles.py` の関数を使う）。
2. 各バンドルで、組 D1（r3 のおとりの登録）と組 D2（ST002901 の麻酔薬の登録）を別々に登録して、
   r3 の閾値の 1 点だけで走らせる。
3. §10.11 の条件 1〜3 を機械で判定し、条件 4（D1 の A・B）の確かめの対象を並べて書く。

読むだけ: `datasets/reanalysis/<accession>/pipeline/`、ms-data-parser のコード（main）、期待の JSON。
書くのは `datasets/reanalysis/<accession>/xenobiotic_calibration/` の下の `<pol>/bundle_r4/`・`<pol>/build_record_r4.json`・
`<pol>/screen_r4_d1.tsv`・`<pol>/screen_r4_d2.tsv`・`artifacts_r4_d1/`・`artifacts_r4_d2/` と、結果の JSON
（`MTBLS3233/xenobiotic_calibration/validation_r4.json`）だけ。

使い方: `MASSBANK_CONTEXT_INDEX=<索引> C:/Python314/python.exe scripts/build_xenobiotic_validation_r4.py`
"""
from __future__ import annotations

import csv
import importlib.util
import json
import shutil
import sys
from collections import Counter
from pathlib import Path

sys.dont_write_bytecode = True

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

_spec = importlib.util.spec_from_file_location(
    'r3_bundles', REPO / 'scripts' / 'build_xenobiotic_calibration_bundles.py')
r3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(r3)  # ms-data-parser を sys.path に入れ、_reproduce・_curation を持つ

from massbank_context.core import config  # noqa: E402
from massbank_context.tools.evidence_tools import load_evidence_bundle  # noqa: E402
from massbank_context.tools.xenobiotic_tools import (  # noqa: E402
    register_xenobiotics, screen_xenobiotic_features)
from massbank_context.xenobiotic.tables import THRESHOLDS, THRESHOLDS_VERSION  # noqa: E402

EXPECTED_VERSION = 'xeno-thresholds-2026-11-05-r3-selected'
EXPECTED_POINT = (0.12, 10)
EXPECTED_INDEX_BUILT_AT = '2026-08-30T15:11:59+00:00'
C_LIMIT = 20

RUNS = (
    {'accession': 'MTBLS3233', 'pol': 'pos', 'polarity': 'positive', 'pipeline': 'v7/pos',
     'restored_result': 'res_7bdad2e1926140a39bc18f8a5cca542c', 'test_group': 'RA'},
    {'accession': 'MTBLS3233', 'pol': 'neg', 'polarity': 'negative', 'pipeline': 'v7/neg',
     'restored_result': 'res_ea83f2e2e05b42c9951d04d77e3607c8', 'test_group': 'RA'},
    {'accession': 'MTBLS6039', 'pol': 'pos', 'polarity': 'positive', 'pipeline': 'pos',
     'legacy_result': 'res_32c0700a49924acfb070d800bb732a4c', 'test_group': 'Prostate_cancer'},
    {'accession': 'MTBLS6039', 'pol': 'neg', 'polarity': 'negative', 'pipeline': 'neg',
     'legacy_result': 'res_576ebe49bd024b718641ace506aeb82f', 'test_group': 'Prostate_cancer'},
)
PPM = 10.0  # 両 run の判断ログ（工程 5）: 全域 ±10 ppm
OUTPUT = r3.DATASETS / 'MTBLS3233' / 'xenobiotic_calibration' / 'validation_r4.json'


def registries() -> dict[str, list[dict]]:
    """D1 は期待の JSON の MTBLS9365、D2 は ST002901 の X3〜X7（変えずに写す）。"""
    runs = json.loads((REPO / 'tests' / 'fixtures' / 'xenobiotic' / 'calibration_expectations.json')
                      .read_text(encoding='utf-8'))['runs']
    decoy = next(r for r in runs if r['accession'] == 'MTBLS9365')['entries']
    st = next(r for r in runs if r['accession'] == 'ST002901' and r['polarity'] == 'positive')['entries']
    d2 = [e for e in st if e['entry_id'] in ('X3', 'X4', 'X5', 'X6', 'X7')]
    if [e['name'] for e in decoy] != ['amlodipine', 'fluoxetine'] or len(d2) != 5:
        raise SystemExit('期待の JSON の登録が §10.11 と違う')
    return {'d1': decoy, 'd2': d2}


def build_bundle(run: dict) -> tuple[Path, dict]:
    root = next((r3.DATASETS / run['accession'] / 'pipeline' / run['pipeline']).glob('pipeline_*'))
    out = r3.DATASETS / run['accession'] / 'xenobiotic_calibration' / run['pol']
    bundle_dir = out / 'bundle_r4'
    if bundle_dir.exists():
        raise SystemExit(f'{bundle_dir} が既にある（r4 は 1 回だけ。消すかどうかは人が決める）')
    ds = r3.session_handoff.load_run_dataset(root)
    record = {'accession': run['accession'], 'pol': run['pol'], 'pipeline': str(root),
              'source_verification': ds.source_verification}
    if 'restored_result' in run:
        result = ds.results[run['restored_result']]
        record['result_from'] = {'restored': run['restored_result'], 'group_a': result['a'],
                                 'group_b': result['b']}
    else:
        result, record['result_from'] = r3._reproduce(ds, root, run['legacy_result'])
    if result['b'] != run['test_group']:
        raise SystemExit(f'{run["accession"]} {run["pol"]}: test 群が {result["b"]}（§10.11 は {run["test_group"]}）')
    curation = r3._curation(ds)
    arf2, dcl = r3._alignment_files(ds, None, 'pipeline_path', str(root))
    summary = r3.export_dataset_evidence(ds, result, bundle_dir, curation=curation,
                                         polarity=run['polarity'], arf2_path=arf2, dcl_path=dcl)
    provenance = json.loads((bundle_dir / 'manifest.json').read_text(encoding='utf-8'))['provenance']
    record.update(n_features=summary['n_features'], n_samples=summary['n_samples'],
                  intensities_stage=summary['intensities_stage'],
                  curation_state=(curation or {}).get('state'), arf2=str(arf2), dcl=str(dcl),
                  spectrum_links_status=provenance.get('spectrum_links_status'),
                  spectrum_links_reason=provenance.get('spectrum_links_reason'),
                  msms_link_status_counts=provenance.get('msms_link_status_counts'),
                  id_space_check={k: v for k, v in (provenance.get('id_space_check') or {}).items()
                                  if k != 'examples'})
    out.mkdir(parents=True, exist_ok=True)
    (out / 'build_record_r4.json').write_text(json.dumps(record, ensure_ascii=False, indent=1) + '\n',
                                              encoding='utf-8')
    if record['spectrum_links_status'] != 'verified_alignment':
        raise SystemExit(f'{bundle_dir}: spectrum_links_status が {record["spectrum_links_status"]}（§10.11 で止める）')
    return bundle_dir, record


def _tsv(path: str) -> list[dict]:
    with open(path, encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def screen(run: dict, bundle_dir: Path, name: str, entries: list[dict]) -> dict:
    config.ARTIFACTS_PATH = bundle_dir.parent.parent / f'artifacts_r4_{name}' / run['pol']
    shutil.rmtree(config.ARTIFACTS_PATH, ignore_errors=True)
    loaded = json.loads(load_evidence_bundle(str(bundle_dir / 'manifest.json')))
    if 'error' in loaded:
        raise SystemExit(f'{bundle_dir}: {loaded["error"]}')
    registered = json.loads(register_xenobiotics(loaded['analysis_id'], entries))
    if 'error' in registered:
        raise SystemExit(f'{bundle_dir}: {registered["error"]}')
    summary = json.loads(screen_xenobiotic_features(
        loaded['analysis_id'], run['test_group'], PPM,
        output_path=str(bundle_dir.parent / f'screen_r4_{name}.tsv')))
    if 'error' in summary:
        raise SystemExit(f'{run["accession"]} {name}: {summary["error"]}')
    if summary['versions']['thresholds'] != EXPECTED_VERSION or summary['index']['built_at'] != EXPECTED_INDEX_BUILT_AT:
        raise SystemExit(f'版か索引が §10.11 と違う: {summary["versions"]} {summary["index"]}')
    rows = _tsv(summary['output_path'])
    strong = [r for r in rows if r['level'] in ('A', 'B')]
    return {'summary': summary, 'rows_ab': strong}


def judge(name: str, result: dict) -> dict:
    levels = result['summary']['levels']
    failures = []
    if name == 'd2':
        for level in ('A', 'B'):
            if levels.get(level, 0):
                failures.append(f'D2 の {level} が {levels[level]}（条件 1 は 0）')
    if levels.get('C', 0) > C_LIMIT:
        failures.append(f'{name.upper()} の C が {levels["C"]}（上限 {C_LIMIT}）')
    review = [] if name == 'd2' else sorted({r['feature_id'] for r in result['rows_ab']})
    return {'failures': failures, 'needs_review': review}


def main():
    if (THRESHOLDS['msms_shared_min'], THRESHOLDS['msms_top_n']) != EXPECTED_POINT \
            or THRESHOLDS_VERSION != EXPECTED_VERSION:
        raise SystemExit(f'閾値が §10.11 と違う: {THRESHOLDS} {THRESHOLDS_VERSION}')
    csv.field_size_limit(10 ** 8)
    print('index', config.INDEX_PATH, flush=True)
    groups = registries()
    results = []
    for run in RUNS:
        bundle_dir, record = build_bundle(run)
        print(json.dumps({'built': run['accession'], 'pol': run['pol'],
                          'links': record['spectrum_links_status'],
                          'n_features': record['n_features']}, ensure_ascii=False), flush=True)
        for name, entries in groups.items():
            result = screen(run, bundle_dir, name, entries)
            verdict = judge(name, result)
            row = {'accession': run['accession'], 'polarity': run['polarity'], 'registry': name,
                   'levels': result['summary']['levels'],
                   'significant_levels': result['summary']['significant_levels'],
                   'per_entry': result['summary']['per_entry'], **verdict,
                   'rows_ab': result['rows_ab'], 'summary': result['summary']}
            results.append(row)
            print(json.dumps({k: row[k] for k in ('accession', 'polarity', 'registry', 'levels',
                                                  'significant_levels', 'failures', 'needs_review')},
                             ensure_ascii=False), flush=True)
    OUTPUT.write_text(json.dumps(results, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print('written', OUTPUT, flush=True)


if __name__ == '__main__':
    main()
