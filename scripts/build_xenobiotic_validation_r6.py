"""外来の分子の走査の独立の検証 r6（設計書 §10.13）。

r5 では親との共起の規則を試す組が無かったので、陽イオン・DDA・非投与の MTBLS5071 pos（ラットの培養細胞）を
おとりの専用に MS-DIAL で処理し（`xenobiotic_calibration/pos/pipeline/`。再解析の `pipeline/` とは別）、1 回だけ確かめる。

1. その処理の出力から、ms-data-parser の main の `export_dataset_evidence` でバンドルを `bundle_r6/` に作る。
2. 組 D1（アムロジピン・フルオキセチン）と組 D2（ST002901 の麻酔薬）を別々に登録して、版 r5 の 1 点だけで走らせる。
3. §10.13 の条件 1〜3 を機械で判定し、規則が試されたか（共起を数えた組の数）を記録する。

読むだけ: `datasets/reanalysis/MTBLS5071/xenobiotic_calibration/pos/pipeline/`、ms-data-parser のコード（main）、期待の JSON。
書くのは `datasets/reanalysis/MTBLS5071/xenobiotic_calibration/` の下の `pos/bundle_r6/`・`pos/build_record_r6.json`・
`pos/screen_r6_d1.tsv`・`pos/screen_r6_d2.tsv`・`artifacts_r6_d1/`・`artifacts_r6_d2/`・`validation_r6.json` だけ。

使い方: `MASSBANK_CONTEXT_INDEX=<索引> .venv/Scripts/python.exe scripts/build_xenobiotic_validation_r6.py`
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

EXPECTED_VERSION = 'xeno-thresholds-2026-11-05-r5'
EXPECTED_CODETECT = 0.5
EXPECTED_POINT = (0.12, 10)
EXPECTED_INDEX_BUILT_AT = '2026-08-30T15:11:59+00:00'
C_LIMIT = 20

RUNS = (
    {'accession': 'MTBLS5071', 'pol': 'pos', 'polarity': 'positive',
     'reference_group': 'alcohol', 'test_group': 'fasudil', 'd2_strict': True},
)
PPM = 10.0  # r3〜r5 のおとりと同じ（MTBLS5071 は再解析の工程 5 の値がまだ無い）
OUTPUT = r3.DATASETS / 'MTBLS5071' / 'xenobiotic_calibration' / 'validation_r6.json'


def registries() -> dict[str, list[dict]]:
    """D1 は期待の JSON の MTBLS9365、D2 は ST002901 の X3〜X7（変えずに写す）。"""
    runs = json.loads((REPO / 'tests' / 'fixtures' / 'xenobiotic' / 'calibration_expectations.json')
                      .read_text(encoding='utf-8'))['runs']
    decoy = next(r for r in runs if r['accession'] == 'MTBLS9365')['entries']
    st = next(r for r in runs if r['accession'] == 'ST002901' and r['polarity'] == 'positive')['entries']
    d2 = [e for e in st if e['entry_id'] in ('X3', 'X4', 'X5', 'X6', 'X7')]
    if [e['name'] for e in decoy] != ['amlodipine', 'fluoxetine'] or len(d2) != 5:
        raise SystemExit('期待の JSON の登録が §10.13 と違う')
    return {'d1': decoy, 'd2': d2}


def build_bundle(run: dict) -> tuple[Path, dict]:
    out = r3.DATASETS / run['accession'] / 'xenobiotic_calibration' / run['pol']
    roots = sorted((out / 'pipeline').glob('pipeline_*'))
    if len(roots) != 1:
        raise SystemExit(f'おとりの処理の pipeline が 1 つでない: {roots}')
    root = roots[0]
    bundle_dir = out / 'bundle_r6'
    if bundle_dir.exists():
        raise SystemExit(f'{bundle_dir} が既にある（r6 は 1 回だけ。消すかどうかは人が決める）')
    ds = r3.session_handoff.load_run_dataset(root)
    record = {'accession': run['accession'], 'pol': run['pol'], 'pipeline': str(root),
              'source_verification': ds.source_verification}
    found = [(rid, res) for rid, res in ds.results.items()
             if res.get('a') == run['reference_group'] and res.get('b') == run['test_group']]
    if len(found) != 1:
        raise SystemExit(f'{run["reference_group"]} 対 {run["test_group"]} の結果が 1 つでない: {[r for r, _ in found]}')
    result_id, result = found[0]
    record['result_from'] = {'restored': result_id, 'group_a': result['a'], 'group_b': result['b']}
    if result['b'] != run['test_group']:
        raise SystemExit(f'{run["accession"]} {run["pol"]}: test 群が {result["b"]}（§10.13 は {run["test_group"]}）')
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
    (out / 'build_record_r6.json').write_text(json.dumps(record, ensure_ascii=False, indent=1) + '\n',
                                              encoding='utf-8')
    if record['spectrum_links_status'] != 'verified_alignment':
        raise SystemExit(f'{bundle_dir}: spectrum_links_status が {record["spectrum_links_status"]}（§10.13 で止める）')
    return bundle_dir, record


def _tsv(path: str) -> list[dict]:
    with open(path, encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def screen(run: dict, bundle_dir: Path, name: str, entries: list[dict]) -> dict:
    config.ARTIFACTS_PATH = bundle_dir.parent.parent / f'artifacts_r6_{name}' / run['pol']
    shutil.rmtree(config.ARTIFACTS_PATH, ignore_errors=True)
    loaded = json.loads(load_evidence_bundle(str(bundle_dir / 'manifest.json')))
    if 'error' in loaded:
        raise SystemExit(f'{bundle_dir}: {loaded["error"]}')
    registered = json.loads(register_xenobiotics(loaded['analysis_id'], entries))
    if 'error' in registered:
        raise SystemExit(f'{bundle_dir}: {registered["error"]}')
    summary = json.loads(screen_xenobiotic_features(
        loaded['analysis_id'], run['test_group'], PPM,
        output_path=str(bundle_dir.parent / f'screen_r6_{name}.tsv')))
    if 'error' in summary:
        raise SystemExit(f'{run["accession"]} {name}: {summary["error"]}')
    if summary['versions']['thresholds'] != EXPECTED_VERSION or summary['index']['built_at'] != EXPECTED_INDEX_BUILT_AT:
        raise SystemExit(f'版か索引が §10.13 と違う: {summary["versions"]} {summary["index"]}')
    rows = _tsv(summary['output_path'])
    strong = [r for r in rows if r['level'] in ('A', 'B')]
    # 質量差のある一致だけで MS/MS の共有があり、共起を数えた組（B にならなかったものを含む）
    codetect = [r for r in rows if any(n.startswith('parent_codetected:') for n in r['notes'].split(';'))]
    blocked = sorted({r['feature_id'] for r in codetect if r['level'] != 'B'})
    return {'summary': summary, 'rows_ab': strong, 'rows_codetect': codetect,
            'codetect_blocked_features': blocked}


def judge(run: dict, name: str, result: dict) -> dict:
    """§10.13 の条件。D1 は全 run で、D2 は ST003841 だけで A・B が 0。C はすべてで 20 以内。"""
    levels = result['summary']['levels']
    failures = []
    strict = name == 'd1' or run['d2_strict']
    if strict:
        for level in ('A', 'B'):
            if levels.get(level, 0):
                failures.append(f'{name.upper()} の {level} が {levels[level]}（条件は 0）')
    if levels.get('C', 0) > C_LIMIT:
        failures.append(f'{name.upper()} の C が {levels["C"]}（上限 {C_LIMIT}）')
    recorded = [] if strict else sorted({r['feature_id'] for r in result['rows_ab']})
    return {'counted': strict, 'failures': failures, 'recorded_ab': recorded,
            'rule_exercised': bool(result['rows_codetect'])}


def main():
    if (THRESHOLDS['msms_shared_min'], THRESHOLDS['msms_top_n']) != EXPECTED_POINT \
            or THRESHOLDS.get('codetect_min_fraction') != EXPECTED_CODETECT \
            or THRESHOLDS_VERSION != EXPECTED_VERSION:
        raise SystemExit(f'閾値が §10.13 と違う: {THRESHOLDS} {THRESHOLDS_VERSION}')
    if OUTPUT.exists():
        raise SystemExit(f'{OUTPUT} が既にある')
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
            verdict = judge(run, name, result)
            row = {'accession': run['accession'], 'polarity': run['polarity'], 'registry': name,
                   'levels': result['summary']['levels'],
                   'significant_levels': result['summary']['significant_levels'],
                   'per_entry': result['summary']['per_entry'], **verdict,
                   'codetect_blocked_features': result['codetect_blocked_features'],
                   'rows_ab': result['rows_ab'], 'rows_codetect': result['rows_codetect'],
                   'summary': result['summary']}
            results.append(row)
            print(json.dumps({k: row[k] for k in ('accession', 'polarity', 'registry', 'levels',
                                                  'significant_levels', 'counted', 'failures',
                                                  'recorded_ab', 'codetect_blocked_features',
                                                  'rule_exercised')},
                             ensure_ascii=False), flush=True)
    OUTPUT.write_text(json.dumps(results, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print('written', OUTPUT, flush=True)


if __name__ == '__main__':
    main()
