"""期待の JSON の各 run を、`msms_shared_min` × `msms_top_n` の格子で走らせて結果を書く（spec §8）。

書くのは `--output`（既定は標準出力だけ）と、run の `xenobiotic_calibration/artifacts/` の下だけ。
索引は `MASSBANK_CONTEXT_INDEX` のものを読み取り専用で開く。

使い方:
    MASSBANK_CONTEXT_INDEX=<索引> C:/Python314/python.exe scripts/calibrate_xenobiotic_thresholds.py \
        --shared 0.3 --top-n 10 --output <JSON>
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
import time
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from massbank_context.core import config  # noqa: E402
from massbank_context.tools.evidence_tools import load_evidence_bundle  # noqa: E402
from massbank_context.tools.xenobiotic_tools import (  # noqa: E402
    register_xenobiotics, screen_xenobiotic_features)
from massbank_context.xenobiotic.tables import THRESHOLDS  # noqa: E402

EXPECTATIONS = REPO / 'tests' / 'fixtures' / 'xenobiotic' / 'calibration_expectations.json'
ORDER = {'A': 0, 'B': 1, 'C': 2, 'D': 3, '': 4}


def feature_levels(tsv_path: str) -> dict[str, str]:
    """特徴ごとの最も強い水準（水準の無い特徴は空）。"""
    best: dict[str, str] = {}
    with open(tsv_path, encoding='utf-8', newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            fid, level = row['feature_id'], row['level']
            if fid not in best or ORDER[level] < ORDER[best[fid]]:
                best[fid] = level
    return best


def prepare(run: dict) -> str:
    """artifacts を空から作り直して読み込み・登録する。analysis_id を返す。"""
    bundle = Path(run['bundle_dir'])
    config.ARTIFACTS_PATH = bundle.parent.parent / 'artifacts' / bundle.parent.name
    shutil.rmtree(config.ARTIFACTS_PATH, ignore_errors=True)
    loaded = json.loads(load_evidence_bundle(str(bundle / 'manifest.json')))
    if 'error' in loaded:
        raise SystemExit(f'{bundle}: {loaded["error"]}')
    registered = json.loads(register_xenobiotics(loaded['analysis_id'], run['entries']))
    if 'error' in registered:
        raise SystemExit(f'{bundle}: {registered["error"]}')
    return loaded['analysis_id']


def screen(run: dict, analysis_id: str, output_path: Path) -> tuple[dict, dict[str, str]]:
    summary = json.loads(screen_xenobiotic_features(
        analysis_id, run['test_group'], run['ppm'], ppm_bands=run.get('ppm_bands'),
        output_path=str(output_path)))
    if 'error' in summary:
        raise SystemExit(f'{run["accession"]}: {summary["error"]}')
    return summary, feature_levels(summary['output_path'])


def judge(run: dict, levels: dict[str, str]) -> dict:
    """期待との突き合わせ。満たさなかった項目を ``failures`` に並べる。"""
    failures = [{'feature_id': item['feature_id'], 'label': item['label'],
                 'min_level': item['min_level'], 'level': levels.get(item['feature_id'], '')}
                for item in run['expect']
                if ORDER[levels.get(item['feature_id'], '')] > ORDER[item['min_level']]]
    counts = Counter(level for level in levels.values() if level)
    for key, level in (('max_A', 'A'), ('max_B', 'B'), ('max_C', 'C')):
        if key in run and counts[level] > run[key]:
            failures.append({'limit': key, 'allowed': run[key], 'count': counts[level]})
    return {'counts': dict(counts), 'failures': failures,
            'expected': {item['label']: levels.get(item['feature_id'], '') for item in run['expect']}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shared', type=float, nargs='+', default=[THRESHOLDS['msms_shared_min']])
    parser.add_argument('--top-n', type=int, nargs='+', default=[THRESHOLDS['msms_top_n']])
    parser.add_argument('--output')
    args = parser.parse_args()
    runs = json.loads(EXPECTATIONS.read_text(encoding='utf-8'))['runs']
    print('index', config.INDEX_PATH, flush=True)
    results = []
    for run in runs:
        analysis_id = prepare(run)
        base = Path(run['bundle_dir']).parent
        for top_n in args.top_n:
            for shared in args.shared:
                THRESHOLDS.update(msms_shared_min=shared, msms_top_n=top_n)
                started = time.time()
                summary, levels = screen(run, analysis_id, base / 'screen.tsv')
                verdict = judge(run, levels)
                row = {'accession': run['accession'], 'polarity': run['polarity'],
                       'msms_shared_min': shared, 'msms_top_n': top_n, **verdict,
                       'seconds': round(time.time() - started, 1)}
                results.append({**row, 'summary': summary})
                print(json.dumps(row, ensure_ascii=False), flush=True)
    if args.output:
        Path(args.output).write_text(json.dumps(results, ensure_ascii=False, indent=1) + '\n',
                                     encoding='utf-8')


if __name__ == '__main__':
    main()
