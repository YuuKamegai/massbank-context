"""期待の JSON の各 run を、`msms_shared_min` × `msms_top_n` の格子で走らせて結果を書く（spec §8）。

書くのは `--output`（既定は標準出力だけ）と、run の `xenobiotic_calibration/` の `artifacts_r3/` と
`<pos|neg>/screen_r3.tsv`（`--flat-ppm` は `screen_flat_r3.tsv`）だけ。r2 の成果物は上書きしない。
索引は `MASSBANK_CONTEXT_INDEX` のものを読み取り専用で開く。

使い方:
    MASSBANK_CONTEXT_INDEX=<索引> .venv/Scripts/python.exe scripts/calibrate_xenobiotic_thresholds.py \
        --shared 0.3 --top-n 10 --output <JSON>

`--flat-ppm 10` は帯を持つ run（おとり）だけを帯なしの一律の許容幅で走らせる。
`--select <JSON>` は格子の結果から ``select_point`` の規則で閾値を選んで表示する（走査はしない）。
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
from massbank_context.xenobiotic import screen as xeno_screen  # noqa: E402
from massbank_context.xenobiotic.tables import THRESHOLDS, THRESHOLDS_VERSION  # noqa: E402

EXPECTATIONS = REPO / 'tests' / 'fixtures' / 'xenobiotic' / 'calibration_expectations.json'
ORDER = {'A': 0, 'B': 1, 'C': 2, 'D': 3, '': 4}
DEFAULT_POINT = (THRESHOLDS['msms_shared_min'], THRESHOLDS['msms_top_n'])
#: 前回の格子で未達が最少の 3 を取り、前回の決定で使った値。再較正ではこれを保つ（設計書 §10.4）。
KEPT_TOP_N = 10


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
    config.ARTIFACTS_PATH = bundle.parent.parent / 'artifacts_r3' / bundle.parent.name
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
            'expected': {item['label']: levels.get(item['feature_id'], '') for item in run['expect']},
            'recorded_misses': {item['feature_id']: levels.get(item['feature_id'], '')
                                for item in run.get('recorded_misses', [])}}


def apply_point(shared: float, top_n: int) -> None:
    """閾値を差し替え、版も差し替える（格子の値に決定の版を書かない）。"""
    THRESHOLDS.update(msms_shared_min=shared, msms_top_n=top_n)
    version = THRESHOLDS_VERSION
    if (shared, top_n) != DEFAULT_POINT:
        version = f'{THRESHOLDS_VERSION}+grid-{shared}-{top_n}'
    xeno_screen.THRESHOLDS_VERSION = version


def select_point(points: dict[tuple[float, int], tuple[int, bool]]) -> tuple[float, int] | None:
    """閾値を選ぶ。``points`` は (shared, top_n) -> (未達の数, おとりの条件を満たすか)。

    未達が最少でおとりの条件を満たす点のうち、top_n は ``KEPT_TOP_N`` を保つ。shared は、その
    top_n の行の格子で同点の点が連続する区間の中央に最も近い値を取る。上端は既知の陽性 1 点に
    過適合するため。区間が複数なら点の多い方、同数なら上の方。中央から等距離なら上の値
    （水準 B を絞る側）。選べる点が無ければ None。
    """
    feasible = {key: unmet for key, (unmet, passed) in points.items() if passed}
    if not feasible:
        return None
    least = min(feasible.values())
    row = sorted(shared for shared, top_n in points if top_n == KEPT_TOP_N)
    runs: list[list[float]] = []
    previous = False
    for shared in row:
        tied = feasible.get((shared, KEPT_TOP_N)) == least
        if tied and previous:
            runs[-1].append(shared)
        elif tied:
            runs.append([shared])
        previous = tied
    if not runs:
        return None
    chosen = max(reversed(runs), key=len)
    centre = (chosen[0] + chosen[-1]) / 2
    return min(reversed(chosen), key=lambda shared: round(abs(shared - centre), 9)), KEPT_TOP_N


def grid_points(results: list[dict]) -> dict[tuple[float, int], tuple[int, bool]]:
    """格子の結果を (未達の数, おとりの条件を満たすか) に畳む。"""
    points: dict[tuple[float, int], tuple[int, bool]] = {}
    for row in results:
        key = (row['msms_shared_min'], row['msms_top_n'])
        unmet, passed = points.get(key, (0, True))
        limits = [f for f in row['failures'] if 'limit' in f]
        points[key] = (unmet + len(row['failures']) - len(limits), passed and not limits)
    return points


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shared', type=float, nargs='+', default=[THRESHOLDS['msms_shared_min']])
    parser.add_argument('--top-n', type=int, nargs='+', default=[THRESHOLDS['msms_top_n']])
    parser.add_argument('--output')
    parser.add_argument('--flat-ppm', type=float)
    parser.add_argument('--select')
    args = parser.parse_args()
    if args.select:
        points = grid_points(json.loads(Path(args.select).read_text(encoding='utf-8')))
        print(json.dumps({'selected': select_point(points),
                          'points': [[*key, *value] for key, value in sorted(points.items())]}))
        return
    runs = json.loads(EXPECTATIONS.read_text(encoding='utf-8'))['runs']
    tsv_name = 'screen_r3.tsv'
    if args.flat_ppm is not None:
        runs = [{**run, 'ppm': args.flat_ppm, 'ppm_bands': None} for run in runs if run.get('ppm_bands')]
        tsv_name = 'screen_flat_r3.tsv'
    print('index', config.INDEX_PATH, flush=True)
    results = []
    for run in runs:
        analysis_id = prepare(run)
        base = Path(run['bundle_dir']).parent
        for top_n in args.top_n:
            for shared in args.shared:
                apply_point(shared, top_n)
                started = time.time()
                summary, levels = screen(run, analysis_id, base / tsv_name)
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
