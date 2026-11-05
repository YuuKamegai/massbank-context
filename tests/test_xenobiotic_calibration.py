"""較正の結果（spec §8・§10）を固定する。

`calibration_expectations.json` は閾値を動かす前に確定したもので、結果に合わせて変えない。
満たせなかった項目は ``KNOWN_UNMET`` に理由つきで残し、増えても減っても落ちるようにする
（満たせなかったことを隠さず、直ったことも見逃さない）。

**realdata のテストの通過は spec §8 の合格ではない。** 未達を固定して通しているだけで、§8 の判定は
不合格のまま（設計書 §10.1）。
"""
import csv
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

from massbank_context.core import config
from massbank_context.tools.evidence_tools import load_evidence_bundle
from massbank_context.tools.xenobiotic_tools import register_xenobiotics, screen_xenobiotic_features
from massbank_context.xenobiotic import screen as xeno_screen
from massbank_context.xenobiotic.tables import THRESHOLDS, THRESHOLDS_VERSION

EXPECTATIONS = Path(__file__).parent / 'fixtures' / 'xenobiotic' / 'calibration_expectations.json'
SCRIPTS = Path(__file__).parents[1] / 'scripts'
ORDER = {'A': 0, 'B': 1, 'C': 2, 'D': 3, '': 4}
RUNS = json.loads(EXPECTATIONS.read_text(encoding='utf-8'))['runs']

# 設計書 §10 の決定。
DESIGN_THRESHOLDS = {'msms_shared_min': 0.15, 'msms_top_n': 10, 'bias_reference_max': 1,
                     'rho_min': 0.6, 'rho_min_n': 5}
DESIGN_VERSION = 'xeno-thresholds-2026-11-05-r3'

# 設計書 §10 の「満たせなかった点」。(accession, polarity) -> {特徴 ID: 走査で付いた水準}
KNOWN_UNMET: dict[tuple[str, str], dict[str, str]] = {
    # 6-Propyl-1,3,5-triazine-2,4-diamine。質量の一致は戻ったが共有 0.047、検出 2/5 例で D に留まる
    ('ST003841', 'positive'): {'6126': 'D'},
    # 13981 フルオキセチン: [M+Na]+ の記載で質量が落ち、共有 0。
    # 9187 2-Ethyl-2-(p-tolyl)malonamide: 質量の一致は付加体で落ち、共有 0.05、検出 1/8 例
    # 13249 ノルフルオキセチン: 順位 1 だが PRECURSOR_ONLY で A でなく、共有 0 で B にも届かない
    ('ST002901', 'positive'): {'13981': '', '9187': '', '13249': 'C'},
}

# 設計書 §10.0 の、おとりの上限を超えた件数。(accession, polarity) -> {上限の鍵: 件数}
# 陰イオンの C は MS/MS の閾値に依らないので、格子のどの点でも 32 になる。
KNOWN_OVER_LIMIT: dict[tuple[str, str], dict[str, int]] = {
    ('MTBLS9365', 'negative'): {'max_C': 32},
}


def _feature_levels(path: str) -> dict[str, str]:
    best: dict[str, str] = {}
    with open(path, encoding='utf-8', newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            fid, level = row['feature_id'], row['level']
            if fid not in best or ORDER[level] < ORDER[best[fid]]:
                best[fid] = level
    return best


def test_thresholds_match_design_record():
    assert THRESHOLDS == DESIGN_THRESHOLDS
    assert THRESHOLDS_VERSION == DESIGN_VERSION


def _calibrate_script():
    path = SCRIPTS / 'calibrate_xenobiotic_thresholds.py'
    spec = importlib.util.spec_from_file_location('calibrate_xenobiotic_thresholds', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_grid_point_replaces_the_version_too(monkeypatch):
    monkeypatch.setattr(xeno_screen, 'THRESHOLDS', dict(THRESHOLDS))
    script = _calibrate_script()
    monkeypatch.setattr(script, 'THRESHOLDS', xeno_screen.THRESHOLDS)
    monkeypatch.setattr(xeno_screen, 'THRESHOLDS_VERSION', THRESHOLDS_VERSION)
    script.apply_point(0.3, 20)
    assert xeno_screen.THRESHOLDS['msms_shared_min'] == 0.3
    assert xeno_screen.THRESHOLDS_VERSION == f'{THRESHOLDS_VERSION}+grid-0.3-20'
    script.apply_point(THRESHOLDS['msms_shared_min'], THRESHOLDS['msms_top_n'])
    assert xeno_screen.THRESHOLDS_VERSION == THRESHOLDS_VERSION


def _points(row, top_n=10):
    # shared -> (未達の数, おとりの条件)
    return {(shared, top_n): value for shared, value in row.items()}



def test_selection_takes_the_centre_of_the_tied_run():
    script = _calibrate_script()
    row = {0.02: (3, False), 0.05: (3, False), 0.08: (3, True), 0.1: (3, True), 0.12: (3, True),
           0.15: (3, True), 0.2: (5, True)}
    assert script.select_point(_points(row)) == (0.12, 10)     # 中央 0.115 に最も近い


def test_selection_breaks_a_centre_tie_towards_the_stricter_value():
    script = _calibrate_script()
    # 中央 0.1 が格子に無く 0.08 と 0.12 が等距離なら、B を絞る上の値を取る
    row = {0.05: (4, True), 0.08: (3, True), 0.12: (3, True), 0.15: (4, True)}
    assert script.select_point(_points(row)) == (0.12, 10)


def test_selection_takes_the_only_point():
    script = _calibrate_script()
    row = {0.08: (4, True), 0.1: (3, True), 0.12: (4, True)}
    assert script.select_point(_points(row)) == (0.1, 10)


def test_selection_ignores_points_that_fail_the_decoys_and_other_top_n():
    script = _calibrate_script()
    points = {**_points({0.08: (2, False), 0.1: (3, True)}), **_points({0.02: (3, True)}, top_n=5)}
    assert script.select_point(points) == (0.1, 10)


def test_selection_is_none_without_a_point():
    script = _calibrate_script()
    assert script.select_point(_points({0.1: (3, False)})) is None
    assert script.select_point({(0.1, 5): (3, True)}) is None


def test_expectations_script_refuses_to_overwrite_without_force(tmp_path):
    before = EXPECTATIONS.read_bytes()
    done = subprocess.run([sys.executable, str(SCRIPTS / 'build_xenobiotic_calibration_expectations.py')],
                          cwd=tmp_path, capture_output=True, text=True, encoding='utf-8')
    assert done.returncode != 0 and '--force' in done.stderr
    assert 'ModuleNotFoundError' not in done.stderr
    assert EXPECTATIONS.read_bytes() == before


def test_bundles_script_sets_bytecode_flag_before_imports():
    source = (SCRIPTS / 'build_xenobiotic_calibration_bundles.py').read_text(encoding='utf-8')
    assert 0 <= source.index('sys.dont_write_bytecode = True') < source.index('from metabolomix')


@pytest.mark.realdata
@pytest.mark.parametrize('run', RUNS, ids=lambda run: f"{run['accession']}-{run['polarity']}")
def test_calibration_expectations(run, tmp_path, monkeypatch):
    manifest = Path(run['bundle_dir']) / 'manifest.json'
    if not manifest.is_file():
        pytest.skip(f'バンドルが無い: {manifest}')
    if not config.INDEX_PATH.is_file():
        pytest.skip(f'索引が無い: {config.INDEX_PATH}（MASSBANK_CONTEXT_INDEX で指定する）')
    monkeypatch.setattr(config, 'ARTIFACTS_PATH', tmp_path / 'artifacts')

    analysis_id = json.loads(load_evidence_bundle(str(manifest)))['analysis_id']
    json.loads(register_xenobiotics(analysis_id, run['entries']))
    summary = json.loads(screen_xenobiotic_features(
        analysis_id, run['test_group'], run['ppm'], ppm_bands=run.get('ppm_bands'),
        output_path=str(tmp_path / 'screen.tsv')))
    levels = _feature_levels(summary['output_path'])

    unmet = {item['feature_id']: levels.get(item['feature_id'], '') for item in run['expect']
             if ORDER[levels.get(item['feature_id'], '')] > ORDER[item['min_level']]}
    assert unmet == KNOWN_UNMET.get((run['accession'], run['polarity']), {})

    counts = summary['levels']
    over = {key: counts[level] for key, level in (('max_A', 'A'), ('max_B', 'B'), ('max_C', 'C'))
            if key in run and counts[level] > run[key]}
    assert over == KNOWN_OVER_LIMIT.get((run['accession'], run['polarity']), {})
