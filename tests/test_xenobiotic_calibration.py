"""較正の結果（spec §8・§10）を固定する。

`calibration_expectations.json` は閾値を動かす前に確定したもので、結果に合わせて変えない。
満たせなかった項目は ``KNOWN_UNMET`` に理由つきで残し、増えても減っても落ちるようにする
（満たせなかったことを隠さず、直ったことも見逃さない）。
"""
import csv
import json
from pathlib import Path

import pytest

from massbank_context.core import config
from massbank_context.tools.evidence_tools import load_evidence_bundle
from massbank_context.tools.xenobiotic_tools import register_xenobiotics, screen_xenobiotic_features
from massbank_context.xenobiotic.tables import THRESHOLDS, THRESHOLDS_VERSION

EXPECTATIONS = Path(__file__).parent / 'fixtures' / 'xenobiotic' / 'calibration_expectations.json'
ORDER = {'A': 0, 'B': 1, 'C': 2, 'D': 3, '': 4}
RUNS = json.loads(EXPECTATIONS.read_text(encoding='utf-8'))['runs']

# 設計書 §10 の決定。
DESIGN_THRESHOLDS = {'msms_shared_min': 0.15, 'msms_top_n': 10, 'bias_reference_max': 1,
                     'rho_min': 0.6, 'rho_min_n': 5}
DESIGN_VERSION = 'xeno-thresholds-2026-11-05'

# 設計書 §10 の「満たせなかった点」。(accession, polarity) -> {特徴 ID: 走査で付いた水準}
KNOWN_UNMET: dict[tuple[str, str], dict[str, str]] = {
    # 6-Propyl-1,3,5-triazine-2,4-diamine。付加体の矛盾で質量が落ち、共有 0.047、検出 2/5 例
    ('ST003841', 'positive'): {'6126': ''},
    # 13981 フルオキセチン: [M+Na]+ の記載で質量が落ち、共有 0。
    # 9187 2-Ethyl-2-(p-tolyl)malonamide: 質量の一致が無く、共有 0.05、検出 1/8 例
    ('ST002901', 'positive'): {'13981': '', '9187': ''},
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
    for key, level in (('max_A', 'A'), ('max_B', 'B'), ('max_C', 'C')):
        if key in run:
            assert counts[level] <= run[key]
