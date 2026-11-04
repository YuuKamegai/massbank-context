"""外来の分子の MCP ツール（登録・走査・行の取得）の公開層。"""
import json
from pathlib import Path

import xenobiotic_helpers
from massbank_context.core import config
from massbank_context.tools.evidence_tools import load_evidence_bundle
from massbank_context.tools.xenobiotic_tools import (
    get_xenobiotic_screen, register_xenobiotics, screen_xenobiotic_features)
from xenobiotic_helpers import KETAMINE_KEY

KET = dict(entry_id='X1', role='drug', name='ketamine', inchikey=KETAMINE_KEY,
           formula='C13H16ClNO', received_by='test')


def _load_v3(tmp_path, monkeypatch):
    monkeypatch.setattr(config, 'ARTIFACTS_PATH', tmp_path / 'artifacts')
    monkeypatch.setattr(config, 'INDEX_PATH',
                        xenobiotic_helpers.make_index(tmp_path / 'index' / 'massbank.sqlite'))
    manifest = xenobiotic_helpers.make_v3_bundle(tmp_path / 'bundle')
    return json.loads(load_evidence_bundle(str(manifest)))['analysis_id']


def test_tools_round_trip(tmp_path, monkeypatch):
    aid = _load_v3(tmp_path, monkeypatch)
    reg = json.loads(register_xenobiotics(aid, [KET]))
    assert reg['registration_id'].startswith('xenobiotics:')
    summary = json.loads(screen_xenobiotic_features(aid, 'Drug', 5.0))
    assert summary['levels']['A'] == 1 and Path(summary['output_path']).is_file()
    rows = json.loads(get_xenobiotic_screen(aid, level='A'))
    assert rows['rows'].count('\n') == 1


def test_errors_are_machine_readable(tmp_path, monkeypatch):
    aid = _load_v3(tmp_path, monkeypatch)
    error = json.loads(screen_xenobiotic_features(aid, 'Drug', 5.0))['error']
    assert error['code'] == 'missing_state'
    error = json.loads(screen_xenobiotic_features(aid, 'Nobody', 5.0))['error']
    assert error['code'] == 'unknown_test_group'
    # カーソルの不正・古さもツールの経路で機械可読になる
    json.loads(register_xenobiotics(aid, [KET]))
    json.loads(screen_xenobiotic_features(aid, 'Drug', 5.0))
    error = json.loads(get_xenobiotic_screen(aid, cursor='でたらめ'))['error']
    assert error['code'] == 'invalid_cursor'
    page = json.loads(get_xenobiotic_screen(aid, limit=1))
    assert page['next_cursor']
    json.loads(screen_xenobiotic_features(aid, 'Drug', 5.0, max_shifts=1))
    error = json.loads(get_xenobiotic_screen(aid, cursor=page['next_cursor']))['error']
    assert error['code'] == 'stale_cursor'
