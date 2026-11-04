"""走査（水準・行・TSV・要約・行の取得）。"""
import csv
import hashlib
import sqlite3
from dataclasses import replace

import pytest

import xenobiotic_helpers
from massbank_context.core import config
from massbank_context.xenobiotic import msms, screen
from massbank_context.xenobiotic.registry import XenobioticError, register
from massbank_context.xenobiotic.screen import (
    LEVELS, SCREEN_COLUMNS, PairEvidence, assign_level, get_screen_rows, resolve_polarity,
    resolve_test_side, run_screen)
from massbank_context.xenobiotic.tables import THRESHOLDS, THRESHOLDS_VERSION
from xenobiotic_helpers import KETAMINE_KEY, registered_analysis

KET = dict(entry_id='X1', role='drug', name='ketamine', inchikey=KETAMINE_KEY,
           formula='C13H16ClNO', received_by='test')
HNK = dict(entry_id='X2', role='drug_metabolite', name='HNK', formula='C12H14ClNO2',
           parent_entry_id='X1', received_by='test')


def _read_tsv(path):
    with open(path, encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def _default_output(registered, polarity='positive'):
    _, aid = registered
    return config.ARTIFACTS_PATH / 'xenobiotic' / f"{aid.split(':', 1)[1][:16]}_{polarity}.tsv"


def _feature_levels(rows):
    best = {}
    for row in rows:
        current = best.get(row['feature_id'], '')
        if row['level'] and (not current or row['level'] < current):
            current = row['level']
        best[row['feature_id']] = current
    return best


def _row(rows, feature_id, entry_id):
    return next(r for r in rows if r['feature_id'] == feature_id and r['entry_id'] == entry_id)


def _analysis(tmp_path, monkeypatch, entries=(KET, HNK), **kwargs):
    store, aid = registered_analysis(tmp_path, monkeypatch, **kwargs)
    register(store, aid, list(entries))
    return store, aid


def _alter_bundle(monkeypatch, store, change):
    """read_view が返すバンドルだけを差し替える（契約の検証では作れない入力のため）。"""
    original = store.read_view

    def read_view(analysis_id):
        bundle, revision, heads = original(analysis_id)
        return change(bundle), revision, heads

    monkeypatch.setattr(store, 'read_view', read_view)


def _with_features(monkeypatch, *extra, change=None):
    features = [change(f) if change else f for f in xenobiotic_helpers._FEATURES]
    monkeypatch.setattr(xenobiotic_helpers, '_FEATURES', list(extra) + features)


@pytest.fixture
def registered(tmp_path, monkeypatch):
    return _analysis(tmp_path, monkeypatch)


@pytest.fixture
def screened_store(registered):
    run_screen(*registered, test_group='Drug', ppm=5)
    return registered


@pytest.fixture
def screened(registered):
    summary = run_screen(*registered, test_group='Drug', ppm=5)
    return summary, _read_tsv(summary['output_path'])


@pytest.fixture
def analysis_v2(tmp_path, monkeypatch):
    return _analysis(tmp_path, monkeypatch, version=2)


def test_columns_and_levels_are_fixed():
    assert SCREEN_COLUMNS == [
        'feature_id', 'mz', 'rt', 'log2fc', 'q_value', 'significant', 'level', 'entry_id', 'role',
        'relation', 'theoretical_mz', 'ppm', 'msms_shared_score', 'msms_shared_n',
        'msms_shared_formulas', 'forward_rank', 'detected_test', 'detected_reference',
        'rho_parent', 'notes']
    assert LEVELS == ('A', 'B', 'C', 'D')


@pytest.mark.parametrize('pair, received_by, expected', [
    (PairEvidence(True, True, True, 1, True, False), 'test', 'A'),
    (PairEvidence(True, True, True, 2, None, False), 'test', 'B'),
    (PairEvidence(False, True, False, None, True, False), 'test', 'C'),
    (PairEvidence(False, False, True, None, True, False), 'test', 'C'),
    (PairEvidence(False, True, False, None, None, False), 'test', 'D'),
    (PairEvidence(False, True, False, None, True, False), 'both', 'D'),
    (PairEvidence(False, True, False, None, None, True), 'both', 'C'),
    (PairEvidence(False, False, False, None, True, True), 'test', None),
])
def test_assign_level(pair, received_by, expected):
    assert assign_level(pair, received_by) == expected


def test_screen_levels_on_fixture(screened):
    summary, rows = screened
    level = _feature_levels(rows)
    assert level == {'k_parent': 'A', 'k_hnk': 'B', 'k_frag': 'C', 'x_msms': 'C',
                     'mass_only': 'D', 'endo': '', 'no_mz': ''}
    hnk_x2 = _row(rows, 'k_hnk', 'X2')
    assert hnk_x2['level'] == 'C' and float(hnk_x2['rho_parent']) == pytest.approx(1.0)
    assert 'msms_unavailable:no_inchikey' in hnk_x2['notes']
    assert 'parent_feature:k_parent' in hnk_x2['notes']
    assert _row(rows, 'x_msms', 'X1')['relation'] == 'msms_only'
    assert summary['levels'] == {'A': 1, 'B': 1, 'C': 2, 'D': 1, 'none': 2}
    assert summary['versions'] == {'adduct_table': 'xeno-adducts-1', 'shift_table': 'xeno-shifts-1',
                                   'thresholds': THRESHOLDS_VERSION}
    assert 'rows' not in summary and len({r['feature_id'] for r in rows}) == 7


def test_summary_reports_entries_and_unavailable(screened):
    summary, rows = screened
    assert summary['n_features'] == 7 and summary['n_rows'] == len(rows)
    assert summary['significant_levels'] == summary['levels']
    assert summary['per_entry'] == {'X1': {'A': 1, 'B': 1}, 'X2': {'A': 0, 'B': 0}}
    assert summary['unavailable'] == {
        'msms': {'X2': 'no_inchikey'}, 'n_features_without_msms': 4, 'pattern': None,
        'detection_unavailable_reason': None}
    assert summary['counts']['no_mz'] == 1 and summary['counts']['adduct_conflict_dropped'] == 0
    assert summary['thresholds'] == THRESHOLDS
    assert summary['polarity'] == 'positive' and summary['test_group'] == 'Drug'
    assert summary['params']['ppm'] == 5 and summary['params']['max_shifts'] == 2
    assert summary['next_tools'] == ['get_xenobiotic_screen']


def test_feature_without_mz_is_kept(screened):
    _, rows = screened
    assert 'no_mz' in _row(rows, 'no_mz', '')['notes']


def test_unrelated_feature_is_blank_after_level(screened):
    _, rows = screened
    row = _row(rows, 'endo', '')
    assert row['mz'] == '300.2' and row['significant'] == 'true'
    assert all(row[c] == '' for c in SCREEN_COLUMNS[6:-1])


def test_rows_follow_feature_order_then_strength(screened):
    _, rows = screened
    order = list(dict.fromkeys(r['feature_id'] for r in rows))
    assert order == ['k_parent', 'k_hnk', 'k_frag', 'x_msms', 'mass_only', 'endo', 'no_mz']
    hnk = [r['level'] for r in rows if r['feature_id'] == 'k_hnk']
    assert hnk == sorted(hnk) and hnk[0] == 'B'


def test_relation_row_repeats_pair_values_and_rounds(screened):
    _, rows = screened
    row = _row(rows, 'k_parent', 'X1')
    assert row['relation'] == 'parent|[M+H]+' and row['forward_rank'] == '1'
    assert (row['detected_test'], row['detected_reference']) == ('6', '0')
    assert len(row['ppm'].split('.')[1]) <= 4 and len(row['msms_shared_score'].split('.')[1]) <= 4
    assert len(row['theoretical_mz'].split('.')[1]) <= 6
    assert 'C7H6Cl+' in row['msms_shared_formulas']
    assert int(row['msms_shared_n']) >= 3
    pairs = [r for r in rows if r['feature_id'] == 'k_hnk' and r['entry_id'] == 'X1']
    assert len({(r['level'], r['msms_shared_score'], r['detected_test']) for r in pairs}) == 1


def test_test_side_follows_argument(registered):
    summary = run_screen(*registered, test_group='Vehicle', ppm=5)
    rows = _read_tsv(summary['output_path'])
    row = _row(rows, 'k_frag', 'X1')
    assert (row['detected_test'], row['detected_reference']) == ('0', '6')
    assert row['level'] == 'D'


def test_tsv_bytes_match_recorded_hash(registered, tmp_path):
    target = tmp_path / 'deep' / 'er' / 'out.tsv'
    summary = run_screen(*registered, test_group='Drug', ppm=5, output_path=str(target))
    raw = target.read_bytes()
    assert summary['output_path'] == str(target)
    assert summary['output_sha256'] == hashlib.sha256(raw).hexdigest()
    assert not raw.startswith(b'\xef\xbb\xbf') and b'\r' not in raw
    assert raw.split(b'\n')[0].decode().split('\t') == SCREEN_COLUMNS


def test_snapshot_keeps_rows_and_versions(registered):
    store, aid = registered
    summary = run_screen(store, aid, test_group='Drug', ppm=5)
    payload = store.get_snapshot(aid, summary['screen_revision'])
    assert set(payload) == {'params', 'versions', 'thresholds', 'registration_id', 'summary',
                            'columns', 'rows', 'output_path', 'output_sha256'}
    assert payload['columns'] == SCREEN_COLUMNS and len(payload['rows']) == summary['n_rows']
    assert payload['versions'] == summary['versions']


def test_v2_without_polarity_requires_argument(analysis_v2):
    with pytest.raises(XenobioticError) as info:
        run_screen(*analysis_v2, test_group='Drug', ppm=5)
    assert info.value.code == 'polarity_required'


def test_v2_screen_marks_pattern_unavailable(analysis_v2):
    summary = run_screen(*analysis_v2, test_group='Drug', ppm=5, polarity='positive')
    assert summary['unavailable']['pattern'] == 'no_intensities'
    assert summary['levels']['C'] == 0
    rows = _read_tsv(summary['output_path'])
    assert all('pattern_unavailable:no_intensities' in r['notes'] for r in rows)
    assert all(r['detected_test'] == '' and r['rho_parent'] == '' for r in rows)


def test_unknown_detection_is_not_counted_as_undetected(tmp_path, monkeypatch):
    registered = _analysis(tmp_path, monkeypatch, detection='unknown')
    summary = run_screen(*registered, test_group='Drug', ppm=5)
    assert summary['unavailable']['pattern'] == 'detection_unknown'
    assert summary['unavailable']['detection_unavailable_reason'] == 'sample_axis_mismatch'
    assert summary['levels']['C'] == 0
    rows = _read_tsv(summary['output_path'])
    assert all(r['detected_test'] == '' and r['detected_reference'] == '' for r in rows)
    assert all('pattern_unavailable:detection_unknown' in r['notes'] for r in rows)
    assert _row(rows, 'k_hnk', 'X2')['rho_parent'] == ''


@pytest.mark.parametrize('kwargs, code', [
    (dict(test_group='Placebo'), 'unknown_test_group'),
    (dict(test_group='Drug', polarity='negative'), 'polarity_mismatch'),
    (dict(test_group='Drug', polarity='POSITIVE'), 'invalid_xenobiotic_request'),
    (dict(test_group='Drug', max_shifts=3), 'invalid_xenobiotic_request'),
    (dict(test_group='Drug', max_shifts=True), 'invalid_xenobiotic_request'),
    (dict(test_group='Drug', ppm=0), 'invalid_xenobiotic_request'),
    (dict(test_group='Drug', ppm=float('inf')), 'invalid_xenobiotic_request'),
    (dict(test_group='Drug', ppm_bands=[{'ppm': -1}]), 'invalid_xenobiotic_request'),
    (dict(test_group='Drug', ppm_bands=['x']), 'invalid_xenobiotic_request'),
])
def test_screen_argument_errors(registered, kwargs, code):
    with pytest.raises(XenobioticError) as info:
        run_screen(*registered, **dict(dict(ppm=5), **kwargs))
    assert info.value.code == code


def test_arguments_are_checked_before_registration(tmp_path, monkeypatch):
    store, aid = registered_analysis(tmp_path, monkeypatch)
    with pytest.raises(XenobioticError) as info:
        run_screen(store, aid, test_group='Placebo', ppm=5)
    assert info.value.code == 'unknown_test_group'
    with pytest.raises(XenobioticError) as info:
        run_screen(store, aid, test_group='Drug', ppm=5)
    assert info.value.code == 'missing_state'


def test_unknown_analysis_keeps_store_error(registered):
    with pytest.raises(ValueError, match='解析ID'):
        run_screen(registered[0], 'analysis:nope', test_group='Drug', ppm=5)


def test_resolve_polarity_and_side():
    v3 = {'comparison': {'polarity': 'negative', 'group_a': 'P', 'group_b': 'Q'}}
    assert resolve_polarity(v3, None) == 'negative' and resolve_polarity(v3, 'negative') == 'negative'
    assert resolve_polarity({'comparison': {}}, 'positive') == 'positive'
    assert resolve_test_side(v3, 'P') == 'a' and resolve_test_side(v3, 'Q') == 'b'


def test_max_shifts_zero_drops_two_shift_relation(registered):
    run_screen(*registered, test_group='Drug', ppm=5, max_shifts=0)
    rows = _read_tsv(_default_output(registered))
    assert _row(rows, 'k_hnk', 'X1')['relation'] == 'msms_only'


def test_ppm_bands_override_tolerance(registered):
    summary = run_screen(*registered, test_group='Drug', ppm=5,
                         ppm_bands=[{'min_mz': None, 'max_mz': None, 'ppm': 0.001}])
    rows = _read_tsv(summary['output_path'])
    assert _row(rows, 'k_parent', 'X1')['relation'] == 'msms_only'
    assert summary['params']['ppm_bands'] == [{'min_mz': None, 'max_mz': None, 'ppm': 0.001}]


def test_entry_without_reference_spectra_is_unavailable(tmp_path, monkeypatch):
    other = dict(KET, entry_id='X3', inchikey='AAAAAAAAAAAAAA-AAAAAAAAAA-A')
    registered = _analysis(tmp_path, monkeypatch, entries=[other])
    summary = run_screen(*registered, test_group='Drug', ppm=5)
    assert summary['unavailable']['msms'] == {'X3': 'no_reference_spectra'}
    row = _row(_read_tsv(summary['output_path']), 'k_parent', 'X3')
    assert 'msms_unavailable:no_reference_spectra' in row['notes']
    assert row['msms_shared_score'] == '' and row['forward_rank'] == ''


def test_no_msms_is_distinct_from_zero_score(screened):
    _, rows = screened
    row = _row(rows, 'k_frag', 'X1')
    assert 'no_msms' in row['notes'] and row['msms_shared_score'] == ''
    assert 'no_msms' not in _row(rows, 'k_parent', 'X1')['notes']


def test_empty_spectrum_is_reported(registered, monkeypatch):
    store, aid = registered
    _alter_bundle(monkeypatch, store, lambda b: replace(
        b, spectra={k: replace(s, peaks=()) for k, s in b.spectra.items()}))
    summary = run_screen(store, aid, test_group='Drug', ppm=5)
    rows = _read_tsv(summary['output_path'])
    row = _row(rows, 'k_parent', 'X1')
    assert 'empty_spectrum' in row['notes'] and row['msms_shared_score'] == ''
    assert summary['unavailable']['n_features_without_msms'] == 7
    assert summary['counts']['empty_spectrum'] == 3 and summary['counts']['no_msms'] == 4


def test_dangling_spectrum_link_is_treated_as_absent(registered, monkeypatch):
    store, aid = registered
    _alter_bundle(monkeypatch, store, lambda b: replace(b, spectra={}))
    summary = run_screen(store, aid, test_group='Drug', ppm=5)
    rows = _read_tsv(summary['output_path'])
    assert 'no_msms' in _row(rows, 'k_parent', 'X1')['notes']
    assert summary['counts']['no_msms'] == 7


def test_unrecognized_adduct_is_noted(registered, monkeypatch):
    store, aid = registered
    _alter_bundle(monkeypatch, store, lambda b: replace(
        b, features=tuple(replace(f, adduct='[M+Zz]+') if f.feature_id == 'endo' else f
                          for f in b.features)))
    summary = run_screen(store, aid, test_group='Drug', ppm=5)
    rows = _read_tsv(summary['output_path'])
    assert 'adduct_unrecognized' in _row(rows, 'endo', '')['notes']
    assert 'adduct_unrecognized' not in _row(rows, 'k_parent', 'X1')['notes']
    assert summary['counts']['adduct_unrecognized'] == 1


def test_parent_candidate_excludes_the_feature_itself(tmp_path, monkeypatch):
    # 親の水準 A の特徴と同じ特徴にだけ一致する both の代謝物。自分との ρ=1 で C にしない。
    twin = dict(entry_id='X3', role='drug_metabolite', name='isomer', formula='C13H16ClNO',
                parent_entry_id='X1', received_by='both')
    registered = _analysis(tmp_path, monkeypatch, entries=[KET, twin])
    rows = _read_tsv(run_screen(*registered, test_group='Drug', ppm=5)['output_path'])
    row = _row(rows, 'k_parent', 'X3')
    assert row['level'] == 'D' and row['rho_parent'] == ''
    assert 'parent_feature' not in row['notes']


def test_entry_without_parent_compares_shifted_match_with_own_a_feature(tmp_path, monkeypatch):
    _with_features(monkeypatch, change=lambda f: (
        f[:3] + ((6, 6), [1, 2, 3, 4, 5, 6]) if f[0] == 'mass_only' else f))
    both = dict(KET, received_by='both')
    registered = _analysis(tmp_path, monkeypatch, entries=[both])
    rows = _read_tsv(run_screen(*registered, test_group='Drug', ppm=5)['output_path'])
    row = _row(rows, 'mass_only', 'X1')
    assert row['level'] == 'C' and row['relation'] == 'oxidation|[M+H]+'
    assert 'parent_feature:k_parent' in row['notes']
    # 質量差なしの一致は parent_entry_id が無ければ連動を見ない
    assert _row(rows, 'k_frag', 'X1')['rho_parent'] == ''


def test_highest_rho_parent_feature_is_taken(tmp_path, monkeypatch):
    twin = ('k_parent_b', 238.0990, [125.0153, 207.0571, 220.0888, 238.0990], (6, 0),
            [100, 300, 200, 500, 400, 600])
    _with_features(monkeypatch, twin)
    registered = _analysis(tmp_path, monkeypatch)
    rows = _read_tsv(run_screen(*registered, test_group='Drug', ppm=5)['output_path'])
    assert _feature_levels(rows)['k_parent_b'] == 'A'
    row = _row(rows, 'k_hnk', 'X2')
    assert 'parent_feature:k_parent' in row['notes'].split(';')
    assert float(row['rho_parent']) == pytest.approx(1.0)


def test_short_rho_does_not_hide_linked_parent_feature(tmp_path, monkeypatch):
    # 2 試料だけの ρ=1 が先にあっても、n の足りる親の特徴を取る
    twin = ('k_parent_b', 238.0990, [125.0153, 207.0571, 220.0888, 238.0990], (2, 0),
            [100, 200, 300, 400, 500, 600])
    _with_features(monkeypatch, twin)
    registered = _analysis(tmp_path, monkeypatch)
    rows = _read_tsv(run_screen(*registered, test_group='Drug', ppm=5)['output_path'])
    assert 'parent_feature:k_parent' in _row(rows, 'k_hnk', 'X2')['notes'].split(';')


def test_failed_snapshot_leaves_no_tsv(registered, monkeypatch):
    store, aid = registered

    def refuse(*args, **kwargs):
        raise ValueError('stale_revision: 処理中に更新されました')

    monkeypatch.setattr(store, 'put_snapshot', refuse)
    with pytest.raises(ValueError, match='^stale_revision'):
        run_screen(store, aid, test_group='Drug', ppm=5)
    assert list(_default_output(registered).parent.iterdir()) == []


def test_index_is_closed_when_matching_fails(registered, monkeypatch):
    opened = []

    def spy():
        conn, meta = msms.open_index()
        opened.append(conn)
        return conn, meta

    def explode(*args):
        raise RuntimeError('照合の失敗')

    monkeypatch.setattr(screen, 'open_index', spy)
    monkeypatch.setattr(screen, 'forward_rank', explode)
    with pytest.raises(RuntimeError):
        run_screen(*registered, test_group='Drug', ppm=5)
    with pytest.raises(sqlite3.ProgrammingError):
        opened[0].execute('SELECT 1')
    assert not _default_output(registered).exists()


def test_get_rows_pages_and_filters(screened_store):
    store, aid = screened_store
    first = get_screen_rows(store, aid, level='A', limit=1)
    assert first['n_matched'] == 1 and first['rows'].splitlines()[0].split('\t') == SCREEN_COLUMNS
    page = get_screen_rows(store, aid, limit=2)
    assert page['next_cursor'] and get_screen_rows(store, aid, cursor=page['next_cursor'], limit=2)['rows']
    assert get_screen_rows(store, aid, feature_ids=['nope'])['unknown_feature_ids'] == ['nope']


def test_get_rows_walks_every_row_once(screened_store):
    store, aid = screened_store
    seen, cursor = [], None
    while True:
        page = get_screen_rows(store, aid, cursor=cursor, limit=3)
        seen += page['rows'].splitlines()[1:]
        cursor = page['next_cursor']
        if cursor is None:
            break
    assert len(seen) == page['n_matched'] and len(set(seen)) == len(seen)


def test_get_rows_filters_by_feature(screened_store):
    store, aid = screened_store
    page = get_screen_rows(store, aid, feature_ids=['endo', 'nope'])
    assert page['n_matched'] == 1 and page['unknown_feature_ids'] == ['nope']
    assert page['next_cursor'] is None
    assert page['screen_revision'] == store.head_revision(aid, 'xenobiotic_screen')


def test_get_rows_requires_a_screen(registered):
    with pytest.raises(XenobioticError) as info:
        get_screen_rows(*registered)
    assert info.value.code == 'missing_state'


@pytest.mark.parametrize('kwargs', [
    dict(limit=0), dict(limit=201), dict(limit=True), dict(level='E'), dict(level='')])
def test_get_rows_argument_errors(screened_store, kwargs):
    with pytest.raises(XenobioticError) as info:
        get_screen_rows(*screened_store, **kwargs)
    assert info.value.code == 'invalid_xenobiotic_request'


@pytest.mark.parametrize('cursor', ['***', 'e30=', 'WzFd', 'eyJyZXYiOiJ4In0='])
def test_broken_cursor_is_invalid(screened_store, cursor):
    with pytest.raises(ValueError, match='^invalid_cursor'):
        get_screen_rows(*screened_store, cursor=cursor)


def test_cursor_is_stale_after_rescreen(screened_store):
    store, aid = screened_store
    cursor = get_screen_rows(store, aid, limit=1)['next_cursor']
    run_screen(store, aid, test_group='Drug', ppm=10)
    with pytest.raises(ValueError, match='^stale_cursor'):
        get_screen_rows(store, aid, cursor=cursor)


def test_cursor_is_stale_after_reregistration(screened_store):
    store, aid = screened_store
    cursor = get_screen_rows(store, aid, limit=1)['next_cursor']
    register(store, aid, [KET])
    run_screen(store, aid, test_group='Drug', ppm=5)
    with pytest.raises(ValueError, match='^stale_cursor'):
        get_screen_rows(store, aid, cursor=cursor)
