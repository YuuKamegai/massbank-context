"""外来の分子の登録（検証・索引の情報・保存）。"""
import pytest

from massbank_context.core import config
from massbank_context.xenobiotic.registry import XenobioticError, load_registration, register
from xenobiotic_helpers import registered_analysis

KET = dict(entry_id='X1', role='drug', name='ketamine', inchikey='YQEZLKZALYSWHR-UHFFFAOYSA-N',
           formula='C13ClH16NO', received_by='test')


@pytest.fixture
def analysis(tmp_path, monkeypatch):
    return registered_analysis(tmp_path, monkeypatch)


def _tsv_rows(text):
    lines = text.strip('\n').split('\n')
    header = lines[0].split('\t')
    return {row[0]: dict(zip(header, row)) for row in (line.split('\t') for line in lines[1:])}


def test_register_reports_index_and_normalizes_formula(analysis):
    store, aid = analysis
    result = register(store, aid, [KET, dict(entry_id='X2', role='drug_metabolite', name='HNK',
                                             formula='C12H14ClNO2', parent_entry_id='X1',
                                             received_by='test')])
    rows = _tsv_rows(result['entries'])
    assert rows['X1']['formula'] == 'C13H16ClNO' and rows['X1']['n_records'] == '2'
    assert rows['X1']['polarities'] == 'POSITIVE:2' and rows['X1']['precursor_types'] == '[M+H]+:2'
    assert rows['X1']['warnings'] == ''
    assert rows['X2']['warnings'] == 'no_inchikey'
    entries, registration_id = load_registration(store, aid)
    assert [e.entry_id for e in entries] == ['X1', 'X2'] and registration_id == result['registration_id']
    assert result['n_entries'] == 2


@pytest.mark.parametrize('change', [
    {'entry_id': ''}, {'role': 'poison'}, {'received_by': 'all'}, {'formula': 'C13Xx'},
    {'inchikey': 'not-a-key'}, {'parent_entry_id': 'X9'}, {'parent_entry_id': 'X1'}, {'extra': 1}])
def test_invalid_entries_raise(analysis, change):
    store, aid = analysis
    with pytest.raises(XenobioticError):
        register(store, aid, [dict(KET, **change)])


def test_duplicate_entry_id_and_empty_list_raise(analysis):
    store, aid = analysis
    for entries in ([KET, KET], []):
        with pytest.raises(XenobioticError):
            register(store, aid, entries)


def test_mass_mismatch_warning(analysis):
    store, aid = analysis
    rows = _tsv_rows(register(store, aid, [dict(KET, formula='C13H18ClNO')])['entries'])
    assert 'formula_mismatch_index:C13H16ClNO' in rows['X1']['warnings']
    assert 'mass_mismatch_index' in rows['X1']['warnings']


def test_not_in_index_warning(analysis):
    store, aid = analysis
    rows = _tsv_rows(register(store, aid, [
        dict(KET, entry_id='A', inchikey='AAAAAAAAAAAAAA-BBBBBBBBBB-C')])['entries'])
    assert rows['A']['warnings'] == 'not_in_index'


def test_reregistration_replaces_all(analysis):
    store, aid = analysis
    register(store, aid, [KET, dict(KET, entry_id='X2')])
    result = register(store, aid, [KET])
    entries, registration_id = load_registration(store, aid)
    assert [e.entry_id for e in entries] == ['X1'] and registration_id == result['registration_id']


def test_missing_registration_is_missing_state(analysis):
    store, aid = analysis
    with pytest.raises(XenobioticError) as info:
        load_registration(store, aid)
    assert info.value.code == 'missing_state'


def test_missing_index_raises(analysis, monkeypatch, tmp_path):
    monkeypatch.setattr(config, 'INDEX_PATH', tmp_path / 'none.sqlite')
    store, aid = analysis
    with pytest.raises(ValueError, match='^missing_index'):
        register(store, aid, [KET])
