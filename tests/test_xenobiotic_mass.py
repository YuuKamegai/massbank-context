"""xenobiotic/tables.py と mass.py の単体テスト。"""
from types import SimpleNamespace

import pytest

from massbank_context.core.evidence_models import Feature
from massbank_context.core.formula import parse_formula
from massbank_context.xenobiotic.mass import (
    find_mass_relations, offset_ppm, tolerance_ppm, validate_ppm_bands, validate_ppm_offset)
from massbank_context.xenobiotic.tables import ion_forms


def _feature(fid, mz, adduct=None):
    return Feature(fid, fid, mz, 1.0, None, None, None, None, None, False, adduct)


def _mz(forms, relation):
    return next(f.mz for f in forms if f.relation == relation)


def test_metformin_positive_ions():
    forms = ion_forms(parse_formula('C4H11N5'), 'positive', 0)
    assert _mz(forms, 'parent|[M+H]+') == pytest.approx(130.108722, abs=1e-5)
    assert _mz(forms, 'parent|[M+Na]+') == pytest.approx(152.090666, abs=1e-5)
    assert _mz(forms, 'parent|[2M+H]+') == pytest.approx(259.210167, abs=1e-5)
    assert _mz(forms, 'parent|[M+2H]2+') == pytest.approx(65.557999, abs=1e-5)
    assert _mz(forms, 'parent|[M+H-NH3]+') == pytest.approx(113.082173, abs=1e-5)


def test_metformin_negative_ions():
    forms = ion_forms(parse_formula('C4H11N5'), 'negative', 0)
    assert _mz(forms, 'parent|[M-H]-') == pytest.approx(128.094169, abs=1e-5)
    assert _mz(forms, 'parent|[M+Cl]-') == pytest.approx(164.070847, abs=1e-5)
    assert _mz(forms, 'parent|[M+HCOO]-') == pytest.approx(174.099648, abs=1e-5)


def test_two_shift_metabolites_of_ketamine():
    forms = ion_forms(parse_formula('C13H16ClNO'), 'positive', 2)
    assert _mz(forms, 'demethylation,oxidation|[M+H]+') == pytest.approx(240.078583, abs=1e-5)
    assert _mz(forms, 'demethylation,dehydrogenation|[M+H]+') == pytest.approx(222.068018, abs=1e-5)
    assert not any(f.shifts == ('demethylation', 'methylation') for f in forms)


def test_shift_never_yields_negative_counts():
    forms = ion_forms(parse_formula('CH4'), 'positive', 2)          # 塩素なし・炭素 1
    assert not any('dechlorination_hydroxylation' in f.shifts for f in forms)
    assert not any('deethylation' in f.shifts for f in forms)


def test_loss_needs_the_atoms_in_the_composition():
    # メトホルミン C4H11N5 に酸素は無い
    forms = ion_forms(parse_formula('C4H11N5'), 'positive', 0)
    assert not any('-H2O' in f.ion or '-CO2' in f.ion for f in forms)
    assert any(f.ion == '[M+H-NH3]+' for f in forms)


def test_loss_counts_atoms_of_molecules_and_adduct():
    # ケタミンは酸素 1 つなので −CO2 は 2 分子のときだけ作れる
    forms = ion_forms(parse_formula('C13H16ClNO'), 'positive', 0)
    assert any(f.ion == '[M+H-H2O]+' for f in forms)
    assert not any(f.ion == '[M+H-CO2]+' for f in forms)
    assert any(f.ion == '[2M+H-CO2]+' for f in forms)
    # 付加体が持ち込む原子も数える
    negative = ion_forms(parse_formula('C4H11N5'), 'negative', 0)
    assert any(f.ion == '[M+HCOO-H2O]-' for f in negative)
    assert not any(f.ion == '[M-H-H2O]-' for f in negative)


def test_loss_counts_atoms_added_by_the_shift():
    forms = ion_forms(parse_formula('C4H11N5'), 'positive', 1)
    assert any(f.relation == 'dioxidation|[M+H-H2O]+' for f in forms)


def test_ammonium_loss_still_collapses_into_protonated():
    forms = ion_forms(parse_formula('C4H11N5'), 'positive', 0)
    assert not any(f.ion == '[M+NH4-NH3]+' for f in forms)
    assert sum(abs(f.mz - 130.108722) < 1e-4 for f in forms) == 1


def test_identical_ion_forms_collapse():
    forms = ion_forms(parse_formula('C13H16ClNO'), 'positive', 2)
    at_parent = [f for f in forms if abs(f.mz - 238.099318) < 1e-4]
    assert [f.relation for f in at_parent] == ['parent|[M+H]+']      # [M+NH4-NH3]+ や hydrolysis|[M+H-H2O]+ は畳まれる


def test_max_shifts_limits_combinations():
    assert all(len(f.shifts) <= 1 for f in ion_forms(parse_formula('C4H11N5'), 'positive', 1))


def test_find_relations_orders_by_shift_count_and_ppm():
    entries = [SimpleNamespace(entry_id='X1', composition=parse_formula('C13H16ClNO'))]
    features = [_feature('k_parent', 238.0990), _feature('k_hnk', 240.0786), _feature('no_mz', None)]
    relations, counts, _ = find_mass_relations(features, entries, 'positive', 5.0, [], 2)
    assert relations[('k_parent', 'X1')][0].relation == 'parent|[M+H]+'
    assert relations[('k_parent', 'X1')][0].ppm == pytest.approx(-1.34, abs=0.05)
    hnk = [r.relation for r in relations[('k_hnk', 'X1')]]
    assert 'demethylation,oxidation|[M+H]+' in hnk
    assert [len(r.shifts) for r in relations[('k_hnk', 'X1')]] == sorted(len(r.shifts) for r in relations[('k_hnk', 'X1')])
    assert counts['no_mz'] == 1


def test_adduct_column_drops_conflicting_relations():
    features = [_feature('f', 260.0813, adduct='[M+H]+')]           # ケタミンの [M+Na]+ の位置
    entries = [SimpleNamespace(entry_id='X1', composition=parse_formula('C13H16ClNO'))]
    relations, counts, _ = find_mass_relations(features, entries, 'positive', 10.0, [], 0)
    assert ('f', 'X1') not in relations and counts['adduct_conflict_dropped'] == 1


def test_adducts_are_compared_without_the_source_loss():
    # [M+H-H2O]+ の付加体は [M+H]+（spec §4.1）。ケタミンの −H2O の位置
    entries = [SimpleNamespace(entry_id='X1', composition=parse_formula('C13H16ClNO'))]
    kept, counts, _ = find_mass_relations(
        [_feature('f', 220.0888, adduct='[M+H]1+')], entries, 'positive', 10.0, [], 0)
    assert [r.relation for r in kept[('f', 'X1')]] == ['parent|[M+H-H2O]+']
    assert counts['adduct_conflict_dropped'] == 0


def test_feature_adduct_with_a_loss_is_compared_without_it():
    entries = [SimpleNamespace(entry_id='X1', composition=parse_formula('C13H16ClNO'))]
    kept, counts, _ = find_mass_relations(
        [_feature('f', 238.0990, adduct='[M+H-H2O]1+')], entries, 'positive', 10.0, [], 0)
    assert [r.relation for r in kept[('f', 'X1')]] == ['parent|[M+H]+']
    assert counts['adduct_conflict_dropped'] == 0 and counts['adduct_unrecognized'] == 0


def test_dropped_relations_are_counted_per_pair():
    entries = [SimpleNamespace(entry_id=f'X{i}', composition=parse_formula('C13H16ClNO')) for i in range(2)]
    relations, counts, dropped = find_mass_relations(
        [_feature('f', 220.0888, adduct='[M+Na]+')], entries, 'positive', 10.0, [], 0)
    assert not relations
    assert dropped == {('f', 'X0'): 1, ('f', 'X1'): 1} and counts['adduct_conflict_dropped'] == 2


def test_low_score_adduct_is_not_compared_and_is_counted():
    entries = [SimpleNamespace(entry_id='X1', composition=parse_formula('C13H16ClNO'))]
    feature = _feature('f', 238.0990, adduct='[M+Na]1+')
    relations, counts, dropped = find_mass_relations(
        [feature], entries, 'positive', 10.0, [], 0, frozenset({'f'}))
    assert relations[('f', 'X1')][0].ion == '[M+H]+'
    assert counts['adduct_low_score_ignored'] == 1
    assert counts['adduct_conflict_dropped'] == 0 and not dropped


def test_adduct_of_other_features_is_still_compared():
    entries = [SimpleNamespace(entry_id='X1', composition=parse_formula('C13H16ClNO'))]
    relations, counts, _ = find_mass_relations(
        [_feature('f', 238.0990, adduct='[M+Na]1+')], entries, 'positive', 10.0, [], 0,
        frozenset({'g'}))
    assert ('f', 'X1') not in relations and counts['adduct_conflict_dropped'] == 1
    assert counts['adduct_low_score_ignored'] == 0


def test_ms_dial_adduct_notation_is_normalized():
    entries = [SimpleNamespace(entry_id='X1', composition=parse_formula('C13H16ClNO'))]
    dropped, counts, _ = find_mass_relations(
        [_feature('f', 260.0813, adduct='[M+H]1+')], entries, 'positive', 10.0, [], 0)
    assert ('f', 'X1') not in dropped
    assert counts['adduct_conflict_dropped'] == 1 and counts['adduct_unrecognized'] == 0
    kept, counts, _ = find_mass_relations(
        [_feature('f', 260.0813, adduct='[M+Na]1+')], entries, 'positive', 10.0, [], 0)
    assert kept[('f', 'X1')][0].relation == 'parent|[M+Na]+'
    assert counts['adduct_conflict_dropped'] == 0


def test_unrecognized_adduct_is_kept_and_counted_per_feature():
    entries = [SimpleNamespace(entry_id=f'X{i}', composition=parse_formula('C13H16ClNO')) for i in range(2)]
    relations, counts, _ = find_mass_relations(
        [_feature('f', 238.0990, adduct='[M+Foo]+')], entries, 'positive', 5.0, [], 0)
    assert ('f', 'X0') in relations and ('f', 'X1') in relations
    assert counts['adduct_unrecognized'] == 1


def test_ppm_bands_refuse_unknown_keys():
    with pytest.raises(ValueError, match='未知のキー'):
        validate_ppm_bands([{'min_mz': 200, 'max': 300, 'ppm': 5}])


def test_ppm_offset_is_normalized_like_ms_data_parser():
    assert validate_ppm_offset(None) is None
    assert validate_ppm_offset({'bands': [{'max_mz': 300, 'offset': 1}, {'max_mz': 200, 'offset': -8.3}]}) == {
        'default': 0.0, 'bands': [{'max_mz': 200.0, 'offset': -8.3}, {'max_mz': 300.0, 'offset': 1.0}]}


@pytest.mark.parametrize('raw', [
    [], {'default': 1, 'extra': 0}, {'default': True}, {'default': float('nan')},
    {'bands': [{'max_mz': 200}]}, {'bands': [{'max_mz': 200, 'offset': 1, 'ppm': 5}]},
    {'bands': [{'max_mz': None, 'offset': 1}]}, {'bands': ['x']},
])
def test_ppm_offset_refuses_malformed_specs(raw):
    with pytest.raises(ValueError, match='ppm_offset'):
        validate_ppm_offset(raw)


def test_ppm_offset_band_is_chosen_by_theoretical_mz():
    spec = validate_ppm_offset({'default': 2.0, 'bands': [{'max_mz': 200, 'offset': -8.3}]})
    assert offset_ppm(199.9, spec) == -8.3
    assert offset_ppm(200.0, spec) == 2.0      # 境目は次の帯（上流と同じ `<`）
    assert offset_ppm(150.0, None) == 0.0


def test_ppm_offset_moves_the_window_center_but_not_the_reported_ppm():
    # ケタミン [M+H]+ の理論 m/z は 238.09931。238.1012 は +7.9 ppm、238.0990 は −1.3 ppm。
    entries = [SimpleNamespace(entry_id='X1', composition=parse_formula('C13H16ClNO'))]
    features = [_feature('shifted', 238.1012), _feature('centered', 238.0990)]
    plain, _, _ = find_mass_relations(features, entries, 'positive', 5.0, [], 0)
    assert ('shifted', 'X1') not in plain and ('centered', 'X1') in plain
    moved, _, _ = find_mass_relations(features, entries, 'positive', 5.0, [], 0,
                                      ppm_offset=validate_ppm_offset({'default': 8.0}))
    assert ('centered', 'X1') not in moved
    shifted = [r for r in moved[('shifted', 'X1')] if r.ion == '[M+H]+']
    assert shifted[0].ppm == pytest.approx(7.9, abs=0.1)   # 行の ppm は実測 − 理論のまま


def test_ppm_bands_choose_by_theoretical_mz():
    bands = validate_ppm_bands([{'min_mz': None, 'max_mz': 200, 'ppm': 10}, {'min_mz': 200, 'max_mz': None, 'ppm': 5}])
    assert tolerance_ppm(150.0, 3.0, bands) == 10 and tolerance_ppm(200.0, 3.0, bands) == 5
    with pytest.raises(ValueError):
        validate_ppm_bands([{'min_mz': None, 'max_mz': 300, 'ppm': 10}, {'min_mz': 200, 'max_mz': None, 'ppm': 5}])
