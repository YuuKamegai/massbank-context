"""core/formula.py の単体テスト。"""
import pytest

from massbank_context.core.formula import (
    FormulaError,
    hill_formula,
    monoisotopic_mass,
    parse_formula,
    shift_composition,
)


def test_parse_accepts_msdial_order():
    assert parse_formula('C12ClH14NO2') == {'C': 12, 'Cl': 1, 'H': 14, 'N': 1, 'O': 2}


def test_parse_sums_duplicates():
    assert parse_formula('CH3CH3') == {'C': 2, 'H': 6}


def test_monoisotopic_mass_metformin():
    assert monoisotopic_mass(parse_formula('C4H11N5')) == pytest.approx(129.101445, abs=1e-6)


def test_hill_order():
    assert hill_formula(parse_formula('C12ClH14NO2')) == 'C12H14ClNO2'
    assert hill_formula(parse_formula('ClNa')) == 'ClNa'          # 炭素なしはアルファベット順


@pytest.mark.parametrize('bad', ['', 'C4H11X', 'C4(H2)', 'C4H11N5+', 'c4h11'])
def test_parse_rejects(bad):
    with pytest.raises(FormulaError):
        parse_formula(bad)


def test_shift_composition_refuses_negative():
    base = parse_formula('C4H11N5')
    assert shift_composition(base, {'O': 1, 'H': 1}, {'Cl': 1}) is None
    assert shift_composition(base, {}, {'C': 1, 'H': 2}) == {'C': 3, 'H': 9, 'N': 5}


def test_shift_composition_drops_zero():
    assert shift_composition({'C': 1, 'H': 2}, {}, {'H': 2}) == {'C': 1}
