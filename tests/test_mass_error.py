from massbank_context.core.models import Peak, SpectrumRecord
from massbank_context.index.mass_error import (
    MASS_ERROR_PERCENTILES,
    choose_bin_width,
    measure_mass_error,
)


def _record(accession: str, inchikey: str, mzs: list[float],
            instrument_type: str = "LC-ESI-QTOF") -> SpectrumRecord:
    return SpectrumRecord(
        accession=accession, inchikey=inchikey, precursor_mz=300.0,
        precursor_type="[M+H]+", polarity="POSITIVE",
        instrument_type=instrument_type, ms_type="MS2", resolution_value=None,
        collision_energy=None, splash=None, license=None,
        molecular_formula=None, exact_mass=None,
        peaks=tuple(Peak(mz, 100.0, None, None) for mz in mzs),
    )


def test_measures_dispersion_within_same_inchikey():
    """同じ InChIKey の 2 レコードで同じフラグメントが 0.004 ずれている場合。"""
    records = [
        _record("A1", "IK1", [100.000, 200.000]),
        _record("A2", "IK1", [100.004, 200.004]),
    ]
    profile = measure_mass_error(records)
    assert "LC-ESI-QTOF" in profile
    deltas = dict(profile["LC-ESI-QTOF"])
    assert abs(deltas[95.0] - 0.004) < 1e-6


def test_ignores_singleton_inchikeys():
    """レコードが 1 件しかない InChIKey からは散らばりを測れない。"""
    profile = measure_mass_error([_record("A1", "IK1", [100.0])])
    assert profile == {}


def test_separates_instrument_classes():
    records = [
        _record("A1", "IK1", [100.000], "LC-ESI-QTOF"),
        _record("A2", "IK1", [100.002], "LC-ESI-QTOF"),
        _record("B1", "IK2", [100.000], "LC-ESI-QQ"),
        _record("B2", "IK2", [100.030], "LC-ESI-QQ"),
    ]
    profile = measure_mass_error(records)
    assert set(profile) == {"LC-ESI-QTOF", "LC-ESI-QQ"}
    assert dict(profile["LC-ESI-QTOF"])[95.0] < dict(profile["LC-ESI-QQ"])[95.0]


def test_distant_peaks_are_not_clustered():
    """0.05 Da 以上離れたピークは同一フラグメントとみなさない。"""
    records = [
        _record("A1", "IK1", [100.000]),
        _record("A2", "IK1", [180.000]),
    ]
    assert measure_mass_error(records) == {}


def test_choose_bin_width_uses_max_across_classes():
    """索引は 1 つのビン幅で作るので、最も散らばる装置に合わせる。"""
    profile = {
        "LC-ESI-QTOF": [(95.0, 0.004)],
        "LC-ESI-ITFT": [(95.0, 0.002)],
    }
    assert choose_bin_width(profile) == 0.004


def test_choose_bin_width_falls_back_when_empty():
    assert choose_bin_width({}) == 0.01


def test_percentiles_are_declared():
    assert MASS_ERROR_PERCENTILES == (50.0, 90.0, 95.0, 99.0)
