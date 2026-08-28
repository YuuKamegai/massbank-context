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


def test_zero_spread_falls_back():
    """散らばりがちょうど 0 のとき fallback を返す（幅 0 のビンは to_bin が拒否する）。

    MassBank には同一化合物の同一ピーク列を持つレコードが頻出するため、
    この経路は実際に踏まれる。
    """
    records = [
        _record("A1", "IK1", [100.0, 200.0]),
        _record("A2", "IK1", [100.0, 200.0]),  # 完全に同一
    ]
    profile = measure_mass_error(records)
    # 測定自体は成立し、散らばり 0 が記録される
    assert dict(profile["LC-ESI-QTOF"])[95.0] == 0.0
    # しかしビン幅としては採用されず fallback になる
    assert choose_bin_width(profile) == 0.01


def test_percentiles_are_declared():
    assert MASS_ERROR_PERCENTILES == (50.0, 90.0, 95.0, 99.0)


def test_two_fragments_in_one_record_are_not_instrument_error():
    """1 レコード内の別々のフラグメント同士を装置誤差として数えない。

    測っているのは「同じフラグメントを別々に測ったときのばらつき」なので、
    標本には必ず 2 つ以上のレコード由来の m/z が要る。既定の
    cluster_tolerance 0.05 は高分解能で測ろうとしている誤差の 10 倍ほどあり、
    50 mDa 以内に並ぶ別フラグメントは実在する。1 レコードから作られた
    クラスタを許すと、その距離がそのまま「装置誤差」に化ける。

    choose_bin_width は全装置クラスの 95 パーセンタイルの最大値を採るため、
    汚染された 1 クラスが索引全体のビンを太らせ、照合を緩めて IDF を薄める。
    """
    records = [
        _record("A1", "IK1", [100.000, 100.030]),   # 同一レコード内の 2 本
        _record("A2", "IK1", [500.000]),            # 群の条件を満たすための 2 件目
    ]
    assert measure_mass_error(records) == {}
