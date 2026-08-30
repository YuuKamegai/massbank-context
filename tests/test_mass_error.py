from massbank_context.core.models import Peak, SpectrumRecord
from massbank_context.index.mass_error import (
    MASS_ERROR_PERCENTILES,
    InstrumentErrorProfile,
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


def test_measures_deviation_between_records_of_the_same_compound():
    """標本は「同じフラグメントを別レコードで測った差」である。"""
    records = [
        _record("A1", "IK1", [100.000, 200.000]),
        _record("A2", "IK1", [100.004, 200.000]),
    ]
    profile = measure_mass_error(records)["LC-ESI-QTOF"]
    assert profile.n_records == 2
    assert profile.n_samples == 2
    # 値はヒストグラムの刻みに量子化され、そのビンの上端を返す
    assert 0.004 <= profile.value_at(99.0) <= 0.0045


def test_ignores_singleton_inchikeys():
    """1 レコードしか無い化合物からは散らばりを測れない。"""
    assert measure_mass_error([_record("A1", "IK1", [100.0])]) == {}


def test_separates_instrument_classes():
    records = [
        _record("A1", "IK1", [100.000], instrument_type="LC-ESI-QTOF"),
        _record("A2", "IK1", [100.004], instrument_type="LC-ESI-QTOF"),
        _record("B1", "IK1", [100.000], instrument_type="LC-ESI-ITFT"),
        _record("B2", "IK1", [100.001], instrument_type="LC-ESI-ITFT"),
    ]
    profile = measure_mass_error(records)
    assert set(profile) == {"LC-ESI-QTOF", "LC-ESI-ITFT"}
    assert (profile["LC-ESI-ITFT"].value_at(99.0)
            < profile["LC-ESI-QTOF"].value_at(99.0))


def test_two_fragments_in_one_record_are_not_instrument_error():
    """1 レコード内の別々のフラグメント同士を装置誤差として数えない。

    測っているのは「同じフラグメントを別々に測ったときのばらつき」なので、
    標本は必ずレコードを跨いだ対でなければならない。
    """
    records = [
        _record("A1", "IK1", [100.000, 100.030]),   # 同一レコード内の 2 本
        _record("A2", "IK1", [500.000]),            # 群の条件を満たすための 2 件目
    ]
    assert measure_mass_error(records) == {}


def test_orphan_peak_does_not_steal_an_already_paired_peak():
    """相手を持たないピークが、既に対になっているピークに寄生しない。

    最近傍を片方向で採ると、100.030 は相手が居ないまま 100.002 を掴み、
    28 mDa が「装置誤差」として標本に入る。対応は相互最近傍で決める。
    """
    records = [
        _record("A1", "IK1", [100.000, 100.030]),
        _record("A2", "IK1", [100.002]),
    ]
    profile = measure_mass_error(records)["LC-ESI-QTOF"]
    assert profile.n_samples == 1
    assert profile.value_at(99.0) <= 0.0025


def test_decoy_fragments_do_not_inflate_the_estimate():
    """対応しないフラグメントの偶然の近接で推定値が太らない。

    これが実データで踏んだ症状の再現である。旧実装は標本を集める窓
    （cluster_tolerance）を変えると出力がその 85〜90% に追随した。つまり
    測っていたのは装置誤差ではなく窓の幅そのものだった。窓の中に散らばる
    無関係なフラグメントは一様な背景として差し引く。
    """
    records = []
    for i in range(300):
        base = 100.0 + i
        records.append(_record(f"A{i}", f"IK{i}", [base, base + 200.0]))
        # 1 本は真の対応（偏差 2 mDa）、もう 1 本は対応しない別フラグメント
        records.append(_record(f"B{i}", f"IK{i}", [base + 0.002, base + 200.030]))
    profile = measure_mass_error(records)["LC-ESI-QTOF"]
    assert profile.n_samples == 600
    assert profile.value_at(95.0) <= 0.003


def test_reports_record_and_sample_counts():
    """標本数を値と一緒に運ぶ（choose_bin_width が採否を決めるのに要る）。"""
    records = [
        _record("A1", "IK1", [100.000, 200.000]),
        _record("A2", "IK1", [100.004, 200.001]),
        _record("A3", "IK1", [100.002, 200.000]),
    ]
    profile = measure_mass_error(records)["LC-ESI-QTOF"]
    assert profile.instrument_class == "LC-ESI-QTOF"
    assert profile.n_records == 3
    # 隣接ペア 2 組 × 2 本
    assert profile.n_samples == 4


def test_zero_spread_falls_back():
    """散らばりがちょうど 0 のとき fallback を返す（幅 0 のビンは to_bin が拒否する）。

    MassBank には同一化合物の同一ピーク列を持つレコード（同じスペクトルの
    重複登録）が頻出するため、この経路は実際に踏まれる。重複からは装置誤差を
    測れないので、量子化して 0 でない値を返してはならない。
    """
    records = [
        _record("A1", "IK1", [100.0, 200.0]),
        _record("A2", "IK1", [100.0, 200.0]),  # 完全に同一
    ]
    profile = measure_mass_error(records)
    assert profile["LC-ESI-QTOF"].value_at(95.0) == 0.0
    assert choose_bin_width(profile, min_samples=1) == 0.01


def test_percentiles_are_declared():
    assert MASS_ERROR_PERCENTILES == (50.0, 90.0, 95.0, 99.0)


def _profile(instrument_class: str, value: float, n_records: int,
             n_samples: int = 10_000) -> InstrumentErrorProfile:
    return InstrumentErrorProfile(
        instrument_class=instrument_class, n_records=n_records,
        n_samples=n_samples, percentiles=((95.0, value),))


def test_choose_bin_width_ignores_classes_with_too_few_samples():
    """標本の乏しい装置クラスに索引全体のビン幅を決めさせない。

    実データで踏んだ経路: 45 レコード（索引の 0.05%）の MALDI-TOFTOF の p95 が
    84,469 レコード分の照合窓を決めていた。標本が少なければ p95 は 1〜2 本の
    外れ値そのものであり、装置誤差の推定になっていない。
    """
    profiles = {
        "LC-ESI-QTOF": _profile("LC-ESI-QTOF", 0.004, n_records=46052),
        "MALDI-TOFTOF": _profile("MALDI-TOFTOF", 0.048, n_records=45,
                                 n_samples=12),
    }
    assert choose_bin_width(profiles) == 0.004


def test_choose_bin_width_covers_records_not_the_widest_class():
    """レコード数で重み付けした分位点を採る（クラス横断の最大値ではない）。

    標本数が足りていても、索引の 1% 未満しか無いクラスの散らばりに全体を
    合わせると、残り 99% の照合が不必要に緩くなる。ビン幅は「索引の
    レコードの大半を覆えるか」で決める。
    """
    profiles = {
        "LC-ESI-ITFT": _profile("LC-ESI-ITFT", 0.002, n_records=40_000),
        "LC-ESI-QTOF": _profile("LC-ESI-QTOF", 0.010, n_records=46_000),
        "MALDI-TOFTOF": _profile("MALDI-TOFTOF", 0.048, n_records=800),
    }
    assert choose_bin_width(profiles) == 0.010


def test_choose_bin_width_falls_back_when_empty():
    assert choose_bin_width({}) == 0.01

