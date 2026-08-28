import pytest

from massbank_context.core.binning import (
    axis_entries,
    derive_neutral_losses,
    probe_bins,
    to_bin,
)
from massbank_context.core.classify import NEUTRAL_LOSS_CAPABLE
from massbank_context.core.models import Peak, SpectrumRecord


def _record(**overrides) -> SpectrumRecord:
    base = dict(
        accession="A", inchikey="IK", precursor_mz=230.0552,
        precursor_type="[M+H]+", polarity="POSITIVE",
        instrument_type="LC-ESI-ITFT", ms_type="MS2", resolution_value=15000.0,
        collision_energy=None, splash=None, license=None,
        molecular_formula=None, exact_mass=None,
        peaks=(Peak(86.0349, 11.0, "C2H4N3O+", None),
               Peak(171.0069, 226.0, None, None)),
    )
    base.update(overrides)
    return SpectrumRecord(**base)


def test_to_bin_rounds_to_nearest():
    assert to_bin(86.0349, 0.01) == 8603
    assert to_bin(86.0351, 0.01) == 8604


def test_probe_bins_covers_neighbours():
    assert probe_bins(86.0349, 0.01) == (8602, 8603, 8604)


def test_bin_boundary_is_recoverable():
    """境界を挟んだ 2 つの m/z が、互いの探索範囲に入ること。"""
    a, b = 86.0349, 86.0351          # 別ビンに落ちる近接ペア
    assert to_bin(a, 0.01) != to_bin(b, 0.01)
    assert to_bin(b, 0.01) in probe_bins(a, 0.01)
    assert to_bin(a, 0.01) in probe_bins(b, 0.01)


def test_derive_neutral_losses_for_lc_esi():
    peaks = derive_neutral_losses(_record())
    losses = [round(p.neutral_loss, 4) for p in peaks]
    assert losses == [144.0203, 59.0483]


def test_derive_neutral_losses_skips_gc_ei():
    """GC-EI-MS は precursor を持たないので中性損失は None のまま（spec §6.1）。"""
    rec = _record(instrument_type="EI-B", ms_type="MS",
                  precursor_mz=None, precursor_type=None)
    peaks = derive_neutral_losses(rec)
    assert all(p.neutral_loss is None for p in peaks)


def test_derive_neutral_losses_keeps_existing_values():
    """世代B が読み込んだ値を上書きしない。"""
    rec = _record(peaks=(Peak(135.0434, 999.0, None, 132.0421),))
    peaks = derive_neutral_losses(rec)
    assert peaks[0].neutral_loss == 132.0421


def test_axis_entries_emits_both_axes():
    entries = axis_entries(_record(), bin_width=0.01)
    axes = {axis for axis, _bin, _intensity, _formula in entries}
    assert axes == {"product_ion", "neutral_loss"}
    assert len(entries) == 4      # ピーク 2 本 × 2 軸


def test_axis_entries_omits_neutral_loss_for_gc_ei():
    rec = _record(instrument_type="EI-B", ms_type="MS",
                  precursor_mz=None, precursor_type=None)
    entries = axis_entries(rec, bin_width=0.01)
    axes = {axis for axis, _bin, _intensity, _formula in entries}
    assert axes == {"product_ion"}


def test_negative_neutral_loss_is_dropped():
    """precursor より重いピーク（同位体・多量体）は中性損失軸に載せない。"""
    rec = _record(peaks=(Peak(300.0, 10.0, None, None),))
    entries = axis_entries(rec, bin_width=0.01)
    axes = {axis for axis, _bin, _intensity, _formula in entries}
    assert axes == {"product_ion"}


@pytest.mark.parametrize("bin_width", [0.0, -0.01])
def test_invalid_bin_width_raises(bin_width):
    with pytest.raises(ValueError):
        to_bin(100.0, bin_width)


def test_neutral_loss_entries_carry_no_formula():
    """formula は product ion の組成式。中性損失の行に載せてはいけない。"""
    entries = axis_entries(_record(), bin_width=0.01)
    for axis, _bin, _intensity, formula in entries:
        if axis == "neutral_loss":
            assert formula is None
    # product_ion 側では失われていないこと
    product_formulas = [f for axis, _b, _i, f in entries if axis == "product_ion"]
    assert "C2H4N3O+" in product_formulas


def test_neutral_loss_skipped_when_precursor_missing_despite_capable_class():
    """MALDI-MS2 は NEUTRAL_LOSS_CAPABLE だが precursor を持たない記録が実在する。

    実データ: JEOL の MALDI-TOFTOF は ms_type=MS2 だが precursor_mz_value を持たない。
    """
    rec = _record(instrument_type="MALDI-TOFTOF", ms_type="MS2", precursor_mz=None)
    assert rec.acquisition_class in NEUTRAL_LOSS_CAPABLE   # 第 1 ゲートは通る
    peaks = derive_neutral_losses(rec)
    assert all(p.neutral_loss is None for p in peaks)      # 第 2 ゲートで止まる


def test_neutral_loss_skipped_for_incapable_class_even_with_precursor():
    """GC-EI-MS は precursor を持っていても中性損失軸を使わない。"""
    rec = _record(instrument_type="EI-B", ms_type="MS", precursor_mz=230.0552)
    assert rec.acquisition_class not in NEUTRAL_LOSS_CAPABLE
    peaks = derive_neutral_losses(rec)
    assert all(p.neutral_loss is None for p in peaks)
