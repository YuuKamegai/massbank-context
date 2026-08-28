import pytest

from massbank_context.core.binning import (
    axis_entries,
    derive_neutral_losses,
    probe_bins,
    to_bin,
)
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
