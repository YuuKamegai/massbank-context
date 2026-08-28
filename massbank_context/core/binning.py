"""m/z のビン化と中性損失の算出（spec §6.3 / §6.4）。

固定幅ビンには境界問題がある。近接した 2 つの m/z が境界を挟むと別ビンに落ち、
同一フラグメントが一致しなくなる。ここでは round() で最近傍ビンに落とし、
照会時に ±1 ビンを探ることで境界損失を無くす。索引を二重化する必要がない。

依存グラフの leaf（classify / models のみ）。
"""
from __future__ import annotations

from massbank_context.core.classify import NEUTRAL_LOSS_CAPABLE
from massbank_context.core.models import Peak, SpectrumRecord

PRODUCT_ION = "product_ion"
NEUTRAL_LOSS = "neutral_loss"

#: 中性損失としてこれ未満は扱わない。precursor 自身のピーク（loss≒0）と
#: 同位体ピークを軸に載せないため。
_MIN_NEUTRAL_LOSS = 0.5


def to_bin(mz: float, bin_width: float) -> int:
    """m/z を最近傍のビン番号にする。"""
    if bin_width <= 0:
        raise ValueError(f"bin_width は正の数でなければならない: {bin_width}")
    return round(mz / bin_width)


def probe_bins(mz: float, bin_width: float) -> tuple[int, int, int]:
    """照会時に探るビン番号。境界を挟んだ相手を取りこぼさないための ±1。"""
    center = to_bin(mz, bin_width)
    return (center - 1, center, center + 1)


def derive_neutral_losses(record: SpectrumRecord) -> tuple[Peak, ...]:
    """中性損失が未設定のピークについて precursor_mz - mz を埋める。

    世代B が読み込んだ値は上書きしない。precursor を持たない acquisition_class
    （GC-EI-MS）では何もしない。
    """
    if (record.acquisition_class not in NEUTRAL_LOSS_CAPABLE
            or record.precursor_mz is None):
        return record.peaks

    filled: list[Peak] = []
    for peak in record.peaks:
        if peak.neutral_loss is not None:
            filled.append(peak)
            continue
        loss = record.precursor_mz - peak.mz
        filled.append(Peak(
            mz=peak.mz,
            relative_intensity=peak.relative_intensity,
            formula=peak.formula,
            neutral_loss=loss if loss >= _MIN_NEUTRAL_LOSS else None,
        ))
    return tuple(filled)


def axis_entries(record: SpectrumRecord,
                 bin_width: float) -> list[tuple[str, int, float, str | None]]:
    """レコードを (axis, bin, relative_intensity, formula) の列に展開する。"""
    peaks = derive_neutral_losses(record)
    entries: list[tuple[str, int, float, str | None]] = []
    for peak in peaks:
        entries.append((PRODUCT_ION, to_bin(peak.mz, bin_width),
                        peak.relative_intensity, peak.formula))
        if peak.neutral_loss is not None:
            entries.append((NEUTRAL_LOSS, to_bin(peak.neutral_loss, bin_width),
                            peak.relative_intensity, peak.formula))
    return entries
