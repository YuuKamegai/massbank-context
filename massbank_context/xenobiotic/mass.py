"""質量の筋: 登録した分子の理論 m/z と、特徴の m/z の照合（純関数）。

名前は照合キーにしない。使うのは組成式（理論 m/z）と特徴の m/z・付加体だけ。
"""
from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from collections import Counter
from dataclasses import dataclass
from typing import Sequence

from massbank_context.core.evidence_models import Feature
from massbank_context.xenobiotic.tables import (
    ADDUCTS,
    SOURCE_LOSSES,
    ion_forms,
    ion_label,
)

__all__ = [
    "MassRelation",
    "validate_ppm_bands",
    "tolerance_ppm",
    "adduct_recognized",
    "find_mass_relations",
]


@dataclass(frozen=True)
class MassRelation:
    feature_id: str
    entry_id: str
    shifts: tuple[str, ...]
    ion: str
    theoretical_mz: float
    ppm: float

    @property
    def relation(self) -> str:
        return f"{','.join(self.shifts) or 'parent'}|{self.ion}"


def validate_ppm_bands(ppm_bands: list[dict] | None) -> list[dict]:
    """ppm の帯（半開区間 [min_mz, max_mz)）を検証して min_mz 昇順で返す。誤りは ``ValueError``。"""
    bands: list[dict] = []
    for band in ppm_bands or []:
        ppm = band.get("ppm")
        lo, hi = band.get("min_mz"), band.get("max_mz")
        if not isinstance(ppm, (int, float)) or isinstance(ppm, bool) or not math.isfinite(ppm) or ppm <= 0:
            raise ValueError(f"ppm_bands の ppm は正の有限数: {band!r}")
        for bound in (lo, hi):
            if bound is not None and (not isinstance(bound, (int, float)) or not math.isfinite(bound)):
                raise ValueError(f"ppm_bands の min_mz / max_mz は有限数か null: {band!r}")
        if lo is not None and hi is not None and not lo < hi:
            raise ValueError(f"ppm_bands は min_mz < max_mz: {band!r}")
        bands.append({"min_mz": lo, "max_mz": hi, "ppm": float(ppm)})
    bands.sort(key=lambda b: -math.inf if b["min_mz"] is None else b["min_mz"])
    for prev, nxt in zip(bands, bands[1:]):
        prev_hi = math.inf if prev["max_mz"] is None else prev["max_mz"]
        next_lo = -math.inf if nxt["min_mz"] is None else nxt["min_mz"]
        if prev_hi > next_lo:
            raise ValueError(f"ppm_bands の帯が重なっている: {prev!r} と {nxt!r}")
    return bands


def tolerance_ppm(theoretical_mz: float, ppm: float, bands: list[dict]) -> float:
    """理論 m/z で帯を選ぶ。どの帯にも入らなければ ``ppm``。"""
    for band in bands:
        lo, hi = band["min_mz"], band["max_mz"]
        if (lo is None or theoretical_mz >= lo) and (hi is None or theoretical_mz < hi):
            return band["ppm"]
    return ppm


def _normalize_adduct(text: str) -> str:
    """空白を除き、MS-DIAL の ``]1+`` / ``]1-`` を ``]+`` / ``]-`` にそろえる。"""
    s = "".join(text.split())
    if s.endswith("]1+") or s.endswith("]1-"):
        s = s[:-3] + "]" + s[-1]
    return s


_KNOWN_LABELS = frozenset(
    ion_label(a.name, loss) for a in ADDUCTS for loss in (None, *SOURCE_LOSSES)
)


def adduct_recognized(text: str) -> bool:
    """特徴の付加体の記載が表のイオンとして読めるか（空白だけの記載は記載なしと同じ）。"""
    label = _normalize_adduct(text)
    return not label or label in _KNOWN_LABELS


def find_mass_relations(
    features: Sequence[Feature],
    entries: Sequence,
    polarity: str,
    ppm: float,
    ppm_bands: list[dict],
    max_shifts: int,
) -> tuple[dict[tuple[str, str], list[MassRelation]], Counter]:
    """各 entry の理論 m/z が窓に入る特徴を探す。鍵は (feature_id, entry_id)。

    各リストは (質量差の数, |ppm|) 順で、1 つの (特徴, entry) に複数あれば全部残す。
    特徴に付加体の記載があり、認識できて、関係のイオンと違えば落とす。認識できない
    記載は残して ``adduct_unrecognized`` を特徴の数で数える（``no_mz`` も特徴の数）。
    """
    counts: Counter = Counter()
    indexed: list[tuple[float, Feature]] = []
    adduct_of: dict[str, str | None] = {}
    for feature in features:
        if feature.mz is None:
            counts["no_mz"] += 1
            continue
        indexed.append((feature.mz, feature))
        label = _normalize_adduct(feature.adduct) if feature.adduct else None
        if label and label not in _KNOWN_LABELS:
            counts["adduct_unrecognized"] += 1
            label = None
        adduct_of[feature.feature_id] = label or None
    indexed.sort(key=lambda item: item[0])
    mzs = [mz for mz, _ in indexed]

    result: dict[tuple[str, str], list[MassRelation]] = {}
    for entry in entries:
        for form in ion_forms(entry.composition, polarity, max_shifts):
            tol = tolerance_ppm(form.mz, ppm, ppm_bands)
            lo = bisect_left(mzs, form.mz * (1 - tol * 1e-6))
            hi = bisect_right(mzs, form.mz * (1 + tol * 1e-6))
            for mz, feature in indexed[lo:hi]:
                label = adduct_of[feature.feature_id]
                if label is not None and label != form.ion:
                    counts["adduct_conflict_dropped"] += 1
                    continue
                result.setdefault((feature.feature_id, entry.entry_id), []).append(
                    MassRelation(
                        feature.feature_id, entry.entry_id, form.shifts, form.ion,
                        form.mz, (mz - form.mz) / form.mz * 1e6,
                    )
                )
    for relations in result.values():
        relations.sort(key=lambda r: (len(r.shifts), abs(r.ppm)))
    return result, counts
