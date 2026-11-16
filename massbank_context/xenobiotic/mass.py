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
    "validate_ppm_offset",
    "offset_ppm",
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

    @property
    def has_loss(self) -> bool:
        """イオン源の中の損失つきの一致か（spec §4.1。水準 C の扱いが違う。§5）。"""
        return _BASE_ADDUCT.get(self.ion, self.ion) != self.ion


_BAND_KEYS = {"min_mz", "max_mz", "ppm"}


def validate_ppm_bands(ppm_bands: list[dict] | None) -> list[dict]:
    """ppm の帯（半開区間 [min_mz, max_mz)）を検証して min_mz 昇順で返す。誤りは ``ValueError``。"""
    bands: list[dict] = []
    for band in ppm_bands or []:
        unknown = set(band) - _BAND_KEYS
        if unknown:
            raise ValueError(f"ppm_bands に未知のキー: {sorted(unknown)}（使えるのは {sorted(_BAND_KEYS)}）")
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


def _finite(value) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def validate_ppm_offset(raw: dict | None) -> dict | None:
    """質量の窓の中心のずれを検証して正規化する。誤りは ``ValueError``。

    形は ms-data-parser の ``ppm_offset`` と同じ ``{"default": 0.0, "bands": [{"max_mz", "offset"}]}``
    （上流で決めた値をそのまま渡せるようにする）。中心は「実測 − 理論」の ppm の符号つきの値。
    帯は ``max_mz`` の昇順に並べ、理論 m/z が ``max_mz`` 未満の最初の帯を使う（``offset_ppm``）。
    """
    if raw is None:
        return None
    if not isinstance(raw, dict) or set(raw) - {"default", "bands"} or not _finite(raw.get("default", 0.0)):
        raise ValueError(f"ppm_offset は {{default, bands: [{{max_mz, offset}}]}}: {raw!r}")
    bands = []
    for band in raw.get("bands") or []:
        if (not isinstance(band, dict) or set(band) != {"max_mz", "offset"}
                or not _finite(band["max_mz"]) or not _finite(band["offset"])):
            raise ValueError(f"ppm_offset の帯は {{max_mz, offset}}（どちらも有限数）: {band!r}")
        bands.append({"max_mz": float(band["max_mz"]), "offset": float(band["offset"])})
    return {"default": float(raw.get("default", 0.0)),
            "bands": sorted(bands, key=lambda b: b["max_mz"])}


def offset_ppm(theoretical_mz: float, spec: dict | None) -> float:
    """理論 m/z の帯の中心のずれ（ppm）。境目の m/z は次の帯（上流と同じ ``<``）。"""
    if spec is None:
        return 0.0
    for band in spec["bands"]:
        if theoretical_mz < band["max_mz"]:
            return band["offset"]
    return spec["default"]


def _normalize_adduct(text: str) -> str:
    """空白を除き、MS-DIAL の ``]1+`` / ``]1-`` を ``]+`` / ``]-`` にそろえる。"""
    s = "".join(text.split())
    if s.endswith("]1+") or s.endswith("]1-"):
        s = s[:-3] + "]" + s[-1]
    return s


#: 表のイオンのラベル -> 損失を除いた付加体。spec §4.1 で損失は付加体の後に付くので、
#: [M+H-NH3]+ の付加体は [M+H]+ である。
_BASE_ADDUCT = {
    ion_label(a.name, loss): a.name for a in ADDUCTS for loss in (None, *SOURCE_LOSSES)
}


def adduct_recognized(text: str) -> bool:
    """特徴の付加体の記載が表のイオンとして読めるか（空白だけの記載は記載なしと同じ）。"""
    label = _normalize_adduct(text)
    return not label or label in _BASE_ADDUCT


def find_mass_relations(
    features: Sequence[Feature],
    entries: Sequence,
    polarity: str,
    ppm: float,
    ppm_bands: list[dict],
    max_shifts: int,
    low_score_adducts: frozenset[str] = frozenset(),
    ppm_offset: dict | None = None,
) -> tuple[dict[tuple[str, str], list[MassRelation]], Counter, Counter]:
    """各 entry の理論 m/z が窓に入る特徴を探す。鍵は (feature_id, entry_id)。

    窓は「理論 m/z × (1 + (中心 ± 許容幅) ppm)」。中心は ``ppm_offset``（正規化済み。
    ``validate_ppm_offset``）で、無ければ 0。行の ``ppm`` は中心に依らず実測 − 理論のまま。

    各リストは (質量差の数, |ppm|) 順で、1 つの (特徴, entry) に複数あれば全部残す。
    特徴に付加体の記載があり、認識できて、損失を除いた付加体同士が違えば落とす。認識できない
    記載は残して ``adduct_unrecognized`` を特徴の数で数える（``no_mz`` も特徴の数）。
    ``low_score_adducts`` の特徴の付加体は低得点の注釈の候補で実測の判定ではないので、比べずに
    ``adduct_low_score_ignored`` を特徴の数で数える。
    戻り値は (一致, 数, 組ごとの落とした一致の数)。組の鍵は一致と同じ。
    """
    counts: Counter = Counter()
    dropped: Counter = Counter()
    indexed: list[tuple[float, Feature]] = []
    adduct_of: dict[str, str | None] = {}
    for feature in features:
        if feature.mz is None:
            counts["no_mz"] += 1
            continue
        indexed.append((feature.mz, feature))
        label = _normalize_adduct(feature.adduct) if feature.adduct else None
        if label and feature.feature_id in low_score_adducts:
            counts["adduct_low_score_ignored"] += 1
            label = None
        if label and label not in _BASE_ADDUCT:
            counts["adduct_unrecognized"] += 1
            label = None
        adduct_of[feature.feature_id] = _BASE_ADDUCT[label] if label else None
    indexed.sort(key=lambda item: item[0])
    mzs = [mz for mz, _ in indexed]

    result: dict[tuple[str, str], list[MassRelation]] = {}
    for entry in entries:
        for form in ion_forms(entry.composition, polarity, max_shifts):
            tol = tolerance_ppm(form.mz, ppm, ppm_bands)
            center = offset_ppm(form.mz, ppm_offset)
            lo = bisect_left(mzs, form.mz * (1 + (center - tol) * 1e-6))
            hi = bisect_right(mzs, form.mz * (1 + (center + tol) * 1e-6))
            for mz, feature in indexed[lo:hi]:
                adduct = adduct_of[feature.feature_id]
                if adduct is not None and adduct != form.adduct:
                    counts["adduct_conflict_dropped"] += 1
                    dropped[(feature.feature_id, entry.entry_id)] += 1
                    continue
                result.setdefault((feature.feature_id, entry.entry_id), []).append(
                    MassRelation(
                        feature.feature_id, entry.entry_id, form.shifts, form.ion,
                        form.mz, (mz - form.mz) / form.mz * 1e6,
                    )
                )
    for relations in result.values():
        relations.sort(key=lambda r: (len(r.shifts), abs(r.ppm)))
    return result, counts, dropped
