"""付加体・質量差の表と、組成式から理論 m/z を並べる関数（純関数）。

表を変えたら版（``*_VERSION``）を上げる。結果にこの版が残るので、表の違いで
結果が変わったことを後から追える。
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from massbank_context.core.formula import (
    ELECTRON_MASS,
    monoisotopic_mass,
    parse_formula,
    shift_composition,
)

__all__ = [
    "ADDUCT_TABLE_VERSION",
    "SHIFT_TABLE_VERSION",
    "THRESHOLDS_VERSION",
    "THRESHOLDS",
    "Adduct",
    "ADDUCTS",
    "SOURCE_LOSSES",
    "Shift",
    "SHIFTS",
    "IonForm",
    "ion_forms",
    "ion_label",
]

ADDUCT_TABLE_VERSION = "xeno-adducts-1"
SHIFT_TABLE_VERSION = "xeno-shifts-1"
THRESHOLDS_VERSION = "xeno-thresholds-2026-11-05"

# 実データで較正した値（spec §10）。陽性の期待を最も多く満たし、おとりで A・B が 0 の範囲で最も厳しい値。
THRESHOLDS = {
    "msms_shared_min": 0.15,
    "msms_top_n": 10,
    "bias_reference_max": 1,
    "rho_min": 0.6,
    "rho_min_n": 5,
}

# 理論 m/z がこの幅（Da）未満で一致する形は 1 つに畳む
_SAME_MZ = 1e-6


@dataclass(frozen=True)
class Adduct:
    """m/z = (n·M + mass(add) − mass(remove) − charge·e) / |charge|。charge は符号付き。"""
    name: str
    polarity: str
    n_molecules: int
    charge: int
    add: str
    remove: str


ADDUCTS: tuple[Adduct, ...] = (
    Adduct("[M+H]+", "positive", 1, 1, "H", ""),
    Adduct("[M+Na]+", "positive", 1, 1, "Na", ""),
    Adduct("[M+K]+", "positive", 1, 1, "K", ""),
    Adduct("[M+NH4]+", "positive", 1, 1, "NH4", ""),
    Adduct("[2M+H]+", "positive", 2, 1, "H", ""),
    Adduct("[M+2H]2+", "positive", 1, 2, "H2", ""),
    Adduct("[M-H]-", "negative", 1, -1, "", "H"),
    Adduct("[M+Cl]-", "negative", 1, -1, "Cl", ""),
    Adduct("[M+HCOO]-", "negative", 1, -1, "CHO2", ""),
    Adduct("[M+CH3COO]-", "negative", 1, -1, "C2H3O2", ""),
    Adduct("[2M-H]-", "negative", 2, -1, "", "H"),
)

SOURCE_LOSSES = ("H2O", "NH3", "CO2")


@dataclass(frozen=True)
class Shift:
    name: str
    add: str
    remove: str


SHIFTS: tuple[Shift, ...] = (
    Shift("demethylation", "", "CH2"),
    Shift("oxidation", "O", ""),
    Shift("dioxidation", "O2", ""),
    Shift("deethylation", "", "C2H4"),
    Shift("dehydrogenation", "", "H2"),
    Shift("hydrogenation", "H2", ""),
    Shift("hydrolysis", "H2O", ""),
    Shift("glucuronidation", "C6H8O6", ""),
    Shift("sulfation", "SO3", ""),
    Shift("glutathione_conjugation", "C10H15N3O6S", ""),
    Shift("n_acetylation", "C2H2O", ""),
    Shift("methylation", "CH2", ""),
    Shift("glycine_conjugation", "C2H3NO", ""),
    Shift("taurine_conjugation", "C2H5NO2S", ""),
    Shift("dechlorination_hydroxylation", "OH", "Cl"),
)


@dataclass(frozen=True)
class IonForm:
    shifts: tuple[str, ...]
    ion: str
    mz: float

    @property
    def relation(self) -> str:
        """``<質量差>|<イオン>``。質量差が無ければ ``parent``。"""
        return f"{','.join(self.shifts) or 'parent'}|{self.ion}"


def ion_label(adduct: str, loss: str | None) -> str:
    """付加体名に損失を入れたラベル。損失は最後の ``]`` の前に ``-<損失>`` で入れる。"""
    if not loss:
        return adduct
    head, tail = adduct.rsplit("]", 1)
    return f"{head}-{loss}]{tail}"


def _mass(formula: str) -> float:
    return monoisotopic_mass(parse_formula(formula)) if formula else 0.0


def ion_forms(composition: dict[str, int], polarity: str, max_shifts: int) -> list[IonForm]:
    """組成式から、質量差（最大 ``max_shifts`` 個の組み合わせ）× 付加体 × 損失の理論 m/z を並べる。

    質量差は重複なしの組み合わせ。正味の質量差が 1e-6 Da 未満の組と、元素数が負になる組は作らない。
    理論 m/z が 1e-6 Da 未満で一致する形は 1 つに畳み、質量差の数が少ない → 損失なし →
    付加体の表の順 → 質量差の表の順で最初のものを残す。
    """
    adducts = [(i, a) for i, a in enumerate(ADDUCTS) if a.polarity == polarity]
    losses: list[str | None] = [None, *SOURCE_LOSSES]

    candidates: list[tuple[tuple[int, int, int, int], IonForm]] = []
    for k in range(max_shifts + 1):
        for ci, combo in enumerate(combinations(SHIFTS, k)):
            # 表の順に 1 つずつ適用し、途中で元素数が負になる組は作らない
            shifted: dict[str, int] | None = dict(composition)
            net = 0.0
            for shift in combo:
                add = parse_formula(shift.add) if shift.add else {}
                remove = parse_formula(shift.remove) if shift.remove else {}
                shifted = shift_composition(shifted, add, remove)
                if shifted is None:
                    break
                net += monoisotopic_mass(add) - monoisotopic_mass(remove)
            if shifted is None or (combo and abs(net) < _SAME_MZ):
                continue
            mass = monoisotopic_mass(shifted)
            names = tuple(s.name for s in combo)
            for li, loss in enumerate(losses):
                loss_mass = _mass(loss) if loss else 0.0
                for ai, adduct in adducts:
                    mz = (
                        adduct.n_molecules * mass + _mass(adduct.add) - _mass(adduct.remove)
                        - adduct.charge * ELECTRON_MASS - loss_mass
                    ) / abs(adduct.charge)
                    candidates.append(((k, li, ai, ci), IonForm(names, ion_label(adduct.name, loss), mz)))

    candidates.sort(key=lambda item: item[0])
    kept: list[IonForm] = []
    buckets: dict[int, list[float]] = {}
    for _, form in candidates:
        b = int(form.mz / _SAME_MZ)
        if any(abs(form.mz - other) < _SAME_MZ for nb in (b - 1, b, b + 1) for other in buckets.get(nb, ())):
            continue
        buckets.setdefault(b, []).append(form.mz)
        kept.append(form)
    return kept
