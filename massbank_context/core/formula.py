"""組成式の解析と単一同位体質量の計算（純関数・stdlib のみ）。

core/ は leaf なので上位モジュールを import しない。
"""
from __future__ import annotations

import re

__all__ = [
    "FormulaError",
    "ELEMENT_MASSES",
    "ELECTRON_MASS",
    "parse_formula",
    "monoisotopic_mass",
    "hill_formula",
    "shift_composition",
]


class FormulaError(ValueError):
    """組成式として解釈できない入力。"""


# 各元素の最も豊富な同位体の質量（Da）
ELEMENT_MASSES: dict[str, float] = {
    "H": 1.00782503207,
    "C": 12.0,
    "N": 14.0030740048,
    "O": 15.99491461956,
    "P": 30.97376163,
    "S": 31.97207100,
    "F": 18.99840322,
    "Cl": 34.96885268,
    "Br": 78.9183371,
    "I": 126.904473,
    "Na": 22.9897692809,
    "K": 38.96370668,
    "Si": 27.9769265325,
    "B": 11.0093054,
    "Se": 79.9165213,
}

ELECTRON_MASS = 0.00054857990946

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(text: str) -> dict[str, int]:
    """組成式を {元素: 個数} に直す。元素の並び順は問わず、重複は足し合わせる。

    括弧・電荷・小文字始まり・未知の元素は ``FormulaError``。
    """
    if not text:
        raise FormulaError("組成式が空")
    comp: dict[str, int] = {}
    pos = 0
    while pos < len(text):
        m = _TOKEN.match(text, pos)
        if m is None:
            raise FormulaError(f"組成式を解釈できない: {text!r}")
        el, num = m.group(1), m.group(2)
        if el not in ELEMENT_MASSES:
            raise FormulaError(f"未対応の元素 {el!r}: {text!r}")
        comp[el] = comp.get(el, 0) + (int(num) if num else 1)
        pos = m.end()
    return comp


def monoisotopic_mass(composition: dict[str, int]) -> float:
    """単一同位体質量（中性分子）。"""
    return sum(ELEMENT_MASSES[el] * n for el, n in composition.items())


def hill_formula(composition: dict[str, int]) -> str:
    """Hill 順の組成式。炭素があれば C, H を先頭に残りはアルファベット順、無ければ全てアルファベット順。"""
    items = {el: n for el, n in composition.items() if n > 0}
    if "C" in items:
        order = [el for el in ("C", "H") if el in items]
        order += sorted(el for el in items if el not in ("C", "H"))
    else:
        order = sorted(items)
    return "".join(el + (str(items[el]) if items[el] != 1 else "") for el in order)


def shift_composition(
    base: dict[str, int], add: dict[str, int], remove: dict[str, int]
) -> dict[str, int] | None:
    """base に add を足し remove を引く。どれかの元素が負になれば None。0 の元素は落とす。"""
    out = dict(base)
    for el, n in add.items():
        out[el] = out.get(el, 0) + n
    for el, n in remove.items():
        out[el] = out.get(el, 0) - n
    if any(n < 0 for n in out.values()):
        return None
    return {el: n for el, n in out.items() if n > 0}
