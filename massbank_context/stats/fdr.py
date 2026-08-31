"""Benjamini-Hochberg の FDR 補正。"""
from __future__ import annotations

import math


def benjamini_hochberg(pvalues: list[float | None]) -> list[float | None]:
    """p 値の列を q 値へ変換する。入力の並びのまま返す。"""
    indexed = [
        (i, p) for i, p in enumerate(pvalues)
        if p is not None and not math.isnan(p)
    ]
    m = len(indexed)
    result: list[float | None] = [None] * len(pvalues)
    if m == 0:
        return result
    indexed.sort(key=lambda pair: pair[1])
    running = 1.0
    for rank in range(m, 0, -1):
        index, p = indexed[rank - 1]
        running = min(running, p * m / rank)
        result[index] = min(running, 1.0)
    return result
