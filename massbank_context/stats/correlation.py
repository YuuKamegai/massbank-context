"""順位相関。stdlib のみ・純関数。"""
from __future__ import annotations

import math
from typing import Sequence


def _average_ranks(values: Sequence[float]) -> list[float]:
    """同順位は平均順位（1 始まり）。"""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def spearman_rho(x: Sequence[float], y: Sequence[float]) -> float | None:
    """Spearman の順位相関。同順位は平均順位にして Pearson を取る。

    n < 2 か、どちらかの順位の分散が 0 なら None（相関は定義できない）。
    """
    if len(x) != len(y):
        raise ValueError('x と y の長さが違います')
    if len(x) < 2:
        return None
    rx, ry = _average_ranks(x), _average_ranks(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0 or syy == 0:
        return None
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    return sxy / math.sqrt(sxx * syy)
