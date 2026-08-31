"""Wilcoxon 符号順位検定（中央値 0 に対する一標本・両側）。"""
from __future__ import annotations

import math

MIN_N = 5
EXACT_MAX_N = 20


def _exact_two_sided(ranks: list[float], w_min: float) -> float:
    """符号の全割り当てを数え、W の小さい側の両側確率を返す。"""
    counts = {0.0: 1}
    for rank in ranks:
        nxt: dict[float, int] = {}
        for total, times in counts.items():
            nxt[total] = nxt.get(total, 0) + times
            nxt[total + rank] = nxt.get(total + rank, 0) + times
        counts = nxt
    total_ways = 2 ** len(ranks)
    at_or_below = sum(
        times for total, times in counts.items() if total <= w_min + 1e-9
    )
    return min(1.0, 2.0 * at_or_below / total_ways)


def wilcoxon_signed_rank(values: list[float]) -> tuple[float | None, str]:
    """両側 p 値と、使った方法（exact / normal / insufficient）を返す。"""
    nonzero = [
        v for v in values
        if v is not None and not math.isnan(v) and v != 0.0
    ]
    n = len(nonzero)
    if n < MIN_N:
        return None, "insufficient"

    order = sorted(range(n), key=lambda i: abs(nonzero[i]))
    ranks = [0.0] * n
    position = 0
    has_ties = False
    while position < n:
        end = position
        while (
            end + 1 < n
            and abs(nonzero[order[end + 1]]) == abs(nonzero[order[position]])
        ):
            end += 1
        if end > position:
            has_ties = True
        average = (position + end) / 2 + 1
        for j in range(position, end + 1):
            ranks[order[j]] = average
        position = end + 1

    w_plus = sum(r for r, v in zip(ranks, nonzero) if v > 0)
    w_minus = sum(r for r, v in zip(ranks, nonzero) if v < 0)
    w_min = min(w_plus, w_minus)

    if n <= EXACT_MAX_N and not has_ties:
        return _exact_two_sided(ranks, w_min), "exact"

    mean = n * (n + 1) / 4
    # 同順位時は平均順位の二乗和を使い、順位の重複を分散へ反映する。
    variance = sum(r * r for r in ranks) / 4
    if variance <= 0:
        return None, "insufficient"
    z = (w_min - mean + 0.5) / math.sqrt(variance)
    p = 2.0 * (0.5 * math.erfc(-z / math.sqrt(2)))
    return min(1.0, max(0.0, p)), "normal"
