"""片側 Fisher 正確検定（濃縮方向のみ）。"""
from __future__ import annotations

from math import comb


def fisher_exact_greater(k: int, n: int, K: int, N: int) -> float:
    """「k 以上が当たる」確率を返す（片側・濃縮方向）。"""
    if not (0 <= k <= n <= N) or not (0 <= K <= N):
        raise ValueError(f"2x2 の値が不整合です: k={k} n={n} K={K} N={N}")
    upper = min(n, K)
    if k > upper:
        return 0.0
    total = comb(N, n)
    if total == 0:
        return 1.0
    hits = sum(
        comb(K, i) * comb(N - K, n - i)
        for i in range(k, upper + 1)
        if 0 <= n - i <= N - K
    )
    return hits / total
