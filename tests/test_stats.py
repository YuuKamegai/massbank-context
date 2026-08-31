"""統計の既知解テスト。

参照実装（scipy）が無いため、教科書の手計算値でしか守れない。
ここが静かに狂うと、並び順と結論が黙って変わる。
"""
import math

import pytest

from massbank_context.stats import (
    benjamini_hochberg,
    fisher_exact_greater,
    wilcoxon_signed_rank,
)


class TestFisher:
    def test_tea_tasting(self):
        """Fisher の紅茶の実験。2x2 = [[3,1],[1,3]]。"""
        assert fisher_exact_greater(k=3, n=4, K=4, N=8) == pytest.approx(17 / 70)

    def test_exactly_one_half(self):
        assert fisher_exact_greater(k=2, n=5, K=3, N=10) == pytest.approx(0.5)

    def test_k_zero_gives_one(self):
        assert fisher_exact_greater(k=0, n=5, K=3, N=10) == pytest.approx(1.0)

    def test_all_members_drawn(self):
        assert fisher_exact_greater(k=3, n=3, K=3, N=10) == pytest.approx(1 / 120)


class TestBenjaminiHochberg:
    def test_all_equal_after_correction(self):
        q = benjamini_hochberg([0.01, 0.02, 0.03, 0.04, 0.05])
        assert q == pytest.approx([0.05, 0.05, 0.05, 0.05, 0.05])

    def test_classic_example(self):
        q = benjamini_hochberg(
            [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205])
        assert q == pytest.approx(
            [0.008, 0.032, 0.0672, 0.0672, 0.0672, 0.08, 0.084571, 0.205],
            abs=1e-6)

    def test_order_is_preserved(self):
        q = benjamini_hochberg([0.05, 0.01])
        assert q[0] > q[1]

    def test_none_passes_through(self):
        q = benjamini_hochberg([0.01, None, 0.02])
        assert q[1] is None
        assert q[0] == pytest.approx(0.02)


class TestWilcoxon:
    def test_all_positive_n5_exact(self):
        p, method = wilcoxon_signed_rank([1.0, 2.0, 3.0, 4.0, 5.0])
        assert method == "exact"
        assert p == pytest.approx(0.0625)

    def test_all_positive_n6_exact(self):
        p, method = wilcoxon_signed_rank([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        assert method == "exact"
        assert p == pytest.approx(0.03125)

    def test_one_negative(self):
        p, method = wilcoxon_signed_rank([-1.0, 2.0, 3.0, 4.0, 5.0])
        assert method == "exact"
        assert p == pytest.approx(0.125)

    def test_zeros_are_dropped(self):
        p, method = wilcoxon_signed_rank([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        assert p == pytest.approx(0.0625)

    def test_too_few_returns_insufficient(self):
        p, method = wilcoxon_signed_rank([1.0, 2.0])
        assert p is None
        assert method == "insufficient"

    def test_large_n_uses_normal(self):
        values = [float(i) for i in range(1, 30)]
        p, method = wilcoxon_signed_rank(values)
        assert method == "normal"
        assert 0.0 < p < 0.001

    def test_ties_use_tie_adjusted_normal_variance(self):
        """同順位では平均順位の二乗和に基づく分散を使う。"""
        p, method = wilcoxon_signed_rank([1.0, 1.0, 1.0, 1.0, 1.0])
        assert method == "normal"
        # ranks=[3,3,3,3,3], variance=sum(r*r)/4=11.25, p=0.0368884257...
        assert p == pytest.approx(0.03688842570704988)
