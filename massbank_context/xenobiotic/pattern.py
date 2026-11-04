"""出方の筋。群の偏りと、親との連動を試料ごとの値から見る。"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from massbank_context.core.evidence_models import Sample
from massbank_context.stats.correlation import spearman_rho

__all__ = ["SampleSides", "sample_sides", "detection_counts", "group_bias",
           "parent_rho", "parent_linked"]


@dataclass(frozen=True)
class SampleSides:
    test: tuple[str, ...]
    reference: tuple[str, ...]


def sample_sides(samples: Sequence[Sample], test_side: str) -> SampleSides:
    """test_side（'a' / 'b'）を test、もう一方を reference に分ける。group が空の試料は入れない。"""
    if test_side not in ('a', 'b'):
        raise ValueError("test_side は 'a' か 'b' です")
    other = 'b' if test_side == 'a' else 'a'
    return SampleSides(
        test=tuple(s.sample_id for s in samples if s.group == test_side),
        reference=tuple(s.sample_id for s in samples if s.group == other))


def _count(row: dict, ids: Sequence[str]) -> int | None:
    flags = [row[sid][1] for sid in ids if sid in row]
    if all(flag is None for flag in flags):
        return None  # 検出が不明なら「0 件」ではなく「分からない」
    return sum(1 for flag in flags if flag is True)


def detection_counts(row: dict[str, tuple[float | None, bool | None]],
                     sides: SampleSides) -> tuple[int | None, int | None]:
    """(test の検出数, reference の検出数)。その側の detected が全部 None なら None。"""
    return _count(row, sides.test), _count(row, sides.reference)


def group_bias(counts: tuple[int | None, int | None], sides: SampleSides,
               received_by: str, thresholds: dict) -> bool | None:
    """received_by 側にだけ出ているか。both か数が None なら判定できないので None。"""
    n_test, n_ref = counts
    if received_by not in ('test', 'reference') or n_test is None or n_ref is None:
        return None
    limit = thresholds['bias_reference_max']
    if received_by == 'test':
        return n_test > len(sides.test) / 2 and n_ref <= limit
    return n_ref > len(sides.reference) / 2 and n_test <= limit


def parent_rho(child: dict, parent: dict, sides: SampleSides) -> tuple[float | None, int]:
    """test 側の試料のうち両方に値があり、どちらも detected が False でない試料の ρ と n。"""
    xs, ys = [], []
    for sid in sides.test:
        if sid not in child or sid not in parent:
            continue
        (cv, cd), (pv, pd) = child[sid], parent[sid]
        if cv is None or pv is None or cd is False or pd is False:
            continue
        xs.append(cv)
        ys.append(pv)
    return spearman_rho(xs, ys), len(xs)


def parent_linked(rho: float | None, n: int, thresholds: dict) -> bool:
    """n と ρ が閾値を満たすとき連動とみなす。"""
    return rho is not None and n >= thresholds['rho_min_n'] and rho >= thresholds['rho_min']
