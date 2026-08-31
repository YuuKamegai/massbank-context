# パスウェイ活性同定 実装計画（massbank-context 側）

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** InChIKey 付きの差次的結果を受け取り、PubChem のパスウェイを引いて、濃縮と方向を数値で返す 4 本の MCP ツールを作る。

**Architecture:** `stats/` に外部依存の無い純関数（Fisher・BH・Wilcoxon・集計）を置き、`experiment/` が契約 TSV を読み、`kg/` が SPARQL とキャッシュを持ち、`report/` が事実と仮説を分けた骨組みを組む。依存は `tools/ → report/ → {experiment/, kg/, stats/} → core/` の一方向。`stats/` は plain な dict と list だけを受け取り、`core/` にも依存しない。

**Tech Stack:** Python 3.14 / FastMCP / sqlite3 / urllib（**新しい依存は足さない**。現在の依存は `rdflib` と `mcp` の 2 つだけ）

**Spec:** `docs/superpowers/specs/2026-08-31-pathway-activity-design.md`

**関連計画:** 契約を書き出す側は `C:/Users/yuu18/Lipidmix_with_LLM/docs/superpowers/plans/2026-08-31-differential-export-for-pathway-activity.md`。**本計画はそちらの完成を待たない**（契約 TSV は合成 fixture で試験する）。

## Global Constraints

- Python は `C:/Python314/python.exe`（`python` は PATH に無いことがある）
- テストは `C:/Python314/python.exe -m pytest tests -q` をリポジトリルートから実行する
- コード・コメント・docstring・コミットメッセージはすべて**日本語**
- **依存を増やさない**。`rdflib` と `mcp` のみ。numpy も scipy も使わない
- JSON は `core/serialization.py` の `json_payload()` で返す。`json.dumps(..., indent=2)` を書かない
- 行が並ぶ一覧は `as_tsv()`（列名 1 回）、float は `round_floats()` で丸める
- 全ツールに `structured_output=False` を付けて `server.py` で登録する
- 可変状態（`config.INDEX_PATH` / `session`）は**モジュール修飾**で参照する。`from ... import INDEX_PATH` はスナップショット束縛になりテストの差し替えが効かない
- 打ち切り・スキップ・除外は必ず件数や理由を payload に出す（黙って減らさない）
- **契約の符号**: `log2fc` は正なら `group_b` が高い（上昇）。**massbank-context 側で反転しない**
- 「活性化」「不活性化」の語は `report/` の `hypotheses` の中にしか置かない

---

## File Structure

| ファイル | 責務 |
|---|---|
| `massbank_context/stats/__init__.py`（新規） | 再エクスポートのみ |
| `massbank_context/stats/fisher.py`（新規） | 片側 Fisher 正確検定（超幾何分布の上側和） |
| `massbank_context/stats/fdr.py`（新規） | Benjamini-Hochberg |
| `massbank_context/stats/wilcoxon.py`（新規） | Wilcoxon 符号順位検定（小 n は厳密・それ以外は正規近似） |
| `massbank_context/stats/activity.py`（新規） | パスウェイ単位の集計（純関数。plain dict を受けて plain dict を返す） |
| `massbank_context/experiment/__init__.py`（新規） | 再エクスポートのみ |
| `massbank_context/experiment/contract.py`（新規） | 契約 TSV のリーダとデータモデル |
| `massbank_context/kg/__init__.py`（新規） | 再エクスポートのみ |
| `massbank_context/kg/cache.py`（新規） | `pathway_cache.sqlite` のスキーマと読み書き |
| `massbank_context/kg/pathway.py`（新規） | PubChem SPARQL アダプタ（HTTP 層は差し替え可能） |
| `massbank_context/report/__init__.py`（新規） | 再エクスポートのみ |
| `massbank_context/report/findings.py`（新規） | 事実と仮説の骨組み |
| `massbank_context/core/errors.py`（修正） | `contract_mismatch` / `incompatible_organism` / `upstream_unavailable` |
| `massbank_context/core/session.py`（修正） | 差次的結果の保持 |
| `massbank_context/core/config.py`（修正） | キャッシュのパスとエンドポイント |
| `massbank_context/tools/pathway_tools.py`（新規） | `load_differential` / `resolve_pathways` / `pathway_activity` |
| `massbank_context/tools/report_tools.py`（新規） | `draft_pathway_findings` |
| `massbank_context/server.py`（修正） | 4 本の登録 |

`stats/` を 4 ファイルに割るのは、それぞれが独立した数値契約を持ち、既知解テストが 1 ファイルに集中すると読めなくなるためである。

---

## Task 1: `stats/` の 3 つの検定

**Files:**
- Create: `massbank_context/stats/__init__.py`, `massbank_context/stats/fisher.py`, `massbank_context/stats/fdr.py`, `massbank_context/stats/wilcoxon.py`
- Test: `tests/test_stats.py`

**Interfaces:**
- Consumes: なし（leaf）
- Produces:
  - `fisher_exact_greater(k: int, n: int, K: int, N: int) -> float`
  - `benjamini_hochberg(pvalues: list[float | None]) -> list[float | None]`
  - `wilcoxon_signed_rank(values: list[float]) -> tuple[float | None, str]`（`(p, method)`、method は `"exact"` / `"normal"` / `"insufficient"`）

- [ ] **Step 1: 既知解テストを書く**

`tests/test_stats.py` を新規作成する。**すべて手計算した値**である。

```python
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
        """Fisher の紅茶の実験。2x2 = [[3,1],[1,3]]。

        N=8, K=4, n=4, k=3。
        p = [C(4,3)C(4,1) + C(4,4)C(4,0)] / C(8,4) = (16+1)/70 = 17/70
        """
        assert fisher_exact_greater(k=3, n=4, K=4, N=8) == pytest.approx(17 / 70)

    def test_exactly_one_half(self):
        """N=10, K=3, n=5, k=2。

        p = [C(3,2)C(7,3) + C(3,3)C(7,2)] / C(10,5)
          = (3*35 + 1*21) / 252 = 126/252 = 0.5
        """
        assert fisher_exact_greater(k=2, n=5, K=3, N=10) == pytest.approx(0.5)

    def test_k_zero_gives_one(self):
        """k=0 なら「k 以上」は全事象なので p=1。"""
        assert fisher_exact_greater(k=0, n=5, K=3, N=10) == pytest.approx(1.0)

    def test_all_members_drawn(self):
        """N=10, K=3, n=3, k=3。p = C(3,3)C(7,0)/C(10,3) = 1/120"""
        assert fisher_exact_greater(k=3, n=3, K=3, N=10) == pytest.approx(1 / 120)


class TestBenjaminiHochberg:
    def test_all_equal_after_correction(self):
        """p=[0.01..0.05], m=5 は全て q=0.05 になる（手計算）。"""
        q = benjamini_hochberg([0.01, 0.02, 0.03, 0.04, 0.05])
        assert q == pytest.approx([0.05, 0.05, 0.05, 0.05, 0.05])

    def test_classic_example(self):
        """教科書例 m=8。単調化（後ろからの累積 min）が効くことを見る。"""
        q = benjamini_hochberg(
            [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205])
        assert q == pytest.approx(
            [0.008, 0.032, 0.0672, 0.0672, 0.0672, 0.08, 0.084571, 0.205],
            abs=1e-6)

    def test_order_is_preserved(self):
        """入力の並びで返す（ソートした順ではない）。"""
        q = benjamini_hochberg([0.05, 0.01])
        assert q[0] > q[1]

    def test_none_passes_through(self):
        """検定できなかった行は None のまま返し、m の分母にも数えない。"""
        q = benjamini_hochberg([0.01, None, 0.02])
        assert q[1] is None
        assert q[0] == pytest.approx(0.02)


class TestWilcoxon:
    def test_all_positive_n5_exact(self):
        """[1,2,3,4,5] は W+=15 で最大。両側 p = 2 * 1/32 = 0.0625"""
        p, method = wilcoxon_signed_rank([1.0, 2.0, 3.0, 4.0, 5.0])
        assert method == "exact"
        assert p == pytest.approx(0.0625)

    def test_all_positive_n6_exact(self):
        """n=6 で全て同符号: 2 * 1/64 = 0.03125"""
        p, method = wilcoxon_signed_rank([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        assert method == "exact"
        assert p == pytest.approx(0.03125)

    def test_one_negative(self):
        """[-1,2,3,4,5] は W-=1。W- <= 1 となる部分集合は {} と {1} の 2 通り。

        両側 p = 2 * 2/32 = 0.125
        """
        p, method = wilcoxon_signed_rank([-1.0, 2.0, 3.0, 4.0, 5.0])
        assert method == "exact"
        assert p == pytest.approx(0.125)

    def test_zeros_are_dropped(self):
        """0 は符号を持たないので落とす（n が減る）。"""
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
```

- [ ] **Step 2: 失敗を確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_stats.py -q`
Expected: FAIL（`ModuleNotFoundError: massbank_context.stats`）。

- [ ] **Step 3: Fisher を実装する**

`massbank_context/stats/fisher.py`:

```python
"""片側 Fisher 正確検定（濃縮方向のみ）。

超幾何分布の上側和をそのまま計算する。母集団は数千のオーダーで、
math.comb は多倍長整数なので桁落ちしない。scipy を持ち込まないための選択でもある。
"""
from __future__ import annotations

from math import comb


def fisher_exact_greater(k: int, n: int, K: int, N: int) -> float:
    """「k 以上が当たる」確率を返す（片側・濃縮方向）。

    引数:
        k: 差次的かつそのパスウェイのメンバーである数
        n: 差次的な数（背景の中で）
        K: そのパスウェイのメンバー数（背景の中で）
        N: 背景の化合物数

    両側にしないのは、パスウェイに差次的化合物が「集まっている」ことだけを
    問うているためである。枯渇（depletion）は主張しない。
    """
    if not (0 <= k <= n <= N) or not (0 <= K <= N):
        raise ValueError(f"2x2 の値が不整合です: k={k} n={n} K={K} N={N}")
    upper = min(n, K)
    if k > upper:
        return 0.0
    total = comb(N, n)
    if total == 0:
        return 1.0
    hits = sum(comb(K, i) * comb(N - K, n - i)
               for i in range(k, upper + 1)
               if 0 <= n - i <= N - K)
    return hits / total
```

- [ ] **Step 4: BH を実装する**

`massbank_context/stats/fdr.py`:

```python
"""Benjamini-Hochberg の FDR 補正。"""
from __future__ import annotations

import math


def benjamini_hochberg(pvalues: list[float | None]) -> list[float | None]:
    """p 値の列を q 値へ変換する。入力の並びのまま返す。

    None と NaN は「検定していない」を意味し、そのまま返して **m の分母にも数えない**。
    数えると、検定していない行の分だけ全体の q 値が甘くなる。
    """
    indexed = [(i, p) for i, p in enumerate(pvalues)
               if p is not None and not math.isnan(p)]
    m = len(indexed)
    result: list[float | None] = [None] * len(pvalues)
    if m == 0:
        return result
    indexed.sort(key=lambda pair: pair[1])
    running = 1.0
    # 後ろ（大きい p）から累積 min を取って単調化する。
    for rank in range(m, 0, -1):
        index, p = indexed[rank - 1]
        running = min(running, p * m / rank)
        result[index] = min(running, 1.0)
    return result
```

- [ ] **Step 5: Wilcoxon を実装する**

`massbank_context/stats/wilcoxon.py`:

```python
"""Wilcoxon 符号順位検定（中央値 0 に対する一標本・両側）。

n が小さいうちは厳密分布を数え上げ、大きくなったら正規近似に切り替える。
同順位（|値| が等しい）があると厳密分布の前提が崩れるため、そのときも近似を使う。
どちらを使ったかは戻り値の method で必ず開示する（黙って切り替えない）。
"""
from __future__ import annotations

import math

MIN_N = 5          # これ未満は両側 p が 0.05 を切れないので検定しない
EXACT_MAX_N = 20   # 2**20 = 約 100 万通り。ここまでは数え上げる


def _exact_two_sided(ranks: list[float], w_min: float) -> float:
    """符号のあらゆる割り当てを数え上げ、W の小さい側が w_min 以下になる確率×2。"""
    counts = {0.0: 1}
    for rank in ranks:
        nxt: dict[float, int] = {}
        for total, times in counts.items():
            nxt[total] = nxt.get(total, 0) + times
            nxt[total + rank] = nxt.get(total + rank, 0) + times
        counts = nxt
    total_ways = 2 ** len(ranks)
    at_or_below = sum(times for total, times in counts.items()
                      if total <= w_min + 1e-9)
    return min(1.0, 2.0 * at_or_below / total_ways)


def wilcoxon_signed_rank(values: list[float]) -> tuple[float | None, str]:
    """両側 p 値と、使った方法（exact / normal / insufficient）を返す。"""
    nonzero = [v for v in values
               if v is not None and not math.isnan(v) and v != 0.0]
    n = len(nonzero)
    if n < MIN_N:
        return None, "insufficient"

    order = sorted(range(n), key=lambda i: abs(nonzero[i]))
    ranks = [0.0] * n
    position = 0
    has_ties = False
    while position < n:
        end = position
        while (end + 1 < n
               and abs(nonzero[order[end + 1]]) == abs(nonzero[order[position]])):
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
    variance = n * (n + 1) * (2 * n + 1) / 24
    if variance <= 0:
        return None, "insufficient"
    z = (w_min - mean + 0.5) / math.sqrt(variance)
    p = 2.0 * (0.5 * math.erfc(-z / math.sqrt(2)))
    return min(1.0, max(0.0, p)), "normal"
```

- [ ] **Step 6: 再エクスポートを書く**

`massbank_context/stats/__init__.py`:

```python
"""統計。外部依存を持たない純関数だけを置く（core にも依存しない）。"""
from massbank_context.stats.activity import pathway_activity_rows
from massbank_context.stats.fdr import benjamini_hochberg
from massbank_context.stats.fisher import fisher_exact_greater
from massbank_context.stats.wilcoxon import wilcoxon_signed_rank

__all__ = ["benjamini_hochberg", "fisher_exact_greater",
           "pathway_activity_rows", "wilcoxon_signed_rank"]
```

**注**: `activity` は Task 2 で作る。このタスクの間は `__init__.py` から
`activity` の 2 行を外しておき、Task 2 で戻す。

- [ ] **Step 7: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_stats.py -q`
Expected: PASS（14 件）。

- [ ] **Step 8: commit**

```bash
git add massbank_context/stats tests/test_stats.py
git commit -m "feat: 外部依存の無い Fisher・BH・Wilcoxon を足す"
```

---

## Task 2: パスウェイ単位の集計（純関数）

**Files:**
- Create: `massbank_context/stats/activity.py`
- Modify: `massbank_context/stats/__init__.py`
- Test: `tests/test_activity.py`

**Interfaces:**
- Consumes: Task 1 の `fisher_exact_greater` / `benjamini_hochberg` / `wilcoxon_signed_rank`
- Produces:
  - `pathway_activity_rows(compounds: list[dict], memberships: dict[str, list[dict]], min_members: int) -> tuple[list[dict], dict]`
    - `compounds`: `[{"inchikey": str, "log2fc": float, "significant": bool}]`（**背景の化合物だけ**。パスウェイを持たない化合物は呼び出し側で除外済み）
    - `memberships`: `inchikey -> [{"pathway_id": str, "pathway_label": str, "taxid": str}]`
    - 戻り値 `(rows, summary)`。`rows` の各要素は
      `{pathway_id, pathway_label, taxid, n_members, n_significant, n_up, n_down, median_log2fc, direction, fisher_p, fisher_q, tested, wilcoxon_p, wilcoxon_tested}`
    - `summary` は `{"n_background": int, "n_significant": int, "n_pathways_total": int, "n_pathways_tested": int}`

- [ ] **Step 1: 失敗テストを書く**

`tests/test_activity.py` を新規作成する。

```python
"""パスウェイ単位の集計（純関数）。"""
import pytest

from massbank_context.stats.activity import pathway_activity_rows


def _compounds(spec):
    """spec は (inchikey, log2fc, significant) のリスト。"""
    return [{"inchikey": ik, "log2fc": fc, "significant": sig}
            for ik, fc, sig in spec]


def _members(mapping, pathway_id="PWID1", label="Citric Acid Cycle",
             taxid="TAXID10090"):
    return {ik: [{"pathway_id": pathway_id, "pathway_label": label,
                  "taxid": taxid}] for ik in mapping}


class TestDirection:
    def test_up_dominant(self):
        compounds = _compounds([("A", 2.0, True), ("B", 1.5, True),
                                ("C", 0.1, False)])
        rows, _ = pathway_activity_rows(
            compounds, _members(["A", "B", "C"]), min_members=3)
        assert rows[0]["direction"] == "UP_DOMINANT"
        assert rows[0]["n_up"] == 2
        assert rows[0]["n_down"] == 0

    def test_down_dominant(self):
        compounds = _compounds([("A", -2.0, True), ("B", -1.5, True),
                                ("C", 0.1, False)])
        rows, _ = pathway_activity_rows(
            compounds, _members(["A", "B", "C"]), min_members=3)
        assert rows[0]["direction"] == "DOWN_DOMINANT"

    def test_mixed(self):
        compounds = _compounds([("A", 2.0, True), ("B", -1.5, True),
                                ("C", 0.1, False)])
        rows, _ = pathway_activity_rows(
            compounds, _members(["A", "B", "C"]), min_members=3)
        assert rows[0]["direction"] == "MIXED"

    def test_single_member(self):
        compounds = _compounds([("A", 2.0, True), ("B", 0.1, False),
                                ("C", 0.1, False)])
        rows, _ = pathway_activity_rows(
            compounds, _members(["A", "B", "C"]), min_members=3)
        assert rows[0]["direction"] == "SINGLE_MEMBER"

    def test_no_significant_member(self):
        compounds = _compounds([("A", 0.2, False), ("B", 0.1, False),
                                ("C", 0.1, False)])
        rows, _ = pathway_activity_rows(
            compounds, _members(["A", "B", "C"]), min_members=3)
        assert rows[0]["direction"] == "NO_SIGNIFICANT_MEMBER"


class TestMedianUsesAllMembers:
    def test_median_includes_non_significant(self):
        """中央値は有意でないメンバーも含めた全員から取る。"""
        compounds = _compounds([("A", 3.0, True), ("B", 3.0, True),
                                ("C", 0.0, False)])
        rows, _ = pathway_activity_rows(
            compounds, _members(["A", "B", "C"]), min_members=3)
        assert rows[0]["median_log2fc"] == pytest.approx(3.0)
        assert rows[0]["n_members"] == 3


class TestTestingGate:
    def test_small_pathway_is_not_tested_but_kept(self):
        compounds = _compounds([("A", 2.0, True), ("B", 0.1, False),
                                ("C", 0.1, False)])
        memberships = {"A": [{"pathway_id": "PWID_SMALL",
                              "pathway_label": "Small", "taxid": "TAXID10090"}],
                       "B": [], "C": []}
        rows, summary = pathway_activity_rows(compounds, memberships,
                                              min_members=3)
        row = [r for r in rows if r["pathway_id"] == "PWID_SMALL"][0]
        assert row["tested"] is False
        assert row["fisher_p"] is None
        assert row["fisher_q"] is None
        assert summary["n_pathways_total"] == 1
        assert summary["n_pathways_tested"] == 0


class TestOrdering:
    def test_sorted_by_q_then_pathway_id(self):
        compounds = _compounds([("A", 2.0, True), ("B", 2.0, True),
                                ("C", 2.0, True), ("D", 0.0, False),
                                ("E", 0.0, False), ("F", 0.0, False)])
        memberships = {
            "A": [{"pathway_id": "PWID_B", "pathway_label": "B",
                   "taxid": "TAXID10090"}],
            "B": [{"pathway_id": "PWID_B", "pathway_label": "B",
                   "taxid": "TAXID10090"}],
            "C": [{"pathway_id": "PWID_B", "pathway_label": "B",
                   "taxid": "TAXID10090"}],
            "D": [{"pathway_id": "PWID_A", "pathway_label": "A",
                   "taxid": "TAXID10090"}],
            "E": [{"pathway_id": "PWID_A", "pathway_label": "A",
                   "taxid": "TAXID10090"}],
            "F": [{"pathway_id": "PWID_A", "pathway_label": "A",
                   "taxid": "TAXID10090"}],
        }
        rows, _ = pathway_activity_rows(compounds, memberships, min_members=3)
        # 差次的化合物が全部入っている PWID_B が先に来る
        assert rows[0]["pathway_id"] == "PWID_B"
        assert rows[0]["fisher_q"] <= rows[1]["fisher_q"]
```

- [ ] **Step 2: 失敗を確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_activity.py -q`
Expected: FAIL（`ModuleNotFoundError`）。

- [ ] **Step 3: 実装する**

`massbank_context/stats/activity.py`:

```python
"""パスウェイ単位の集計（純関数）。

plain な dict と list だけを受け取り、plain な dict を返す。experiment/ や kg/ の
データモデルを知らないので、stats/ は leaf のままでいられる。

濃縮と方向は**別々の量**であり、合成しない（spec §5.3）。並び順を決めるのは
濃縮の q 値だけで、方向の p 値では並べ替えない。
"""
from __future__ import annotations

from massbank_context.stats.fdr import benjamini_hochberg
from massbank_context.stats.fisher import fisher_exact_greater
from massbank_context.stats.wilcoxon import wilcoxon_signed_rank


def _median(values: list[float]) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2 == 1:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def _direction(n_up: int, n_down: int) -> str:
    """有意メンバーの上昇・低下の内訳から方向ラベルを決める。

    これは検定ではなく**記述**である。恣意的な閾値を置かないため、
    「片方が 0 かどうか」だけで判定する。
    """
    n_significant = n_up + n_down
    if n_significant == 0:
        return "NO_SIGNIFICANT_MEMBER"
    if n_significant == 1:
        return "SINGLE_MEMBER"
    if n_down == 0:
        return "UP_DOMINANT"
    if n_up == 0:
        return "DOWN_DOMINANT"
    return "MIXED"


def pathway_activity_rows(compounds: list[dict],
                          memberships: dict[str, list[dict]],
                          min_members: int) -> tuple[list[dict], dict]:
    """背景の化合物とパスウェイ所属から、パスウェイ単位の行を作る。

    引数:
        compounds: [{"inchikey", "log2fc", "significant"}]。**背景の化合物だけ**を渡す
            （パスウェイを 1 つも持たない化合物は呼び出し側で除外しておく）。
        memberships: inchikey -> [{"pathway_id", "pathway_label", "taxid"}]。
        min_members: この数未満のメンバーしか持たないパスウェイは検定しない。

    戻り値: (rows, summary)。rows は濃縮の q 値昇順、同値は pathway_id 昇順。
    """
    by_inchikey = {c["inchikey"]: c for c in compounds}
    n_background = len(by_inchikey)
    n_significant_total = sum(1 for c in by_inchikey.values()
                              if c.get("significant"))

    grouped: dict[str, dict] = {}
    for inchikey, entries in memberships.items():
        if inchikey not in by_inchikey:
            continue
        for entry in entries:
            pathway_id = entry["pathway_id"]
            bucket = grouped.setdefault(pathway_id, {
                "pathway_id": pathway_id,
                "pathway_label": entry.get("pathway_label", ""),
                "taxid": entry.get("taxid", ""),
                "inchikeys": set(),
            })
            bucket["inchikeys"].add(inchikey)

    rows: list[dict] = []
    for bucket in grouped.values():
        members = [by_inchikey[ik] for ik in bucket["inchikeys"]]
        log2fcs = [m["log2fc"] for m in members if m.get("log2fc") is not None]
        significant = [m for m in members if m.get("significant")]
        n_up = sum(1 for m in significant if (m.get("log2fc") or 0.0) > 0)
        n_down = sum(1 for m in significant if (m.get("log2fc") or 0.0) < 0)
        n_members = len(members)
        tested = n_members >= min_members
        fisher_p = None
        if tested:
            fisher_p = fisher_exact_greater(
                k=len(significant), n=n_significant_total,
                K=n_members, N=n_background)
        wilcoxon_p, method = (wilcoxon_signed_rank(log2fcs) if tested
                              else (None, "insufficient"))
        rows.append({
            "pathway_id": bucket["pathway_id"],
            "pathway_label": bucket["pathway_label"],
            "taxid": bucket["taxid"],
            "n_members": n_members,
            "n_significant": len(significant),
            "n_up": n_up,
            "n_down": n_down,
            "median_log2fc": _median(log2fcs),
            "direction": _direction(n_up, n_down),
            "fisher_p": fisher_p,
            "fisher_q": None,
            "tested": tested,
            "wilcoxon_p": wilcoxon_p,
            "wilcoxon_tested": method != "insufficient",
        })

    for row, q in zip(rows, benjamini_hochberg([r["fisher_p"] for r in rows])):
        row["fisher_q"] = q
    # Wilcoxon の BH は濃縮とは別に掛ける（合成しないため分けて補正する）。
    for row, q in zip(rows, benjamini_hochberg([r["wilcoxon_p"] for r in rows])):
        row["wilcoxon_q"] = q

    rows.sort(key=lambda r: (r["fisher_q"] if r["fisher_q"] is not None else 2.0,
                             r["pathway_id"]))
    summary = {
        "n_background": n_background,
        "n_significant": n_significant_total,
        "n_pathways_total": len(rows),
        "n_pathways_tested": sum(1 for r in rows if r["tested"]),
    }
    return rows, summary
```

- [ ] **Step 4: `__init__.py` に戻す**

Task 1 Step 6 で外した 2 行（`from massbank_context.stats.activity import pathway_activity_rows` と `__all__` の項目）を戻す。

- [ ] **Step 5: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_activity.py tests/test_stats.py -q`
Expected: PASS。

- [ ] **Step 6: commit**

```bash
git add massbank_context/stats tests/test_activity.py
git commit -m "feat: パスウェイ単位の濃縮と方向を集計する純関数を足す"
```

---

## Task 3: エラーエンベロープを 3 つ足す

**Files:**
- Modify: `massbank_context/core/errors.py`
- Test: `tests/test_errors.py`

**Interfaces:**
- Consumes: なし
- Produces: `contract_mismatch(reason, expected, path)` / `incompatible_organism(available, message)` / `upstream_unavailable(endpoint, message)`

- [ ] **Step 1: 失敗テストを書く**

`tests/test_errors.py` を新規作成する。

```python
"""エラーエンベロープ（spec §8・§9）。"""
import json

import pytest

from massbank_context.core import errors


def test_contract_mismatch_names_what_was_expected():
    payload = json.loads(errors.contract_mismatch(
        reason="必須列がありません: inchikey",
        expected="contract_version = 1",
        path="C:/tmp/differential.tsv"))
    assert payload["error"]["code"] == "contract_mismatch"
    assert payload["error"]["expected"] == "contract_version = 1"
    assert payload["error"]["path"].endswith("differential.tsv")


def test_incompatible_organism_lists_choices():
    payload = json.loads(errors.incompatible_organism(
        available=[{"taxid": "TAXID9606", "n_compounds": 67}],
        message="organism を指定してください。"))
    assert payload["error"]["code"] == "incompatible_organism"
    assert payload["error"]["available"][0]["taxid"] == "TAXID9606"


def test_upstream_unavailable_names_endpoint():
    payload = json.loads(errors.upstream_unavailable(
        endpoint="https://rdfportal.org/pubchem/sparql",
        message="接続できません。"))
    assert payload["error"]["code"] == "upstream_unavailable"
    assert "rdfportal" in payload["error"]["endpoint"]


@pytest.mark.parametrize("call", [
    lambda: errors.contract_mismatch("", "x", "y"),
    lambda: errors.incompatible_organism([], ""),
    lambda: errors.upstream_unavailable("", "m"),
])
def test_required_arguments(call):
    """手掛かりの無いエンベロープを黙って出さない。"""
    with pytest.raises(ValueError):
        call()
```

- [ ] **Step 2: 失敗を確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_errors.py -q`
Expected: FAIL（`AttributeError: module ... has no attribute 'contract_mismatch'`）。

- [ ] **Step 3: 実装する**

`massbank_context/core/errors.py` の末尾に追加する。

```python
def contract_mismatch(reason: str, expected: str, path: str) -> str:
    """入力ファイルが契約に合わない（spec §6.1）。

    どの契約のどこが合わなかったかを示さないと、書き出した側を直せない。
    """
    if not reason or not expected:
        raise ValueError("reason / expected は必須です。")
    return json_payload({
        "error": {
            "code": "contract_mismatch",
            "reason": reason,
            "expected": expected,
            "path": path,
        }
    })


def incompatible_organism(available: list[dict], message: str) -> str:
    """organism が未指定、または照会結果に無い（spec §6.3）。

    生物種を跨いだパスウェイを黙って混ぜない。既定で human に寄せることも
    しない — 選べる候補と件数を返して呼び出し側に選ばせる。
    """
    if not available or not message:
        raise ValueError("available / message は必須です。")
    return json_payload({
        "error": {
            "code": "incompatible_organism",
            "available": list(available),
            "message": message,
        }
    })


def upstream_unavailable(endpoint: str, message: str) -> str:
    """外部のエンドポイントに到達できない（spec §5.2）。"""
    if not endpoint or not message:
        raise ValueError("endpoint / message は必須です。")
    return json_payload({
        "error": {
            "code": "upstream_unavailable",
            "endpoint": endpoint,
            "message": message,
        }
    })
```

- [ ] **Step 4: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_errors.py -q`
Expected: PASS（6 件）。

- [ ] **Step 5: commit**

```bash
git add massbank_context/core/errors.py tests/test_errors.py
git commit -m "feat: 契約不一致・生物種不一致・上流到達不能のエラーを足す"
```

---

## Task 4: 契約リーダと `load_differential`

**Files:**
- Create: `massbank_context/experiment/__init__.py`, `massbank_context/experiment/contract.py`
- Create: `massbank_context/tools/pathway_tools.py`
- Modify: `massbank_context/core/session.py`
- Test: `tests/test_contract.py`, `tests/test_pathway_tools.py`

**Interfaces:**
- Consumes: Task 3 の `contract_mismatch`
- Produces:
  - `parse_differential(text: str) -> DifferentialSet`（`ContractError` を送出）
  - `DifferentialSet` は `.meta: dict[str, str]`、`.rows: list[dict]`、`.group_a`, `.group_b`, `.n_features_total`, `.n_unannotated` を持つ
  - `session.differential: DifferentialSet | None`
  - `load_differential(path: str) -> str`

- [ ] **Step 1: リーダの失敗テストを書く**

`tests/test_contract.py` を新規作成する。

```python
"""契約 TSV のリーダ（spec §6.1）。"""
import pytest

from massbank_context.experiment.contract import ContractError, parse_differential

GOOD = """# contract_version = 1
# exported_at = 2026-08-31T12:00:00+00:00
# source_arf = C:/x/AlignmentResult_2026_01_01_00_00_00.arf
# source_arf2 = C:/x/AlignmentResult_2026_01_01_00_00_00.arf2
# group_a = gf_AIN\tn_a = 5
# group_b = gf_HFD\tn_b = 5
# log2fc_sign = positive means group_b is higher
# q_threshold = 0.05\tlog2fc_threshold = 1.0\tlog_transform = true
# n_features_total = 10\tn_with_inchikey = 2\tn_unannotated = 8
spot_id\tname\tname_source\tontology\tinchikey\tinchikey_source\tmsi_level\tmz\trt\tlog2fc\tp_value\tq_value\tmean_a\tmean_b\tsignificant
1\tPC 34:1\tarf2\tPC\tAAAAAAAAAAAAAA-BBBBBBBBBB-C\tarf2\t\t760.5851\t12.34\t1.900000\t0.001\t0.01\t10\t40\ttrue
2\tCitric acid\tarf2\t\tKRKNYBCHXYNGOX-UHFFFAOYSA-N\tarf2\t\t191.0197\t2.10\t-0.500000\t0.4\t0.6\t20\t14\tfalse
"""


def test_parses_meta_and_rows():
    result = parse_differential(GOOD)
    assert result.group_a == "gf_AIN"
    assert result.group_b == "gf_HFD"
    assert result.n_features_total == 10
    assert result.n_unannotated == 8
    assert len(result.rows) == 2
    assert result.rows[0]["inchikey"] == "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
    assert result.rows[0]["log2fc"] == pytest.approx(1.9)
    assert result.rows[0]["significant"] is True
    assert result.rows[1]["significant"] is False


def test_rejects_unknown_contract_version():
    text = GOOD.replace("contract_version = 1", "contract_version = 99")
    with pytest.raises(ContractError) as excinfo:
        parse_differential(text)
    assert "contract_version" in str(excinfo.value)


def test_rejects_missing_sign_meta():
    """符号の向きが書かれていないファイルは読まない。

    向きが分からないまま読むと、上昇と低下が静かに入れ替わる。
    """
    text = "\n".join(line for line in GOOD.splitlines()
                     if not line.startswith("# log2fc_sign"))
    with pytest.raises(ContractError) as excinfo:
        parse_differential(text)
    assert "log2fc_sign" in str(excinfo.value)


def test_rejects_unexpected_sign_direction():
    text = GOOD.replace("positive means group_b is higher",
                        "positive means group_a is higher")
    with pytest.raises(ContractError):
        parse_differential(text)


def test_rejects_missing_column():
    text = GOOD.replace("\tinchikey\t", "\t")
    with pytest.raises(ContractError) as excinfo:
        parse_differential(text)
    assert "inchikey" in str(excinfo.value)


def test_rejects_empty_body():
    text = "\n".join(GOOD.splitlines()[:10])
    with pytest.raises(ContractError):
        parse_differential(text)
```

- [ ] **Step 2: 失敗を確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_contract.py -q`
Expected: FAIL（`ModuleNotFoundError`）。

- [ ] **Step 3: リーダを実装する**

`massbank_context/experiment/contract.py`:

```python
"""ms-data-parser が書き出す差次的結果を読む（spec §6.1 の契約）。

契約に合わないファイルは**読まずに落とす**。特に log2fc の向きは、書かれて
いなければ読まない — 向きが分からないまま読むと、上昇と低下が静かに入れ替わる。
"""
from __future__ import annotations

from dataclasses import dataclass, field

SUPPORTED_CONTRACT_VERSION = "1"
EXPECTED_SIGN = "positive means group_b is higher"

REQUIRED_COLUMNS = [
    "spot_id", "name", "ontology", "inchikey", "mz", "rt",
    "log2fc", "p_value", "q_value", "mean_a", "mean_b", "significant",
]
REQUIRED_META = ["contract_version", "group_a", "group_b", "log2fc_sign",
                 "n_features_total", "n_unannotated"]


class ContractError(ValueError):
    """契約に合わない入力。呼び出し側が contract_mismatch に変換する。"""


@dataclass
class DifferentialSet:
    meta: dict[str, str] = field(default_factory=dict)
    rows: list[dict] = field(default_factory=list)

    @property
    def group_a(self) -> str:
        return self.meta["group_a"]

    @property
    def group_b(self) -> str:
        return self.meta["group_b"]

    @property
    def n_features_total(self) -> int:
        return int(self.meta["n_features_total"])

    @property
    def n_unannotated(self) -> int:
        return int(self.meta["n_unannotated"])

    @property
    def inchikeys(self) -> list[str]:
        """重複を除いた InChIKey。並びは決定的にする。"""
        return sorted({row["inchikey"] for row in self.rows})


def _parse_meta_line(line: str, meta: dict[str, str]) -> None:
    """`# key = value` を、タブ区切りで複数個まで拾う。"""
    for chunk in line.lstrip("#").split("\t"):
        if "=" not in chunk:
            continue
        key, _, value = chunk.partition("=")
        meta[key.strip()] = value.strip()


def _to_float(value: str) -> float | None:
    value = (value or "").strip()
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_differential(text: str) -> DifferentialSet:
    """契約 TSV を読む。合わなければ ContractError。"""
    meta: dict[str, str] = {}
    body: list[str] = []
    for line in text.splitlines():
        if line.startswith("#"):
            _parse_meta_line(line, meta)
        elif line.strip():
            body.append(line)

    missing_meta = [key for key in REQUIRED_META if key not in meta]
    if missing_meta:
        raise ContractError(
            f"メタ行が足りません: {', '.join(missing_meta)}")
    if meta["contract_version"] != SUPPORTED_CONTRACT_VERSION:
        raise ContractError(
            f"未対応の contract_version です: {meta['contract_version']}"
            f"（対応: {SUPPORTED_CONTRACT_VERSION}）")
    if meta["log2fc_sign"] != EXPECTED_SIGN:
        raise ContractError(
            f"log2fc_sign が想定と違います: {meta['log2fc_sign']!r}"
            f"（想定: {EXPECTED_SIGN!r}）。向きが違うファイルは読みません。")
    if len(body) < 2:
        raise ContractError("ヘッダ行かデータ行がありません。")

    columns = body[0].split("\t")
    missing_columns = [c for c in REQUIRED_COLUMNS if c not in columns]
    if missing_columns:
        raise ContractError(
            f"必須列がありません: {', '.join(missing_columns)}")

    index = {name: i for i, name in enumerate(columns)}
    rows: list[dict] = []
    for line in body[1:]:
        cells = line.split("\t")
        if len(cells) != len(columns):
            raise ContractError(
                f"列数が合いません（期待 {len(columns)}、実際 {len(cells)}）: {line[:60]}")
        inchikey = cells[index["inchikey"]].strip()
        if not inchikey:
            raise ContractError(
                "inchikey が空の行があります。書き出し側で落とすべき行です。")
        rows.append({
            "spot_id": cells[index["spot_id"]].strip(),
            "name": cells[index["name"]].strip(),
            "ontology": cells[index["ontology"]].strip(),
            "inchikey": inchikey,
            "mz": _to_float(cells[index["mz"]]),
            "rt": _to_float(cells[index["rt"]]),
            "log2fc": _to_float(cells[index["log2fc"]]),
            "p_value": _to_float(cells[index["p_value"]]),
            "q_value": _to_float(cells[index["q_value"]]),
            "mean_a": _to_float(cells[index["mean_a"]]),
            "mean_b": _to_float(cells[index["mean_b"]]),
            "significant": cells[index["significant"]].strip().lower() == "true",
        })
    return DifferentialSet(meta=meta, rows=rows)
```

`massbank_context/experiment/__init__.py`:

```python
"""実験側の入力（ms-data-parser の契約）を読む層。"""
from massbank_context.experiment.contract import (
    ContractError,
    DifferentialSet,
    parse_differential,
)

__all__ = ["ContractError", "DifferentialSet", "parse_differential"]
```

- [ ] **Step 4: リーダのテストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_contract.py -q`
Expected: PASS（6 件）。

- [ ] **Step 5: セッションを拡張する**

`massbank_context/core/session.py` の `Session` に足す。

```python
        self.differential = None       # experiment.DifferentialSet | None
        self.differential_source: str | None = None
```

`reset()` にも同じ 2 行の初期化を足す。**型注釈のための import はしない**
（`core/` は leaf であり、`experiment/` を import してはならない）。

- [ ] **Step 6: `load_differential` の失敗テストを書く**

`tests/test_pathway_tools.py` を新規作成する。`GOOD` は `tests/test_contract.py` と
同じ内容をこのファイルにも書く（テスト間の暗黙の結合を作らない）。

```python
"""パスウェイ系ツール。"""
import json

import pytest

from massbank_context.core import session as session_module
from massbank_context.tools.pathway_tools import load_differential

GOOD = """# contract_version = 1
# exported_at = 2026-08-31T12:00:00+00:00
# group_a = gf_AIN\tn_a = 5
# group_b = gf_HFD\tn_b = 5
# log2fc_sign = positive means group_b is higher
# q_threshold = 0.05\tlog2fc_threshold = 1.0\tlog_transform = true
# n_features_total = 10\tn_with_inchikey = 2\tn_unannotated = 8
spot_id\tname\tname_source\tontology\tinchikey\tinchikey_source\tmsi_level\tmz\trt\tlog2fc\tp_value\tq_value\tmean_a\tmean_b\tsignificant
1\tPC 34:1\tarf2\tPC\tAAAAAAAAAAAAAA-BBBBBBBBBB-C\tarf2\t\t760.5851\t12.34\t1.900000\t0.001\t0.01\t10\t40\ttrue
2\tCitric acid\tarf2\t\tKRKNYBCHXYNGOX-UHFFFAOYSA-N\tarf2\t\t191.0197\t2.10\t-0.500000\t0.4\t0.6\t20\t14\tfalse
"""


@pytest.fixture(autouse=True)
def clean_session():
    session_module.session.reset()
    yield
    session_module.session.reset()


def test_load_differential_reports_counts(tmp_path):
    path = tmp_path / "differential.tsv"
    path.write_text(GOOD, encoding="utf-8")
    payload = json.loads(load_differential(str(path)))
    assert payload["n_rows"] == 2
    assert payload["n_inchikeys"] == 2
    assert payload["n_unannotated"] == 8
    assert payload["group_a"] == "gf_AIN"
    assert payload["group_b"] == "gf_HFD"
    assert "group_b" in payload["log2fc_sign"]
    assert session_module.session.differential is not None


def test_load_differential_rejects_bad_contract(tmp_path):
    path = tmp_path / "differential.tsv"
    path.write_text(GOOD.replace("contract_version = 1",
                                 "contract_version = 99"), encoding="utf-8")
    payload = json.loads(load_differential(str(path)))
    assert payload["error"]["code"] == "contract_mismatch"
    assert session_module.session.differential is None


def test_load_differential_missing_file(tmp_path):
    payload = json.loads(load_differential(str(tmp_path / "nope.tsv")))
    assert payload["error"]["code"] == "contract_mismatch"
    assert "見つかりません" in payload["error"]["reason"]


def test_payload_does_not_carry_rows(tmp_path):
    """行そのものは payload に載せない（spec §12）。"""
    path = tmp_path / "differential.tsv"
    path.write_text(GOOD, encoding="utf-8")
    payload = json.loads(load_differential(str(path)))
    assert "rows" not in payload
    assert "AAAAAAAAAAAAAA" not in json.dumps(payload)
```

- [ ] **Step 7: 失敗を確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_tools.py -q`
Expected: FAIL（`ModuleNotFoundError: massbank_context.tools.pathway_tools`）。

- [ ] **Step 8: ツールを実装する**

`massbank_context/tools/pathway_tools.py` を新規作成する（このタスクでは 1 本だけ）。

```python
"""パスウェイ系のツール。Phase 2〜3（spec §8）。"""
from __future__ import annotations

from pathlib import Path

from massbank_context.core.errors import contract_mismatch
from massbank_context.core.serialization import json_payload
from massbank_context.core.session import session
from massbank_context.experiment.contract import ContractError, parse_differential

__all__ = ["load_differential"]

_SIGN_NOTE = "log2fc は正なら group_b が高い（上昇）。massbank-context 側で反転しない。"


def load_differential(path: str) -> str:
    """ms-data-parser が書き出した差次的結果を読み込む（spec §6.1）。

    行そのものはセッションに保持し、payload には件数と来歴だけを返す。
    n_unannotated は「注釈が付かず書き出されなかった行」であり、
    **「変化が無かった」ではなく「調べていない」**行である。
    """
    file_path = Path(path)
    if not file_path.is_file():
        return contract_mismatch(
            reason=f"ファイルが見つかりません: {path}",
            expected="ms-data-parser の arf_export_differential が書き出した TSV",
            path=path)
    try:
        result = parse_differential(file_path.read_text(encoding="utf-8"))
    except ContractError as error:
        return contract_mismatch(
            reason=str(error),
            expected="contract_version = 1 / log2fc_sign = "
                     "positive means group_b is higher / 必須列一式",
            path=path)

    session.differential = result
    session.differential_source = str(file_path)
    return json_payload({
        "source": str(file_path),
        "group_a": result.group_a,
        "group_b": result.group_b,
        "n_a": result.meta.get("n_a"),
        "n_b": result.meta.get("n_b"),
        "n_rows": len(result.rows),
        "n_inchikeys": len(result.inchikeys),
        "n_features_total": result.n_features_total,
        "n_unannotated": result.n_unannotated,
        "n_significant": sum(1 for r in result.rows if r["significant"]),
        "log2fc_sign": _SIGN_NOTE,
        "note": ("n_unannotated は注釈が付かず書き出されなかった行数です。"
                 "「変化が無かった」ではなく「調べていない」行です。"),
    })
```

- [ ] **Step 9: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_tools.py tests/test_contract.py -q`
Expected: PASS（10 件）。

- [ ] **Step 10: commit**

```bash
git add massbank_context/experiment massbank_context/tools/pathway_tools.py massbank_context/core/session.py tests/test_contract.py tests/test_pathway_tools.py
git commit -m "feat: 差次的結果の契約リーダと load_differential を足す"
```

---

## Task 5: パスウェイキャッシュ

**Files:**
- Create: `massbank_context/kg/__init__.py`, `massbank_context/kg/cache.py`
- Modify: `massbank_context/core/config.py`
- Test: `tests/test_pathway_cache.py`

**Interfaces:**
- Consumes: なし
- Produces:
  - `open_cache(path) -> sqlite3.Connection`（スキーマを作る）
  - `pending_inchikeys(conn, inchikeys, endpoint, refresh=False) -> list[str]`
  - `record_result(conn, inchikey, endpoint, hits: list[dict]) -> None`（`hits` が空なら `status='empty'`）
  - `record_failure(conn, inchikey, endpoint) -> None`
  - `failed_inchikeys(conn, inchikeys) -> list[str]`
  - `memberships(conn, inchikeys, taxid=None) -> dict[str, list[dict]]`
  - `organism_counts(conn, inchikeys) -> list[dict]`（`[{"taxid", "label", "n_compounds"}]`、件数の降順）
  - `put_organism_labels(conn, labels: dict[str, str]) -> None`
  - `fetched_at_range(conn, inchikeys) -> tuple[str | None, str | None]`
  - `config.PATHWAY_CACHE_PATH`, `config.PUBCHEM_SPARQL_ENDPOINT`

- [ ] **Step 1: 失敗テストを書く**

`tests/test_pathway_cache.py` を新規作成する。

```python
"""パスウェイキャッシュ（spec §6.2）。"""
import pytest

from massbank_context.kg import cache

ENDPOINT = "https://rdfportal.org/pubchem/sparql"


@pytest.fixture
def conn(tmp_path):
    connection = cache.open_cache(tmp_path / "pathway_cache.sqlite")
    yield connection
    connection.close()


def _hit(pathway_id="PWID1", label="Citric Acid Cycle", taxid="TAXID10090"):
    return {"pathway_id": pathway_id, "pathway_label": label, "taxid": taxid}


class TestPending:
    def test_all_pending_at_first(self, conn):
        assert cache.pending_inchikeys(conn, ["A", "B"], ENDPOINT) == ["A", "B"]

    def test_recorded_hit_is_not_pending(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [_hit()])
        assert cache.pending_inchikeys(conn, ["A", "B"], ENDPOINT) == ["B"]

    def test_empty_result_is_not_pending(self, conn):
        """「照会したがヒット 0」は照会済みである。

        これを未照会と混同すると、被覆の分母が壊れる。
        """
        cache.record_result(conn, "A", ENDPOINT, [])
        assert cache.pending_inchikeys(conn, ["A"], ENDPOINT) == []

    def test_failure_stays_pending(self, conn):
        cache.record_failure(conn, "A", ENDPOINT)
        assert cache.pending_inchikeys(conn, ["A"], ENDPOINT) == ["A"]
        assert cache.failed_inchikeys(conn, ["A"]) == ["A"]

    def test_refresh_makes_everything_pending(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [_hit()])
        assert cache.pending_inchikeys(conn, ["A"], ENDPOINT,
                                       refresh=True) == ["A"]

    def test_other_endpoint_is_pending(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [_hit()])
        assert cache.pending_inchikeys(conn, ["A"], "http://other") == ["A"]


class TestMemberships:
    def test_returns_only_requested(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [_hit()])
        cache.record_result(conn, "B", ENDPOINT, [_hit(pathway_id="PWID2")])
        got = cache.memberships(conn, ["A"])
        assert list(got) == ["A"]
        assert got["A"][0]["pathway_id"] == "PWID1"

    def test_filters_by_taxid(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [
            _hit(pathway_id="PWID1", taxid="TAXID10090"),
            _hit(pathway_id="PWID2", taxid="TAXID562"),
        ])
        got = cache.memberships(conn, ["A"], taxid="TAXID10090")
        assert [e["pathway_id"] for e in got["A"]] == ["PWID1"]

    def test_refresh_replaces_previous_hits(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [_hit(pathway_id="PWID1")])
        cache.record_result(conn, "A", ENDPOINT, [_hit(pathway_id="PWID2")])
        got = cache.memberships(conn, ["A"])
        assert [e["pathway_id"] for e in got["A"]] == ["PWID2"]


class TestOrganismCounts:
    def test_counts_compounds_not_rows(self, conn):
        """同じ化合物が同じ生物種の複数パスウェイに属しても 1 と数える。"""
        cache.record_result(conn, "A", ENDPOINT, [
            _hit(pathway_id="PWID1", taxid="TAXID10090"),
            _hit(pathway_id="PWID2", taxid="TAXID10090"),
        ])
        cache.record_result(conn, "B", ENDPOINT, [
            _hit(pathway_id="PWID3", taxid="TAXID562")])
        counts = cache.organism_counts(conn, ["A", "B"])
        by_taxid = {c["taxid"]: c["n_compounds"] for c in counts}
        assert by_taxid == {"TAXID10090": 1, "TAXID562": 1}

    def test_labels_are_filled_when_known(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [_hit(taxid="TAXID10090")])
        cache.put_organism_labels(conn, {"TAXID10090": "Mus musculus"})
        assert cache.organism_counts(conn, ["A"])[0]["label"] == "Mus musculus"
```

- [ ] **Step 2: 失敗を確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_cache.py -q`
Expected: FAIL（`ModuleNotFoundError: massbank_context.kg`）。

- [ ] **Step 3: config にパスとエンドポイントを足す**

`massbank_context/core/config.py` の末尾に追加する。

```python
# パスウェイのキャッシュ。索引本体とは別ファイルにする（索引の再ビルドで消えない）。
PATHWAY_CACHE_PATH = Path(os.environ.get(
    "MASSBANK_CONTEXT_PATHWAY_CACHE", BASE_DIR / "index" / "pathway_cache.sqlite"))

PUBCHEM_SPARQL_ENDPOINT = os.environ.get(
    "MASSBANK_CONTEXT_PUBCHEM_SPARQL", "https://rdfportal.org/pubchem/sparql")
```

- [ ] **Step 4: キャッシュを実装する**

`massbank_context/kg/cache.py`:

```python
"""パスウェイのローカルキャッシュ（spec §6.2）。

lookup_log が要点である。「照会したがヒット 0」と「まだ照会していない」を
区別できないと被覆の分母が壊れる。pathway_hit に行が無いことだけでは
両者が同じに見えるため、照会した事実そのものを別テーブルに記録する。

日数ベースの TTL は持たない。PubChem のパスウェイは滅多に変わらない一方、
日付でキャッシュが勝手に切れると同じ解析が日によって違う結果を出す。
"""
from __future__ import annotations

import sqlite3
from datetime import datetime, timezone
from pathlib import Path

CACHE_CONTRACT_VERSION = "1"

_SCHEMA = """
CREATE TABLE IF NOT EXISTS pathway_hit (
    inchikey TEXT NOT NULL,
    pathway_id TEXT NOT NULL,
    pathway_label TEXT,
    taxid TEXT,
    PRIMARY KEY (inchikey, pathway_id, taxid)
);
CREATE INDEX IF NOT EXISTS idx_pathway_hit_inchikey ON pathway_hit(inchikey);
CREATE TABLE IF NOT EXISTS lookup_log (
    inchikey TEXT NOT NULL,
    endpoint TEXT NOT NULL,
    fetched_at TEXT NOT NULL,
    status TEXT NOT NULL,
    PRIMARY KEY (inchikey, endpoint)
);
CREATE TABLE IF NOT EXISTS organism (
    taxid TEXT PRIMARY KEY,
    label TEXT
);
CREATE TABLE IF NOT EXISTS cache_meta (
    key TEXT PRIMARY KEY,
    value TEXT
);
"""


def open_cache(path) -> sqlite3.Connection:
    """キャッシュを開く（無ければ作る）。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(path))
    conn.executescript(_SCHEMA)
    conn.execute(
        "INSERT OR IGNORE INTO cache_meta(key, value) VALUES ('contract_version', ?)",
        (CACHE_CONTRACT_VERSION,))
    conn.commit()
    return conn


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def pending_inchikeys(conn: sqlite3.Connection, inchikeys: list[str],
                      endpoint: str, refresh: bool = False) -> list[str]:
    """まだ照会が成功していない InChIKey を、入力の並びのまま返す。

    status が failed のものは pending に含める（次回の再試行対象）。
    """
    if refresh:
        return list(inchikeys)
    done = {row[0] for row in conn.execute(
        "SELECT inchikey FROM lookup_log WHERE endpoint = ? "
        "AND status IN ('ok', 'empty')", (endpoint,))}
    return [ik for ik in inchikeys if ik not in done]


def failed_inchikeys(conn: sqlite3.Connection,
                     inchikeys: list[str]) -> list[str]:
    """直近の照会が失敗している InChIKey。"""
    failed = {row[0] for row in conn.execute(
        "SELECT inchikey FROM lookup_log WHERE status = 'failed'")}
    return [ik for ik in inchikeys if ik in failed]


def record_result(conn: sqlite3.Connection, inchikey: str, endpoint: str,
                  hits: list[dict]) -> None:
    """照会結果を記録する。hits が空でも「照会した」ことを必ず残す。"""
    conn.execute("DELETE FROM pathway_hit WHERE inchikey = ?", (inchikey,))
    conn.executemany(
        "INSERT OR REPLACE INTO pathway_hit"
        "(inchikey, pathway_id, pathway_label, taxid) VALUES (?, ?, ?, ?)",
        [(inchikey, h["pathway_id"], h.get("pathway_label"), h.get("taxid"))
         for h in hits])
    conn.execute(
        "INSERT OR REPLACE INTO lookup_log"
        "(inchikey, endpoint, fetched_at, status) VALUES (?, ?, ?, ?)",
        (inchikey, endpoint, _now(), "ok" if hits else "empty"))
    conn.commit()


def record_failure(conn: sqlite3.Connection, inchikey: str,
                   endpoint: str) -> None:
    """照会に失敗したことを記録する（次回の再試行対象として残す）。"""
    conn.execute(
        "INSERT OR REPLACE INTO lookup_log"
        "(inchikey, endpoint, fetched_at, status) VALUES (?, ?, ?, 'failed')",
        (inchikey, endpoint, _now()))
    conn.commit()


def memberships(conn: sqlite3.Connection, inchikeys: list[str],
                taxid: str | None = None) -> dict[str, list[dict]]:
    """InChIKey ごとのパスウェイ所属。taxid を渡すとその生物種に絞る。"""
    if not inchikeys:
        return {}
    placeholders = ",".join("?" * len(inchikeys))
    sql = ("SELECT inchikey, pathway_id, pathway_label, taxid FROM pathway_hit "
           f"WHERE inchikey IN ({placeholders})")
    params = list(inchikeys)
    if taxid:
        sql += " AND taxid = ?"
        params.append(taxid)
    sql += " ORDER BY inchikey, pathway_id"
    result: dict[str, list[dict]] = {}
    for row in conn.execute(sql, params):
        result.setdefault(row[0], []).append({
            "pathway_id": row[1], "pathway_label": row[2], "taxid": row[3]})
    return result


def organism_counts(conn: sqlite3.Connection,
                    inchikeys: list[str]) -> list[dict]:
    """生物種ごとの**化合物数**（パスウェイ行数ではない）。件数の降順。"""
    if not inchikeys:
        return []
    placeholders = ",".join("?" * len(inchikeys))
    rows = conn.execute(
        "SELECT h.taxid, o.label, COUNT(DISTINCT h.inchikey) "
        "FROM pathway_hit h LEFT JOIN organism o ON o.taxid = h.taxid "
        f"WHERE h.inchikey IN ({placeholders}) "
        "GROUP BY h.taxid, o.label ORDER BY 3 DESC, 1 ASC",
        list(inchikeys)).fetchall()
    return [{"taxid": r[0], "label": r[1] or "", "n_compounds": r[2]}
            for r in rows]


def put_organism_labels(conn: sqlite3.Connection,
                        labels: dict[str, str]) -> None:
    """taxid のラベルを記録する（best effort。引けなくても解析は続く）。"""
    conn.executemany(
        "INSERT OR REPLACE INTO organism(taxid, label) VALUES (?, ?)",
        list(labels.items()))
    conn.commit()


def fetched_at_range(conn: sqlite3.Connection,
                     inchikeys: list[str]) -> tuple[str | None, str | None]:
    """照会日時の最小と最大（来歴として payload に出す）。"""
    if not inchikeys:
        return None, None
    placeholders = ",".join("?" * len(inchikeys))
    row = conn.execute(
        "SELECT MIN(fetched_at), MAX(fetched_at) FROM lookup_log "
        f"WHERE inchikey IN ({placeholders})", list(inchikeys)).fetchone()
    return row[0], row[1]
```

`massbank_context/kg/__init__.py`:

```python
"""InChIKey から外部のメタデータを引く層（spec §5.2）。"""
```

- [ ] **Step 5: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_cache.py -q`
Expected: PASS（11 件）。

- [ ] **Step 6: commit**

```bash
git add massbank_context/kg massbank_context/core/config.py tests/test_pathway_cache.py
git commit -m "feat: パスウェイのローカルキャッシュを足す（未照会とヒット0を区別する）"
```

---

## Task 6: SPARQL アダプタと `resolve_pathways`

**Files:**
- Create: `massbank_context/kg/pathway.py`
- Modify: `massbank_context/tools/pathway_tools.py`
- Test: `tests/test_pathway_adapter.py`, `tests/test_pathway_tools.py`（追記）

**Interfaces:**
- Consumes: Task 5 のキャッシュ、Task 3 の `upstream_unavailable`、Task 4 の `session.differential`
- Produces:
  - `build_pathway_query(inchikeys: list[str]) -> str`
  - `build_organism_label_query(taxids: list[str]) -> str`
  - `PubChemPathwayAdapter(endpoint: str, transport=None)`
    - `.fetch(inchikeys: list[str]) -> dict[str, list[dict]]`（**バッチが落ちたら 1 件ずつ再試行**し、落ちた InChIKey は戻り値に含めず `UpstreamError` も送出しない。失敗は `.failed` に溜まる）
    - `.failed: list[str]`
    - `.fetch_organism_labels(taxids) -> dict[str, str]`（best effort。失敗時は空 dict）
  - `UpstreamError`
  - `resolve_pathways(refresh: bool = False, batch_size: int = 50) -> str`

- [ ] **Step 1: アダプタの失敗テストを書く**

`tests/test_pathway_adapter.py` を新規作成する。**実エンドポイントは叩かない。**

```python
"""PubChem パスウェイのアダプタ。HTTP 層は差し替える。"""
import pytest

from massbank_context.kg.pathway import (
    PubChemPathwayAdapter,
    UpstreamError,
    build_pathway_query,
)

ENDPOINT = "https://rdfportal.org/pubchem/sparql"


def _binding(inchikey, pathway="PWID1", label="Citric Acid Cycle",
             organism="http://rdf.ncbi.nlm.nih.gov/pubchem/taxonomy/TAXID10090"):
    row = {
        "value_inchikey": {"type": "literal", "value": inchikey},
        "pathway": {"type": "uri",
                    "value": f"http://rdf.ncbi.nlm.nih.gov/pubchem/pathway/{pathway}"},
        "pathway_label": {"type": "literal", "value": label},
    }
    if organism:
        row["pathway_organism"] = {"type": "uri", "value": organism}
    return row


class TestQuery:
    def test_includes_all_inchikeys_and_pathway_graph(self):
        query = build_pathway_query(["A", "B"])
        assert '"A"' in query and '"B"' in query
        assert "rdf.ncbi.nlm.nih.gov/pubchem/pathway" in query
        assert "RO_0000057" in query

    def test_rejects_empty(self):
        with pytest.raises(ValueError):
            build_pathway_query([])


class TestFetch:
    def test_parses_bindings(self):
        def transport(endpoint, query):
            return {"results": {"bindings": [_binding("A")]}}

        adapter = PubChemPathwayAdapter(ENDPOINT, transport=transport)
        got = adapter.fetch(["A"])
        assert got["A"][0]["pathway_id"] == "PWID1"
        assert got["A"][0]["pathway_label"] == "Citric Acid Cycle"
        assert got["A"][0]["taxid"] == "TAXID10090"
        assert adapter.failed == []

    def test_missing_organism_becomes_empty_taxid(self):
        def transport(endpoint, query):
            return {"results": {"bindings": [_binding("A", organism=None)]}}

        adapter = PubChemPathwayAdapter(ENDPOINT, transport=transport)
        assert adapter.fetch(["A"])["A"][0]["taxid"] == ""

    def test_batch_failure_retries_individually(self):
        """バッチが落ちたら 1 件ずつ引き直し、どの InChIKey で落ちたかを特定する。"""
        calls = []

        def transport(endpoint, query):
            calls.append(query)
            if '"A"' in query and '"B"' in query:
                raise UpstreamError("batch failed")
            if '"B"' in query:
                raise UpstreamError("B is bad")
            return {"results": {"bindings": [_binding("A")]}}

        adapter = PubChemPathwayAdapter(ENDPOINT, transport=transport)
        got = adapter.fetch(["A", "B"])
        assert list(got) == ["A"]
        assert adapter.failed == ["B"]
        assert len(calls) == 3   # バッチ 1 回 + 個別 2 回

    def test_single_inchikey_failure_is_not_retried_twice(self):
        def transport(endpoint, query):
            raise UpstreamError("down")

        adapter = PubChemPathwayAdapter(ENDPOINT, transport=transport)
        assert adapter.fetch(["A"]) == {}
        assert adapter.failed == ["A"]
```

- [ ] **Step 2: 失敗を確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_adapter.py -q`
Expected: FAIL（`ModuleNotFoundError`）。

- [ ] **Step 3: アダプタを実装する**

`massbank_context/kg/pathway.py`:

```python
"""PubChem のパスウェイを SPARQL で引く（spec §5.2 の SparqlAdapter 実装）。

obo:RO_0000057 は **has participant** であり、その化合物が基質か産物か・
反応の向き・量論を持たない。ここで得られるのは「参加している」だけで、
活性そのものではない（spec §3.1）。

HTTP 層は transport で差し替えられるようにしてある。テストで実エンドポイントを
叩かないためと、腐敗防止テスト（別マーク）だけが実物を見るようにするため。
"""
from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request

DEFAULT_TIMEOUT = 120

_QUERY_TEMPLATE = """
PREFIX sio: <http://semanticscience.org/resource/>
PREFIX vocab: <http://rdf.ncbi.nlm.nih.gov/pubchem/vocabulary#>
PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX up: <http://purl.uniprot.org/core/>
SELECT DISTINCT ?value_inchikey ?pathway ?pathway_label ?pathway_organism
FROM <http://rdf.ncbi.nlm.nih.gov/pubchem/inchikey>
FROM <http://rdf.ncbi.nlm.nih.gov/pubchem/compound>
FROM <http://rdf.ncbi.nlm.nih.gov/pubchem/pathway>
WHERE {
  VALUES ?value_inchikey { %(values)s }
  ?inchikey_node sio:SIO_000300 ?value_inchikey ; sio:SIO_000011 ?cid .
  ?cid a vocab:Compound .
  ?pathway a vocab:Pathway ;
           obo:RO_0000057 ?cid ;
           dcterms:title ?pathway_label .
  OPTIONAL { ?pathway up:organism ?pathway_organism . }
}
"""

_ORGANISM_QUERY_TEMPLATE = """
PREFIX dcterms: <http://purl.org/dc/terms/>
SELECT DISTINCT ?taxonomy ?label
FROM <http://rdf.ncbi.nlm.nih.gov/pubchem/taxonomy>
WHERE {
  VALUES ?taxonomy { %(values)s }
  ?taxonomy dcterms:title ?label .
}
"""


class UpstreamError(RuntimeError):
    """外部エンドポイントへの照会が失敗した。"""


def _values(items: list[str], as_uri: bool = False) -> str:
    if as_uri:
        return " ".join(f"<{item}>" for item in items)
    return " ".join(json.dumps(item) for item in items)


def build_pathway_query(inchikeys: list[str]) -> str:
    """InChIKey から CID を経てパスウェイを引くクエリ。"""
    if not inchikeys:
        raise ValueError("inchikeys が空です。")
    return _QUERY_TEMPLATE % {"values": _values(inchikeys)}


def build_organism_label_query(taxids: list[str]) -> str:
    """taxonomy URI からラベルを引くクエリ。"""
    if not taxids:
        raise ValueError("taxids が空です。")
    uris = [f"http://rdf.ncbi.nlm.nih.gov/pubchem/taxonomy/{t}" for t in taxids]
    return _ORGANISM_QUERY_TEMPLATE % {"values": _values(uris, as_uri=True)}


def http_transport(endpoint: str, query: str) -> dict:
    """既定の HTTP 層。SPARQL results JSON を返す。"""
    request = urllib.request.Request(
        endpoint,
        data=urllib.parse.urlencode({"query": query}).encode(),
        headers={"Accept": "application/sparql-results+json"})
    try:
        with urllib.request.urlopen(request, timeout=DEFAULT_TIMEOUT) as response:
            return json.load(response)
    except (urllib.error.URLError, OSError, ValueError) as error:
        raise UpstreamError(str(error)) from error


def _tail(uri: str) -> str:
    return uri.rsplit("/", 1)[-1] if uri else ""


class PubChemPathwayAdapter:
    """バッチで引き、落ちたら 1 件ずつ引き直すアダプタ。"""

    def __init__(self, endpoint: str, transport=None) -> None:
        self.endpoint = endpoint
        self._transport = transport or http_transport
        self.failed: list[str] = []

    def _query(self, inchikeys: list[str]) -> dict[str, list[dict]]:
        payload = self._transport(self.endpoint, build_pathway_query(inchikeys))
        result: dict[str, list[dict]] = {}
        for binding in payload.get("results", {}).get("bindings", []):
            inchikey = binding["value_inchikey"]["value"]
            result.setdefault(inchikey, []).append({
                "pathway_id": _tail(binding["pathway"]["value"]),
                "pathway_label": binding.get("pathway_label", {}).get("value", ""),
                "taxid": _tail(binding.get("pathway_organism", {}).get("value", "")),
            })
        return result

    def fetch(self, inchikeys: list[str]) -> dict[str, list[dict]]:
        """1 バッチ分を引く。落ちたら 1 件ずつ引き直す。

        バッチ全体を諦めると、1 件の悪い InChIKey が 50 件を巻き添えにする。
        どれが悪いのかを特定できないと、再試行の対象も決められない。
        """
        if not inchikeys:
            return {}
        try:
            return self._query(inchikeys)
        except UpstreamError:
            if len(inchikeys) == 1:
                self.failed.extend(inchikeys)
                return {}
        result: dict[str, list[dict]] = {}
        for inchikey in inchikeys:
            try:
                result.update(self._query([inchikey]))
            except UpstreamError:
                self.failed.append(inchikey)
        return result

    def fetch_organism_labels(self, taxids: list[str]) -> dict[str, str]:
        """taxid のラベルを引く（best effort。失敗しても解析は続ける）。"""
        taxids = [t for t in taxids if t]
        if not taxids:
            return {}
        try:
            payload = self._transport(self.endpoint,
                                      build_organism_label_query(taxids))
        except UpstreamError:
            return {}
        return {_tail(b["taxonomy"]["value"]): b["label"]["value"]
                for b in payload.get("results", {}).get("bindings", [])}
```

- [ ] **Step 4: アダプタのテストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_adapter.py -q`
Expected: PASS（6 件）。

- [ ] **Step 5: `resolve_pathways` の失敗テストを書く**

`tests/test_pathway_tools.py` に追記する。冒頭の import に足す。

```python
from massbank_context.core import config
from massbank_context.tools.pathway_tools import resolve_pathways
```

`clean_session` フィクスチャを次に差し替える（キャッシュも毎回作り直す）。

```python
@pytest.fixture(autouse=True)
def clean_session(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "PATHWAY_CACHE_PATH",
                        tmp_path / "pathway_cache.sqlite")
    session_module.session.reset()
    yield
    session_module.session.reset()
```

テストを追記する。

```python
def _load(tmp_path, text=GOOD):
    path = tmp_path / "differential.tsv"
    path.write_text(text, encoding="utf-8")
    load_differential(str(path))


def _stub_adapter(hits, failed=None):
    class _Stub:
        def __init__(self):
            self.failed = list(failed or [])

        def fetch(self, inchikeys):
            return {ik: hits[ik] for ik in inchikeys if ik in hits}

        def fetch_organism_labels(self, taxids):
            return {"TAXID10090": "Mus musculus"}

    return _Stub()


def test_resolve_pathways_requires_differential():
    payload = json.loads(resolve_pathways())
    assert payload["error"]["code"] == "missing_state"
    assert "load_differential" in payload["error"]["required_tools"]


def test_resolve_pathways_reports_coverage(tmp_path, monkeypatch):
    _load(tmp_path)
    hits = {"KRKNYBCHXYNGOX-UHFFFAOYSA-N": [
        {"pathway_id": "PWID1", "pathway_label": "Citric Acid Cycle",
         "taxid": "TAXID10090"}]}
    monkeypatch.setattr("massbank_context.tools.pathway_tools._make_adapter",
                        lambda: _stub_adapter(hits))
    payload = json.loads(resolve_pathways())
    assert payload["n_queried"] == 2
    assert payload["n_with_pathway"] == 1
    assert payload["pathway_coverage"] == pytest.approx(0.5)
    assert payload["organisms"][0]["taxid"] == "TAXID10090"
    assert payload["organisms"][0]["label"] == "Mus musculus"
    assert payload["n_failed"] == 0


def test_resolve_pathways_reports_failures(tmp_path, monkeypatch):
    _load(tmp_path)
    monkeypatch.setattr(
        "massbank_context.tools.pathway_tools._make_adapter",
        lambda: _stub_adapter({}, failed=["KRKNYBCHXYNGOX-UHFFFAOYSA-N"]))
    payload = json.loads(resolve_pathways())
    assert payload["n_failed"] == 1
    assert "再実行" in payload["note"]


def test_resolve_pathways_reports_upstream_down(tmp_path, monkeypatch):
    """照会した全件が落ちたら、部分的な成功として返さず upstream_unavailable を返す。

    「0 件見つかった」と「エンドポイントが落ちていた」を同じ形で返すと、
    被覆 0% を結果として読まれる。
    """
    _load(tmp_path)
    monkeypatch.setattr(
        "massbank_context.tools.pathway_tools._make_adapter",
        lambda: _stub_adapter({}, failed=["AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                                          "KRKNYBCHXYNGOX-UHFFFAOYSA-N"]))
    payload = json.loads(resolve_pathways())
    assert payload["error"]["code"] == "upstream_unavailable"
    assert "rdfportal" in payload["error"]["endpoint"]


def test_resolve_pathways_is_cached(tmp_path, monkeypatch):
    _load(tmp_path)
    hits = {"KRKNYBCHXYNGOX-UHFFFAOYSA-N": [
        {"pathway_id": "PWID1", "pathway_label": "Citric Acid Cycle",
         "taxid": "TAXID10090"}]}
    monkeypatch.setattr("massbank_context.tools.pathway_tools._make_adapter",
                        lambda: _stub_adapter(hits))
    resolve_pathways()

    def _explode():
        raise AssertionError("2 回目はキャッシュで済むはず")

    monkeypatch.setattr("massbank_context.tools.pathway_tools._make_adapter",
                        _explode)
    payload = json.loads(resolve_pathways())
    assert payload["n_queried"] == 0
    assert payload["n_with_pathway"] == 1
```

- [ ] **Step 6: 失敗を確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_tools.py -q`
Expected: FAIL（`ImportError: cannot import name 'resolve_pathways'`）。

- [ ] **Step 7: `resolve_pathways` を実装する**

`massbank_context/tools/pathway_tools.py` に追記する。import と `__all__` を更新する。

```python
from massbank_context.core import config
from massbank_context.core.errors import (
    contract_mismatch,
    missing_state,
    upstream_unavailable,
)
from massbank_context.kg import cache
from massbank_context.kg.pathway import PubChemPathwayAdapter

__all__ = ["load_differential", "resolve_pathways"]

_BATCH_SIZE = 50


def _make_adapter():
    """アダプタの生成点を 1 か所にする（テストで差し替えられるようにするため）。"""
    return PubChemPathwayAdapter(config.PUBCHEM_SPARQL_ENDPOINT)


def resolve_pathways(refresh: bool = False) -> str:
    """読み込んだ差次的結果の InChIKey にパスウェイを引く（spec §6.2）。

    キャッシュに無いものだけを 50 件ずつ照会する。refresh=True で全件引き直す
    （日数ベースの TTL は持たない。日付でキャッシュが切れると同じ解析が
    日によって違う結果を出すため）。

    照会に失敗した InChIKey が残っている間は pathway_activity を走らせない。
    部分的な母集団で統計を回すと被覆が見かけ上変わり、濃縮の p 値が動く。
    """
    if session.differential is None:
        return missing_state(
            "differential", ["load_differential"],
            "差次的結果が読み込まれていません。先に load_differential を実行してください。")

    inchikeys = session.differential.inchikeys
    endpoint = config.PUBCHEM_SPARQL_ENDPOINT
    conn = cache.open_cache(config.PATHWAY_CACHE_PATH)
    try:
        pending = cache.pending_inchikeys(conn, inchikeys, endpoint,
                                          refresh=refresh)
        n_failed = 0
        if pending:
            adapter = _make_adapter()
            for start in range(0, len(pending), _BATCH_SIZE):
                batch = pending[start:start + _BATCH_SIZE]
                hits = adapter.fetch(batch)
                failed = set(adapter.failed)
                for inchikey in batch:
                    if inchikey in failed:
                        cache.record_failure(conn, inchikey, endpoint)
                    else:
                        cache.record_result(conn, inchikey, endpoint,
                                            hits.get(inchikey, []))
            if len(adapter.failed) >= len(pending):
                # 照会した全件が落ちた = 個々の InChIKey の問題ではなく上流が死んでいる。
                # 「0 件見つかった」と同じ形で返すと、被覆 0% を結果として読まれる。
                return upstream_unavailable(
                    endpoint=endpoint,
                    message=(f"照会した {len(pending)} 件すべてが失敗しました。"
                             "エンドポイントに到達できていない可能性があります。"
                             "復旧後に resolve_pathways を再実行してください。"))
            n_failed = len(cache.failed_inchikeys(conn, inchikeys))
            taxids = [o["taxid"] for o in cache.organism_counts(conn, inchikeys)
                      if o["taxid"] and not o["label"]]
            if taxids:
                cache.put_organism_labels(conn,
                                          adapter.fetch_organism_labels(taxids))
        else:
            n_failed = len(cache.failed_inchikeys(conn, inchikeys))

        found = cache.memberships(conn, inchikeys)
        organisms = cache.organism_counts(conn, inchikeys)
        first, last = cache.fetched_at_range(conn, inchikeys)
    finally:
        conn.close()

    n_with_pathway = len(found)
    n_resolved = len(inchikeys) - n_failed
    coverage = (n_with_pathway / n_resolved) if n_resolved else 0.0
    note = ("パスウェイが付かなかった化合物は「変化が無かった」のではなく"
            "「調べていない」ものです。")
    if n_failed:
        note += (f" 照会に失敗した InChIKey が {n_failed} 件あります。"
                 "resolve_pathways を再実行してください"
                 "（未完了のまま pathway_activity は走りません）。")
    return json_payload(round_floats({
        "endpoint": endpoint,
        "n_inchikeys": len(inchikeys),
        "n_queried": len(pending),
        "n_cached": len(inchikeys) - len(pending),
        "n_with_pathway": n_with_pathway,
        "n_failed": n_failed,
        "pathway_coverage": coverage,
        "organisms": organisms,
        "fetched_at": [first, last],
        "note": note,
    }))
```

`round_floats` を import に足す。

```python
from massbank_context.core.serialization import json_payload, round_floats
```

- [ ] **Step 8: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_tools.py tests/test_pathway_adapter.py -q`
Expected: PASS。

- [ ] **Step 9: 腐敗防止テストを別マークで足す**

`pytest.ini` に足す。

```ini
[pytest]
testpaths = tests
markers =
    network: 実エンドポイントに接続する。既定では走らせない（-m network で実行）
addopts = -m "not network"
```

`tests/test_pathway_upstream.py` を新規作成する。

```python
"""腐敗防止テスト: 実エンドポイントの SPARQL スキーマ追従を検知する。

既定では走らない（pytest.ini の addopts = -m "not network"）。
実行するときは `-m network` を付ける。合成 fixture では上流のスキーマ変更を
検知できないため、これだけは実物を見る（spec §11）。
"""
import pytest

from massbank_context.kg.pathway import PubChemPathwayAdapter

ENDPOINT = "https://rdfportal.org/pubchem/sparql"
# クエン酸。中心代謝なのでパスウェイが必ず付く。
CITRIC_ACID = "KRKNYBCHXYNGOX-UHFFFAOYSA-N"


@pytest.mark.network
def test_real_endpoint_returns_pathways():
    adapter = PubChemPathwayAdapter(ENDPOINT)
    got = adapter.fetch([CITRIC_ACID])
    assert CITRIC_ACID in got, "クエン酸にパスウェイが付かない = スキーマが変わった疑い"
    entry = got[CITRIC_ACID][0]
    assert entry["pathway_id"].startswith("PWID")
    assert entry["pathway_label"]
```

- [ ] **Step 10: 腐敗防止テストが既定で除外されることを確認する**

Run: `C:/Python314/python.exe -m pytest tests -q`
Expected: `test_pathway_upstream.py` が deselect される（全件 PASS）。

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_upstream.py -m network -q`
Expected: PASS（ネットワークがある環境でのみ）。

- [ ] **Step 11: commit**

```bash
git add massbank_context/kg/pathway.py massbank_context/tools/pathway_tools.py pytest.ini tests/test_pathway_adapter.py tests/test_pathway_tools.py tests/test_pathway_upstream.py
git commit -m "feat: PubChem のパスウェイを引くアダプタと resolve_pathways を足す"
```

---

## Task 7: `pathway_activity`

**Files:**
- Modify: `massbank_context/tools/pathway_tools.py`
- Test: `tests/test_pathway_activity_tool.py`

**Interfaces:**
- Consumes: Task 2 の `pathway_activity_rows`、Task 5 のキャッシュ、Task 3 の `incompatible_organism`
- Produces: `pathway_activity(organism: str = "", min_members: int = 3, top_k: int = 30) -> str`

- [ ] **Step 1: 失敗テストを書く**

`tests/test_pathway_activity_tool.py` を新規作成する。

```python
"""pathway_activity ツール（spec §7）。"""
import json

import pytest

from massbank_context.core import config
from massbank_context.core import session as session_module
from massbank_context.kg import cache
from massbank_context.tools.pathway_tools import load_differential, pathway_activity

HEADER = """# contract_version = 1
# group_a = gf_AIN\tn_a = 5
# group_b = gf_HFD\tn_b = 5
# log2fc_sign = positive means group_b is higher
# q_threshold = 0.05\tlog2fc_threshold = 1.0
# n_features_total = 20\tn_with_inchikey = 6\tn_unannotated = 14
spot_id\tname\tname_source\tontology\tinchikey\tinchikey_source\tmsi_level\tmz\trt\tlog2fc\tp_value\tq_value\tmean_a\tmean_b\tsignificant
"""

ROWS = [
    ("1", "IK1", "2.0", "true"), ("2", "IK2", "2.0", "true"),
    ("3", "IK3", "2.0", "true"), ("4", "IK4", "0.1", "false"),
    ("5", "IK5", "0.1", "false"), ("6", "IK6", "0.1", "false"),
]


def _tsv():
    lines = [HEADER.rstrip("\n")]
    for spot, ik, fc, sig in ROWS:
        lines.append("\t".join([spot, f"C{spot}", "arf2", "", ik, "arf2", "",
                                "100.0", "1.0", fc, "0.001", "0.01",
                                "1", "2", sig]))
    return "\n".join(lines) + "\n"


@pytest.fixture(autouse=True)
def prepared(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "PATHWAY_CACHE_PATH",
                        tmp_path / "pathway_cache.sqlite")
    session_module.session.reset()
    path = tmp_path / "differential.tsv"
    path.write_text(_tsv(), encoding="utf-8")
    load_differential(str(path))
    conn = cache.open_cache(config.PATHWAY_CACHE_PATH)
    endpoint = config.PUBCHEM_SPARQL_ENDPOINT
    for ik in ("IK1", "IK2", "IK3"):
        cache.record_result(conn, ik, endpoint, [
            {"pathway_id": "PWID_UP", "pathway_label": "Up pathway",
             "taxid": "TAXID10090"}])
    for ik in ("IK4", "IK5", "IK6"):
        cache.record_result(conn, ik, endpoint, [
            {"pathway_id": "PWID_FLAT", "pathway_label": "Flat pathway",
             "taxid": "TAXID10090"}])
    conn.close()
    yield
    session_module.session.reset()


def test_requires_organism():
    payload = json.loads(pathway_activity())
    assert payload["error"]["code"] == "incompatible_organism"
    assert payload["error"]["available"][0]["taxid"] == "TAXID10090"


def test_rejects_unknown_organism():
    payload = json.loads(pathway_activity(organism="TAXID9999"))
    assert payload["error"]["code"] == "incompatible_organism"


def test_returns_rows_and_header():
    payload = json.loads(pathway_activity(organism="TAXID10090"))
    assert payload["n_background"] == 6
    assert payload["n_no_pathway"] == 0
    assert payload["n_pathways_total"] == 2
    assert payload["group_a"] == "gf_AIN"
    assert "上昇" in payload["log2fc_sign"]
    lines = payload["candidates"].splitlines()
    assert lines[0].split("\t")[0] == "pathway_id"
    assert lines[1].split("\t")[0] == "PWID_UP"
    assert "UP_DOMINANT" in lines[1]


def test_blocks_when_lookup_incomplete():
    conn = cache.open_cache(config.PATHWAY_CACHE_PATH)
    cache.record_failure(conn, "IK1", config.PUBCHEM_SPARQL_ENDPOINT)
    conn.close()
    payload = json.loads(pathway_activity(organism="TAXID10090"))
    assert payload["error"]["code"] == "missing_state"
    assert "resolve_pathways" in payload["error"]["required_tools"]


def test_reports_truncation():
    payload = json.loads(pathway_activity(organism="TAXID10090", top_k=1))
    assert payload["truncated"] is True
    assert payload["n_pathways_total"] == 2
    assert len(payload["candidates"].splitlines()) == 2   # ヘッダ + 1 行
```

- [ ] **Step 2: 失敗を確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_activity_tool.py -q`
Expected: FAIL（`ImportError: cannot import name 'pathway_activity'`）。

- [ ] **Step 3: 実装する**

`massbank_context/tools/pathway_tools.py` に追記する。

```python
from massbank_context.core.errors import (
    contract_mismatch,
    incompatible_organism,
    missing_state,
)
from massbank_context.core.serialization import as_tsv, json_payload, round_floats
from massbank_context.stats.activity import pathway_activity_rows

__all__ = ["load_differential", "resolve_pathways", "pathway_activity"]

_ACTIVITY_COLUMNS = [
    "pathway_id", "pathway_label", "taxid", "n_members", "n_significant",
    "n_up", "n_down", "median_log2fc", "direction", "fisher_p", "fisher_q",
    "tested", "wilcoxon_p", "wilcoxon_tested",
]


def pathway_activity(organism: str = "", min_members: int = 3,
                     top_k: int = 30) -> str:
    """パスウェイ単位の濃縮と方向を返す（spec §7）。

    organism は必須である。生物種を跨いだパスウェイを黙って混ぜず、既定で
    human に寄せることもしない（実測で非哺乳類が 4 割混ざる）。

    背景は「InChIKey が付き、照会が完了し、指定 organism のパスウェイを
    1 つ以上持つ化合物」に限る。パスウェイを持たない化合物を背景に含めると
    全パスウェイが一律に有意寄りになる。除外件数は n_no_pathway で返す。

    並び順は濃縮の q 値だけで決まる。方向の p 値では並べ替えない
    （スコアを合成しない・spec §5.3）。
    """
    if session.differential is None:
        return missing_state(
            "differential", ["load_differential"],
            "差次的結果が読み込まれていません。先に load_differential を実行してください。")

    differential = session.differential
    inchikeys = differential.inchikeys
    conn = cache.open_cache(config.PATHWAY_CACHE_PATH)
    try:
        pending = cache.pending_inchikeys(conn, inchikeys,
                                          config.PUBCHEM_SPARQL_ENDPOINT)
        if pending:
            return missing_state(
                "resolved_pathways", ["resolve_pathways"],
                f"照会が完了していない InChIKey が {len(pending)} 件あります。"
                "部分的な母集団で統計を回すと被覆が変わって濃縮の p 値が動くため、"
                "resolve_pathways を実行してから呼んでください。")
        available = cache.organism_counts(conn, inchikeys)
        if not available:
            return missing_state(
                "resolved_pathways", ["resolve_pathways"],
                "パスウェイが 1 件も引けていません。resolve_pathways を実行してください。")
        if not organism or organism not in {o["taxid"] for o in available}:
            return incompatible_organism(
                available=available,
                message=("organism を指定してください。生物種を跨いだパスウェイを"
                         "混ぜると解釈できません（既定で human に寄せることもしません）。"))
        found = cache.memberships(conn, inchikeys, taxid=organism)
    finally:
        conn.close()

    by_inchikey: dict[str, dict] = {}
    for row in differential.rows:
        inchikey = row["inchikey"]
        current = by_inchikey.get(inchikey)
        # 同じ InChIKey が複数スポットに割れることがある。|log2fc| の大きい方を採る。
        if current is None or abs(row["log2fc"] or 0.0) > abs(current["log2fc"] or 0.0):
            by_inchikey[inchikey] = {
                "inchikey": inchikey,
                "log2fc": row["log2fc"],
                "significant": row["significant"],
            }
    background = [c for ik, c in by_inchikey.items() if ik in found]
    n_no_pathway = len(by_inchikey) - len(background)

    rows, summary = pathway_activity_rows(background, found, min_members)
    truncated = len(rows) > top_k
    shown = rows[:top_k]
    table = as_tsv(_ACTIVITY_COLUMNS,
                   [[round_floats(row.get(column)) for column in _ACTIVITY_COLUMNS]
                    for row in shown])
    return json_payload(round_floats({
        "organism": organism,
        "group_a": differential.group_a,
        "group_b": differential.group_b,
        "log2fc_sign": f"上昇 = {differential.group_b} が高い。",
        "min_members": min_members,
        "n_background": summary["n_background"],
        "n_no_pathway": n_no_pathway,
        "n_significant": summary["n_significant"],
        "n_pathways_total": summary["n_pathways_total"],
        "n_pathways_tested": summary["n_pathways_tested"],
        "truncated": truncated,
        "candidates": table,
        "note": ("direction は検定ではなく記述です。fisher_q と wilcoxon_p は"
                 "別々の量であり、合成しないでください。"),
    }))
```

- [ ] **Step 4: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_pathway_activity_tool.py -q`
Expected: PASS（5 件）。

- [ ] **Step 5: commit**

```bash
git add massbank_context/tools/pathway_tools.py tests/test_pathway_activity_tool.py
git commit -m "feat: パスウェイ単位の濃縮と方向を返す pathway_activity を足す"
```

---

## Task 8: `report/` と登録

**Files:**
- Create: `massbank_context/report/__init__.py`, `massbank_context/report/findings.py`, `massbank_context/tools/report_tools.py`
- Modify: `massbank_context/server.py`
- Test: `tests/test_report.py`, `tests/test_server_registration.py`（新規）

**Interfaces:**
- Consumes: Task 7 の `pathway_activity`（同じ計算を再利用するため、`pathway_activity` の payload を組み立てる内部関数を共有する）
- Produces:
  - `build_findings(facts: dict, rows: list[dict], coverage: float) -> dict`
  - `draft_pathway_findings(organism: str = "", min_members: int = 3, top_k: int = 10) -> str`

- [ ] **Step 1: 骨組みの失敗テストを書く**

`tests/test_report.py` を新規作成する。

```python
"""事実と仮説の分離（spec §10）。"""
import json

from massbank_context.report.findings import (
    INTERPRETATION_WORDS,
    build_findings,
)


def _facts():
    return {
        "comparison": {"group_a": "gf_AIN", "group_b": "gf_HFD",
                       "n_a": "5", "n_b": "5",
                       "log2fc_sign": "正 = gf_HFD が高い（上昇）"},
        "population": {"n_features_total": 20, "n_with_inchikey": 6,
                       "n_unannotated": 14, "n_background": 6,
                       "n_no_pathway": 0},
        "provenance": {"endpoint": "https://rdfportal.org/pubchem/sparql",
                       "fetched_at": ["2026-08-31T00:00:00+00:00", None],
                       "organism": "TAXID10090"},
    }


def _rows():
    return [{"pathway_id": "PWID_UP", "pathway_label": "Up pathway",
             "direction": "UP_DOMINANT", "n_members": 3}]


def test_limits_always_has_five_items():
    findings = build_findings(_facts(), _rows(), coverage=0.5)
    assert len(findings["limits"]) == 5


def test_coverage_is_substituted_into_limits():
    findings = build_findings(_facts(), _rows(), coverage=0.145)
    assert "14.5%" in findings["limits"][0]


def test_hypothesis_slot_is_empty_and_requires_evidence():
    findings = build_findings(_facts(), _rows(), coverage=0.5)
    slot = findings["hypotheses"][0]
    assert slot["claim"] is None
    assert slot["evidence_rows"] == []
    assert "evidence_rows" in slot["required"]


def test_facts_contain_no_interpretation_words():
    """事実側に解釈語が混入しないことを機械的に検査する。

    「活性化」は hypotheses の中にしか置けない。facts に混ざると、
    計算していない主張が計算結果の顔をして出ていく。
    """
    findings = build_findings(_facts(), _rows(), coverage=0.5)
    text = json.dumps(findings["facts"], ensure_ascii=False)
    for word in INTERPRETATION_WORDS:
        assert word not in text


def test_top_pathways_are_carried_verbatim():
    findings = build_findings(_facts(), _rows(), coverage=0.5)
    assert findings["facts"]["top_pathways"] == _rows()
```

- [ ] **Step 2: 失敗を確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_report.py -q`
Expected: FAIL（`ModuleNotFoundError`）。

- [ ] **Step 3: 実装する**

`massbank_context/report/findings.py`:

```python
"""事実と仮説を分けた骨組み（spec §10）。

LLM を呼ばない。空のスロットを返すだけで、埋めるのは LLM である。
「活性化」「不活性化」の語は hypotheses の中にしか置けない。facts 側の語彙は
「上昇（group_b が高い）」「低下」で、群名を必ず伴わせる。
"""
from __future__ import annotations

#: facts 側に現れてはならない語。test で機械的に検査する。
INTERPRETATION_WORDS = ("活性化", "不活性化", "活性が上がった", "活性が下がった")


def _limits(coverage: float) -> list[str]:
    """毎回・固定で出す制約。省略可能にすると、いずれ省略される。"""
    return [
        f"パスウェイ被覆は {coverage * 100:.1f}% です。"
        "パスウェイが付かなかった化合物は「変化が無かった」のではなく「調べていない」ものです。",
        "PubChem のパスウェイは has participant であり、基質か産物か・反応の向き・"
        "量論を持ちません。代謝物の増減からパスウェイの活性は直接は導けません。",
        "化合物同定は MS-DIAL のライブラリ注釈由来で、MS/MS による独立検証は"
        "していません。",
        "メンバー数が閾値未満のパスウェイは検定していません（tested=false）。",
        "異性体は区別していません（MSI Level 2）。",
    ]


def build_findings(facts: dict, rows: list[dict], coverage: float) -> dict:
    """事実・仮説スロット・制約の 3 ブロックを組む。"""
    body = dict(facts)
    body["top_pathways"] = rows
    return {
        "facts": body,
        "hypotheses": [{
            "slot": "pathway_activity",
            "pathway_id": None,
            "claim": None,
            "evidence_rows": [],
            "required": ["pathway_id", "claim", "evidence_rows"],
            "note": ("活性化・不活性化の解釈はここにしか書けません。"
                     "claim には必ず evidence_rows（pathway_id の列）を伴わせてください。"),
        }],
        "limits": _limits(coverage),
    }
```

`massbank_context/report/__init__.py`:

```python
"""事実と仮説を分けた提示の骨組み。"""
from massbank_context.report.findings import INTERPRETATION_WORDS, build_findings

__all__ = ["INTERPRETATION_WORDS", "build_findings"]
```

- [ ] **Step 4: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_report.py -q`
Expected: PASS（5 件）。

- [ ] **Step 5: ツールの失敗テストを書く**

`tests/test_report.py` に追記する。`tests/test_pathway_activity_tool.py` の
`prepared` フィクスチャと同じ前準備をこのファイルにも書く（テスト間の暗黙の結合を作らない）。

```python
import pytest

from massbank_context.core import config
from massbank_context.core import session as session_module
from massbank_context.kg import cache
from massbank_context.tools.pathway_tools import load_differential
from massbank_context.tools.report_tools import draft_pathway_findings

_HEADER = """# contract_version = 1
# group_a = gf_AIN\tn_a = 5
# group_b = gf_HFD\tn_b = 5
# log2fc_sign = positive means group_b is higher
# q_threshold = 0.05\tlog2fc_threshold = 1.0
# n_features_total = 20\tn_with_inchikey = 6\tn_unannotated = 14
spot_id\tname\tname_source\tontology\tinchikey\tinchikey_source\tmsi_level\tmz\trt\tlog2fc\tp_value\tq_value\tmean_a\tmean_b\tsignificant
"""


@pytest.fixture
def prepared(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "PATHWAY_CACHE_PATH",
                        tmp_path / "pathway_cache.sqlite")
    session_module.session.reset()
    lines = [_HEADER.rstrip("\n")]
    spec = [("1", "IK1", "2.0", "true"), ("2", "IK2", "2.0", "true"),
            ("3", "IK3", "2.0", "true"), ("4", "IK4", "0.1", "false"),
            ("5", "IK5", "0.1", "false"), ("6", "IK6", "0.1", "false")]
    for spot, ik, fc, sig in spec:
        lines.append("\t".join([spot, f"C{spot}", "arf2", "", ik, "arf2", "",
                                "100.0", "1.0", fc, "0.001", "0.01",
                                "1", "2", sig]))
    path = tmp_path / "differential.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    load_differential(str(path))
    conn = cache.open_cache(config.PATHWAY_CACHE_PATH)
    for ik in ("IK1", "IK2", "IK3", "IK4", "IK5", "IK6"):
        cache.record_result(conn, ik, config.PUBCHEM_SPARQL_ENDPOINT, [
            {"pathway_id": "PWID_UP", "pathway_label": "Up pathway",
             "taxid": "TAXID10090"}])
    conn.close()
    yield
    session_module.session.reset()


def test_draft_returns_three_blocks(prepared):
    payload = json.loads(draft_pathway_findings(organism="TAXID10090"))
    assert set(payload) == {"facts", "hypotheses", "limits"}
    assert payload["facts"]["comparison"]["group_b"] == "gf_HFD"
    assert payload["facts"]["provenance"]["organism"] == "TAXID10090"
    assert len(payload["limits"]) == 5


def test_draft_propagates_errors(prepared):
    payload = json.loads(draft_pathway_findings())
    assert payload["error"]["code"] == "incompatible_organism"
```

- [ ] **Step 6: `draft_pathway_findings` を実装する**

`massbank_context/tools/report_tools.py`:

```python
"""提示層のツール（spec §10）。"""
from __future__ import annotations

import json

from massbank_context.core.serialization import json_payload, round_floats
from massbank_context.core.session import session
from massbank_context.report.findings import build_findings
from massbank_context.tools.pathway_tools import pathway_activity

__all__ = ["draft_pathway_findings"]


def draft_pathway_findings(organism: str = "", min_members: int = 3,
                           top_k: int = 10) -> str:
    """事実・仮説スロット・制約の骨組みを返す（spec §10）。

    pathway_activity と同じ計算を使い、結果を事実として写す。解釈は書かない。
    「活性化」「不活性化」の語は hypotheses の中にしか置けない。
    """
    raw = pathway_activity(organism=organism, min_members=min_members,
                           top_k=top_k)
    payload = json.loads(raw)
    if "error" in payload:
        return raw

    differential = session.differential
    columns = payload["candidates"].splitlines()[0].split("\t")
    rows = [dict(zip(columns, line.split("\t")))
            for line in payload["candidates"].splitlines()[1:]]
    n_resolved = payload["n_background"] + payload["n_no_pathway"]
    coverage = (payload["n_background"] / n_resolved) if n_resolved else 0.0

    facts = {
        "comparison": {
            "group_a": payload["group_a"],
            "group_b": payload["group_b"],
            "n_a": differential.meta.get("n_a"),
            "n_b": differential.meta.get("n_b"),
            "log2fc_sign": payload["log2fc_sign"],
        },
        "population": {
            "n_features_total": differential.n_features_total,
            "n_with_inchikey": len(differential.inchikeys),
            "n_unannotated": differential.n_unannotated,
            "n_background": payload["n_background"],
            "n_no_pathway": payload["n_no_pathway"],
            "n_pathways_total": payload["n_pathways_total"],
            "n_pathways_tested": payload["n_pathways_tested"],
        },
        "provenance": {
            "source": session.differential_source,
            "organism": payload["organism"],
        },
    }
    return json_payload(round_floats(build_findings(facts, rows, coverage)))
```

- [ ] **Step 7: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_report.py -q`
Expected: PASS（7 件）。

- [ ] **Step 8: server に登録する**

`massbank_context/server.py` の import に足す。

```python
from massbank_context.tools.pathway_tools import (
    load_differential as _load_differential,
    pathway_activity as _pathway_activity,
    resolve_pathways as _resolve_pathways,
)
from massbank_context.tools.report_tools import (
    draft_pathway_findings as _draft_pathway_findings,
)
```

`mcp.tool(...)` の並びの末尾に足す。

```python
mcp.tool(structured_output=False)(_load_differential)
mcp.tool(structured_output=False)(_resolve_pathways)
mcp.tool(structured_output=False)(_pathway_activity)
mcp.tool(structured_output=False)(_draft_pathway_findings)
```

- [ ] **Step 9: 登録テストを書いて確認する**

`tests/test_server_registration.py` を新規作成する。

```python
"""MCP 登録面を固定する（登録漏れの早期検出）。"""
import asyncio

from massbank_context import server

EXPECTED_TOOLS = sorted([
    "draft_pathway_findings",
    "index_status",
    "inspect_spectrum",
    "load_differential",
    "load_spectra",
    "match_spectra",
    "pathway_activity",
    "resolve_pathways",
])


def test_registered_tools():
    tools = asyncio.run(server.mcp.list_tools())
    assert sorted(t.name for t in tools) == EXPECTED_TOOLS


def test_no_structured_output():
    """structured_output=False が全ツールに効いていること（content の二重送出を防ぐ）。"""
    tools = asyncio.run(server.mcp.list_tools())
    assert all(getattr(t, "outputSchema", None) is None for t in tools)
```

Run: `C:/Python314/python.exe -m pytest tests/test_server_registration.py -q`
Expected: PASS（2 件）。

- [ ] **Step 10: 全体テストを実行する**

Run: `C:/Python314/python.exe -m pytest tests -q`
Expected: 全件 PASS（既存 198 件 + 本計画の新規約 60 件）。腐敗防止テストは deselect される。

- [ ] **Step 11: CLAUDE.md を更新する**

`CLAUDE.md` の「現状は Phase 1」の記述を差し替える。

```markdown
現状は Phase 1（索引ビルド + バイアス無し照合）に加え、パスウェイ活性同定の
Phase A〜D（`experiment/` `kg/` `stats/` `report/`、ツール 4 本）まで。
設計は [docs/superpowers/specs/2026-08-31-pathway-activity-design.md](docs/superpowers/specs/2026-08-31-pathway-activity-design.md)。
`bias/` `catalog/` と spec §7.2 のツール 10 本のうち残りは未実装。
```

アーキテクチャの図にも 4 モジュールを足す。

```
core/       設定・セッション状態・データモデル・層の導出・ビン化・直列化・エラー契約
stats/      Fisher / BH-FDR / Wilcoxon / パスウェイ集計（純関数・leaf）
ttl/        MassBank TTL リーダ（世代判定 → gen_a / gen_b）
query/      実験スペクトルのリーダ（.msp / .mgf / .dcl）
experiment/ 差次的結果の契約リーダ（ms-data-parser との唯一の接点）
kg/         InChIKey → パスウェイ（PubChem SPARQL + SQLite キャッシュ）
index/      索引ビルド CLI・SQLite スキーマ・IDF 背景・質量誤差プロファイル
spectral/   照合エンジン
report/     事実と仮説を分けた骨組み
tools/      MCP 公開層（実体はここ、server.py は登録するだけ）
server.py   FastMCP ファサード
```

- [ ] **Step 12: commit**

```bash
git add massbank_context/report massbank_context/tools/report_tools.py massbank_context/server.py CLAUDE.md tests/test_report.py tests/test_server_registration.py
git commit -m "feat: 事実と仮説を分けた骨組みを返し、4 ツールを登録する"
```

---

## 完了条件

- [ ] `C:/Python314/python.exe -m pytest tests -q` が全件 PASS（腐敗防止テストは deselect）
- [ ] `C:/Python314/python.exe -m pytest tests/test_pathway_upstream.py -m network -q` が PASS（ネットワークがある環境で 1 回は確認する）
- [ ] MCP ツールが 8 本登録されている
- [ ] `requirements.txt` が変わっていない（依存を増やしていない）
- [ ] `facts` 側に「活性化」の語が現れないことがテストで守られている
