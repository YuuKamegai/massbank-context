# massbank-context 試用・検証プラン

**作成**: 2026-08-31
**対象**: Phase 1 の 4 ツール（`index_status` / `load_spectra` / `inspect_spectrum` /
`match_spectra`）と索引ビルド CLI
**前提**: ビン幅の推定器を直した直後（0.0481 → 0.0075）の索引で評価する

## 目的

「MCP サーバとして会話から使える」ことと、「返す候補が信用に足る」ことを分けて
確かめる。前者は接続の問題、後者は照合の問題であり、混ぜると原因が切り分けられない。

## Step 0 — MCP 連携（完了済み）

`.mcp.json` をリポジトリ直下に置いた。stdio で `initialize` → `list_tools` →
`call_tool` まで実際に往復させ、次を確認済み:

- サーバ名 `massbank-context`、ツール 4 本が列挙される
- `index_status` が `bin_width 0.0075` / `n_records 84469` / `n_inchikey 8091` /
  `acquisition_classes ["GC-EI-MS","LC-ESI-MS2","MALDI-MS2"]` を返す
- `structuredContent` が付かない（spec §12 の `structured_output=False` が効いている）

新しいセッションでは、プロジェクトスコープの `.mcp.json` の使用許可を一度求められる。

## Step 1 — 既知解テスト（自己照合）

**正解の分かっている問いでしか照合器の健全性は測れない。** 索引の元になった
MassBank レコードそのものをクエリとして投げ、自分が 1 位に来るかを見る。

題材は `scripts/export_selfcheck_msp.py` が作る:

```
python scripts/export_selfcheck_msp.py \
    --ttl-dir C:/Users/yuu18/massbank-context/MassBank-RDF/draft_ttl/MassBank_RDF_ttl \
    --output C:/Users/yuu18/datasets/_selfcheck/massbank_selfcheck \
    --per-class 5 --stride 400
```

生成物は `massbank_selfcheck.msp`（20 件 / 4 装置クラス）と、
`massbank_selfcheck.answer.tsv`（正解表）。

### 合格条件

| 見るもの | 合格 |
|---|---|
| 候補が付いたスペクトル | 20/20 |
| 1 位が**同一 InChIKey** | 19/20 以上 |
| 1 位が同一の骨格（InChIKey 第 1 ブロック） | 20/20 |
| `weighted_score` | 1.0 付近 |
| `band` | すべて `PASS` |

### 読み方の注意

- **1 位の accession が「自分」でなくても不合格ではない。** 同一化合物の CE 違いが
  索引に多数あり、`weighted_score` が 1.0 で並ぶ。同点は accession で決定的に
  並べるため、辞書順で前のレコードが 1 位に来る。判定は InChIKey 単位で行う
- 骨格が同じで立体表記だけ違う（InChIKey 第 1 ブロックが同じで第 2 ブロックが違う）
  ケースは**設計どおりの限界**であり、MSI Level 2 の意味そのもの。不合格にしない
- **実施済み（2026-08-31）: 合格。**
  [Step 1 実施結果](2026-08-31-selfcheck-step1-results.md) に採点表・不一致行・
  同点が別 InChIKey に割れた 3 件を残した。採点には正解表の読み取りが要るため
  MCP セッションからは完結せず、ローカル実行が必要である
- 事前実行（2026-08-31）の実績: 候補 20/20、1 位が同一 InChIKey 19/20、
  残り 1 件も骨格一致（`AMRQXHFXNZFDCH-SECBINFHSA-N` に対し
  `AMRQXHFXNZFDCH-UHFFFAOYSA-N` が 1 位）。全件 score 1.0 / `PASS`

### 落ちたときの読み方

ここが落ちたら実データに進む意味がない。疑う順序は
(1) ビン化（`bin_width` と `probe_bins` の窓）、(2) IDF の数え方、
(3) 1 対 1 一致の制約、(4) `.msp` 書き出し側の不備。

## Step 2 — 実データ 1 ファイル

`C:/Users/yuu18/datasets/20230824_liver/20230824_liver/POS/20230823_liver_gf_AIN_1_POS_202308241734.dcl`

### 手順

1. `load_spectra` — `n_spectra` と `n_skipped_empty` を確認する
2. `match_spectra` は **`spectrum_ids` で 10〜20 本に絞って**呼ぶ
   （既定は読み込んだ全件で、候補行 200 の上限に当たり `truncated` になる）

### 見る数字

`matched_fragments` / `total_query_peaks` / `weighted_score` / `cosine` /
`unexplained_peak_ratio` / `band` / `precursor_delta`、および payload の
`fragment_tolerance`（現在 ±0.0112）。

### 既知の実績（2026-08-31、参照値）

- `n_spectra 917` / `n_skipped_empty 6573` / ピーク中央値 91
- 全 917 本を照合すると候補が付くのは 118 本（12.9%）、
  `band` は PASS 100 / WEAK 12 / PRECURSOR_ONLY 6

ヒット率が 1 割強なのは、索引が高分解能 8,091 化合物であるのに対し肝臓
リピドミクスの主成分（脂質）が MassBank に薄いため。**索引の欠落であって
照合器の失敗ではない。**

### 危険な読み方

`band` が `PRECURSOR_ONLY` の候補は精密質量しか合っていない。`PASS` と
同一視しない。`weighted_score` は IDF 重み付きカバレッジで、cosine とは
別物であり、合成してはならない（spec §5.3）。

## Step 3 — MS-DIAL の注釈との突き合わせ

massbank-context 単体では「合っているか」が分からないので、外部の答え合わせを使う。

1. `ms-data-parser` の `pai2_parser` / `verify_peak_annotation` で、MS-DIAL が
   注釈を付けたピークを取る（`analytical_checks.msms.band` が `PASS` のもの、
   すなわち実スペクトルが見えているものに限る）
2. 同じ precursor のスペクトルを massbank-context に投げ、第 1 候補と突き合わせる

### 不一致の内訳を必ず分ける

1. MassBank に該当化合物が無い（索引の欠落）
2. 異性体で分けられない（MSI Level 2 の限界）
3. 本当の誤り（照合器の問題）

この 3 つを混ぜて「一致率 X%」とだけ言うと、改善すべき場所が分からなくなる。

## Step 4 — 失敗経路の確認

黙って壊れないことの確認。設計上わざと拒否する経路である。

| 操作 | 期待 |
|---|---|
| `load_spectra` せずに `match_spectra` | `missing_state`（次に呼ぶべきツール名付き） |
| `acquisition_class="LC-ESI-MS3"` | `incompatible_index`（別層で黙って代用しない） |
| `spectrum_ids` に存在しない id | `unknown_spectrum_ids` に列挙、`n_spectra` に数えない |
| 全件で `match_spectra` | `truncated: true` と打ち切り理由 |

エラーは `isError` ではなく本文 JSON のエンベロープで返る（使用中の MCP SDK が
`isError=true` の経路で structuredContent を捨てるため）。

## Step 5 — 記録

Step 1〜4 の結果と、見つかった不足を記録する。特に次を残す:

- 「候補が付いたが明らかに違う」ケース — 索引の層・ビン幅・IDF のどれを
  疑うべきかの判断材料になる
- 打ち切り・スキップの件数（`truncated` / `n_skipped_empty` /
  `unknown_spectrum_ids`）が実データでどう出たか

## この検証で分からないこと

- **群間比較**（gf_AIN vs gf_HFD）— massbank-context の外。`ms-data-parser` の
  `arf_differential` の領分
- **疾患・生理活性との関連** — `kg/` 未実装。データの出所はラボ内 private Virtuoso 待ち
- **異性体の区別** — MSI Level 2 止まり
