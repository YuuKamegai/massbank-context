# 試用・検証 Step 1 実施結果（既知解テスト）

**実施日**: 2026-08-31
**対象**: [試用・検証プラン](2026-08-31-massbank-context-trial-plan.md) の Step 1
**索引**: `index/massbank.sqlite`（2026-08-30 15:12Z / bin_width 0.0075 /
84,469 レコード / 8,091 InChIKey / `resolution_class == "high"`）
**題材**: `C:/Users/yuu18/datasets/_selfcheck/massbank_selfcheck.msp`（20 件）と
`massbank_selfcheck.answer.tsv`（正解表）

## 実施条件

`load_spectra(msp)` → `match_spectra(top_k=3)`（`acquisition_class` は既定の
`LC-ESI-MS2`、`precursor_tolerance` 0.01、実効 `fragment_tolerance` 0.01125）。

**注意**: MCP セッション（会話）からは正解表を読めない。massbank-context も
ms-data-parser も汎用のファイル読み取りツールを公開していないためで、これは
設計どおりの制約である。**既知解テストの採点はローカル実行を要する。** 本結果は
同一のツール関数をローカルで直接呼んで採点した。候補側の数値は MCP 経由の実行と
一致しており、経路による差は無い。

## 合格判定

| 見るもの | 合格条件 | 実測 | 判定 |
|---|---|---|---|
| 候補が付いたスペクトル | 20/20 | **20/20** | 合格 |
| 1 位が同一 InChIKey | 19/20 以上 | **19/20** | 合格 |
| 1 位が同一骨格（第 1 ブロック） | 20/20 | **20/20** | 合格 |
| `weighted_score` | 1.0 付近 | **1 位は全 20 件が 1.0** | 合格 |
| `band` | すべて PASS | **1 位は全 20 件が PASS** | 合格 |

補助的な数値:

- `n_spectra` 20 / `n_skipped_empty` 0 / 候補総数 60
- 正解 InChIKey が top-3 に入る: **20/20**
- 1 位が同一 accession（＝自分自身）: 15/20。**これは不合格ではない**（下記）
- `precursor_delta` は 19 件が 0.0、`massbank_selfcheck:6` のみ 0.0001
- band の内訳（候補 60 件全体）: PASS 59 / WEAK 1 / PRECURSOR_ONLY 0
- 候補全体の `weighted_score` 最小は 0.014（`:19` の 3 位、1/108 一致）

## 不一致行（1 件）

| spectrum_id | 正解 | 1 位 | 性質 |
|---|---|---|---|
| `massbank_selfcheck:5` | `MSBNK-UFZ-WANA010601AD6CPH` / `AMRQXHFXNZFDCH-SECBINFHSA-N` | `MSBNK-Eawag-EA013602` / `AMRQXHFXNZFDCH-UHFFFAOYSA-N` | 骨格同一・立体層のみ相違 |

正解レコード自身は 3 位に `weighted_score` 1.0・`cosine` 1.0 で入っている。
1〜3 位が全て 1.0 で並び、accession の辞書順で `Eawag` が先に出た形である。
**MSI Level 2 の設計上の限界**であり、照合器の欠陥ではない。

## 同点が別 InChIKey に割れた 3 件

1 位と同点（`weighted_score` 1.0）で InChIKey が異なる候補が並んだのは 3 件。
**うち 2 件は立体層のみの相違だが、1 件は骨格が別である。**

| spectrum_id | 1 位 | 同点の対抗馬 | 正解 |
|---|---|---|---|
| `:3` | `QNFBKOHHLAWWTC-UHFFFAOYSA-N`（cos 1.0） | `HOEVRHHMDJKUMZ-UHFFFAOYSA-N`（Isofraxidin, cos 0.9753） | 1 位が正解 |
| `:5` | `AMRQXHFXNZFDCH-UHFFFAOYSA-N` | `AMRQXHFXNZFDCH-SECBINFHSA-N` | 対抗馬が正解 |
| `:9` | `UFEODZBUAFNAEU-NLRVBDNBSA-N` | `UFEODZBUAFNAEU-UHFFFAOYSA-N` | 1 位が正解 |

`:3` は両者とも分子式 `C11H10O5` の位置異性体で、15 本のピークが全て一致し
`weighted_score` が同値になった。**1 位が正解になったのは accession の辞書順が
たまたま味方しただけで、設計が保証したものではない**（`MSBNK-N…` < `MSBNK-R…`）。
`cosine` なら 1.0 対 0.9753 で正しく割れるが、cosine は並び替えに使わないと
決めている（spec §6.8 / §5.3、スコアを合成しない）。ここは合格条件を満たしたまま
残る既知の脆さであり、**「同点内では化合物を決めていない」ことを候補提示側で
明示する**のが筋である（Phase 2 以降の課題）。

## 副次的に見つかったデータ品質の問題

索引の `compounds.exact_mass` で `HOEVRHHMDJKUMZ-UHFFFAOYSA-N`（Isofraxidin,
`C11H10O5`）が **222.196** になっている。単同位体質量は 222.0528 で、同じ分子式の
`QNFBKOHHLAWWTC-UHFFFAOYSA-N` は 222.0528 と正しい。上流 TTL 側の値であり
massbank-context の計算ではないが、**将来 exact_mass で絞り込むと静かに落ちる**。
精密質量による絞り込みを実装するときは、分子式から再計算した値と照合すること。

## 結論

Step 1 は合格。照合器の健全性は確認できたので、Step 2（実データ 1 ファイル）へ
進んでよい。持ち越す課題は次の 2 点である。

1. 同点かつ骨格が別のときに「決められていない」ことを payload で明示する
2. `compounds.exact_mass` を信用しない（上流の誤値が実在する）
