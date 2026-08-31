# 到達度調査 — 群間比較からパスウェイ変化まで、今どこまで届くか

**調査日**: 2026-08-31
**問い**: 上流クローンの移動とパス修正、および
[上流比較調査](2026-08-31-massbank-rdf-upstream-comparison.md)を終えた時点で、
**「実験データの群間比較 → 生化学的パスウェイの変化（活性・不活性）の同定」**
という目的に対してどこまで到達しているか。

## 結論

鎖は **「群間比較」と「スペクトル同定」の 2 本が完成**していて、
**その 2 本を繋ぐ環と、化合物 → パスウェイ → 変化の方向、の環が無い。**

ただし本調査で、欠けている環のうち 1 つ（パスウェイの取得元）が
**公開エンドポイントで解けることを実測した**。これは
「Phase 2 の KG は全てラボ内 private Virtuoso 待ち」という従来の理解を、
**疾患・活性に限った話へ縮める**。

## 1. 現に動くもの（実測）

| 確認 | 結果 |
|---|---|
| `pytest tests -q` | 198 passed |
| MCP `index_status`（移動後の実呼び出し） | 応答あり。`.mcp.json` の絶対パスは新配置と整合 |
| 索引 `index/massbank.sqlite` | 433 MB / 84,469 レコード / **8,091 InChIKey** / 4,467,236 ピーク / bin_width 0.0075 |
| 層の内訳 | LC-ESI-MS2 high POS 56,812・NEG 26,516 / GC-EI-MS high 1,085 / MALDI-MS2 56 |
| 化合物クラス | `compound_class` 31,422 行（ClassyFire 4 階層 22,522 行 + massbank 由来 8,900 行） |
| 既知解テスト（Step 1） | 合格。1 位が同一 InChIKey 19/20、同一骨格 20/20 — [結果](../verification/2026-08-31-selfcheck-step1-results.md) |

MCP ツールは spec §7.2 の 10 本中 **4 本**。`kg/` `bias/` `catalog/` `report/` は
ディレクトリごと未実装である。

## 2. 「上流が更新された」の実態

検証を保留する理由として「MassBank RDF の最新版が別リポジトリで更新されていた」が
挙がったが、これは[上流比較調査](2026-08-31-massbank-rdf-upstream-comparison.md)の
結論と食い違う。`MotohiroOGAWA/MassbankRdf` は fork ではなく**共通祖先を持たない
別履歴のアプリ**で、TTL を 1 件も含まない（11.7 万件の TTL・TXT→TTL 変換・QLever 設定は
いずれも `systemsomicslab/MassBank-RDF` 側にある）。

**したがって索引の元データは陳腐化していない。** 検証プランの Step 1〜4 は
そのまま実行できる（実際に Step 1 は本日実行して合格した）。

保留に意味があるとすれば「照合の後段（KG・統計・GUI）をどう作るかの設計判断待ち」で
あって、データの鮮度ではない。

## 3. 目的への到達度（環ごと）

| # | 環 | 状態 | 実体 |
|---|---|---|---|
| 1 | 実験データ読込（`.arf` / `.arf2` / `.pai2` / `.dcl` / `.eic`） | 済 | ms-data-parser |
| 2 | 群間比較 | 済 | `arf_preprocess` → `arf_differential`（Welch + log2FC + FDR）、PCA・volcano。**2 群のみ**（多群 ANOVA は意図的に封鎖） |
| 3 | 差次的フィーチャ → InChIKey | 部分的 | `.arf2` / `.pai2` は InChIKey を持つ（`arf2/reader.py` Key16、`pai2/reader.py` data[29]）が、差次的解析が回る `.arf` 側は `annotation`（名前）だけ。**spot id での結合が要る** |
| 4 | MS/MS で同定を検証 | 済（被覆は薄い） | massbank-context Phase 1。肝リピドミクス実測でヒット率 12.9%（917 本中 118 本）。索引の欠落であって照合器の失敗ではない |
| 5 | 群比較の結果を massbank-context へ渡す口 | **無い** | massbank-context には群・試料の概念が一切ない。`load_spectra` はファイル単位 |
| 6 | 化合物 → パスウェイ | 到達可・**未実装** | §4 を参照 |
| 7 | パスウェイ → 変化の方向（活性・不活性） | **無い** | 設計もされていない。上流アプリにあるのは疾患単位の Fisher + BH であって、パスウェイ単位でも方向でもない |
| 8 | 化合物 → 疾患・生理活性 | ネットワーク的に不可 | HMDB・KNApSAcK の既定 endpoint は `http://10.242.145.97:3030/...`（private IP） |

## 4. 環 6 の実測 — パスウェイは公開ソースで引ける

`https://rdfportal.org/pubchem/sparql` は**認証なしで到達でき、
InChIKey → CID → Pathway が実際に返る**。索引内の化合物での実結果:

```
KRKNYBCHXYNGOX-UHFFFAOYSA-N (Citric acid)  | Citric Acid Cycle
KRKNYBCHXYNGOX-UHFFFAOYSA-N                | Propanoate Metabolism
ODBLHEXUDAPZAU-UHFFFAOYSA-N (Isocitrate)   | Congenital Lactic Acidosis
```

クエリの形は `MassbankRdf` の
`massbank_rdf/services/kg/query_builders/pubchem_query_builder.py`
`build_pubchem_pathway_query`（`FROM <http://rdf.ncbi.nlm.nih.gov/pubchem/pathway>` /
`?pathway a vocab:Pathway ; obo:RO_0000057 ?cid ; dcterms:title ?label`）が
そのまま使える。

**被覆**: 索引の 8,091 InChIKey から 200 件を無作為抽出（seed 0）して
50 件ずつ 4 バッチで照会したところ、パスウェイが付いたのは **29 件 = 14.5%**
（標準誤差 約 2.5%）。中心代謝は当たり、脂質はほぼ当たらない。SMPDB 由来のため
疾患名のパスウェイ（`Congenital Lactic Acidosis` 等）も混じる。

**この被覆は payload に必ず出す必要がある。** 出さないと
「パスウェイが付かない ＝ 変化が無い」と読まれる。85% は「調べていない」である。

一方、疾患・生理活性は公開ソースに無い。HMDB は公開 SPARQL 無し（XML 数 GB の
ダウンロードのみ）、KNApSAcK は公開 API 無し、`MassBank-RDF` の QLever には
MassBank しか載っていない。ここはラボへの確認待ちのまま変わらない。

## 5. 上流アプリから取ってはいけないもの

`MassbankRdf` の照合 SQL は query peak ごとに最近傍 reference peak を選ぶだけで、
**同じ reference peak の再利用を禁じていない**（多対一 → dot product の水増し →
cosine が 1 を超える）。massbank-context で既に潰した欠陥（`969f19a`）なので、
移植対象から外す。

移植する価値があるのは次の 4 つ:

- PubChem / HMDB / KNApSAcK の SPARQL クエリビルダ
- short InChIKey（先頭 14 文字）検索
- KG evidence JSON（LLM 向けに圧縮した証拠の形）
- チェックポイント付きバッチ処理と ZIP 出力

## 6. 次の一手

**環 5（群 → 候補化合物リストの受け渡し）を先に作る**ことを勧める。
理由は、環 6 が公開ソースで解けると分かった今、鎖の中で**唯一「誰も持っていない」
のがこの受け渡し**だからである。環 3・4・6 はいずれも既存部品の接続で済む。

具体的には:

1. `kg/` に PubChem パスウェイのアダプタを置く（`SparqlAdapter` の最初の実装）。
   diseases / activities は spec §8 の `upstream_unavailable` で穴を開けておく
2. InChIKey 群を受けてパスウェイ単位に畳む口を足す。**被覆 14.5% を必ず併記する**
3. 方向（活性・不活性）は、パスウェイ構成メタボライトの log2FC の符号を
   集約する設計になる。仕様が無いので設計から起こす

## 参照

- [上流リポジトリ比較調査](2026-08-31-massbank-rdf-upstream-comparison.md)
- [試用・検証プラン](../verification/2026-08-31-massbank-context-trial-plan.md)
- [Step 1 実施結果](../verification/2026-08-31-selfcheck-step1-results.md)
- 設計書 [spec](../superpowers/specs/2026-08-28-massbank-context-mcp-design.md) §5.2（kg/ の 3 アダプタ）・§7.2（ツール 10 本）・§9（フェーズ）・§11（未確定事項）
