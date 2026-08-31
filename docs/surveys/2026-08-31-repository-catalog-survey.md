# repository_catalog 調査報告

**対象**: `C:/Users/yuu18/repository_catalog/native-smoke.sqlite`
**調査日**: 2026-08-31（DB のクロールは 2026-08-27〜28）
**目的**: spec §7.2 のツール 7 `sample_context(inchikeys[])` が前提にしている
`msdial_repository_catalog` として、この DB が使えるかを判断する。

## 結論

**「化合物 → 検出実績」の方向には使えない。「文脈 → 研究・試料」の方向と、
語彙の供給には今すぐ使える。**

理由は単純で、InChIKey を保持する `metabolite_entity` と、それを研究に結び付ける
`metabolite_observation` が**どちらも 0 行**だから。massbank-context が候補として
返すのは InChIKey なので、結合キーの側にデータが無い。技術ではなくデータの
充填状況の問題であり、spec §7.2 のツール 7 は定義どおりには実装できない。

## 調査方法

読み取り専用・不変で開いて、テーブルの行数と代表値を数えただけ。書き込みは
一切していない。

```python
sqlite3.connect(
    "file:C:/Users/yuu18/repository_catalog/native-smoke.sqlite?mode=ro&immutable=1",
    uri=True)
```

`-shm` が残っているため、書き込みプロセスと競合させないこと。`immutable=1` は
WAL を無視するので、クロール中に読むと取りこぼす（今回は `-wal` が 0 バイト）。

## 実測

DB サイズ 4.0 GB、`schema_info.version = 2`。

### クロールの範囲

| リポジトリ | study 数 | ライセンス |
|---|---|---|
| metabolomics_workbench | 4,508 | CC BY 4.0 |
| metabolights | 875 | EMBL-EBI Terms of Use |
| metabobank | 114 | （空） |
| mb_post | 104 | CC0 1.0 |
| **計** | **5,601** | |

`crawl_run` は 17 件。最新は 2026-08-28、metabolights のみ
`completed_with_errors`（discovered 1,097 / hydrated 806）。

### 中身のあるテーブル

| テーブル | 行数 | 備考 |
|---|---|---|
| `sample_attribute` | 12,061,500 | 元データの属性をそのまま保持 |
| `sample_context` | 2,600,998 | 正規化済みの生物文脈 |
| `sample` | 1,051,091 | |
| `raw_file` | 365,276 | `download_url` 付き。raw 283,395 / converted 72,642 |
| `analysis_unit` | 9,855 | LC-MS 8,035 / GC-MS 1,161 / DI-MS 276 |
| `publication` | 6,169 | DOI・PubMed ID |
| `source_snapshot` / `source_blob` | 5,608 / 5,603 | 原データ 1.55 GB を 64.7 MB に圧縮保持 |
| `study` | 5,601 | `study_fts` で全文検索可 |
| `analysis_unit_evidence` | 29,280 | 分類の根拠 |
| `analysis_unit_context` | 3,979 | すべて `study_topic` |

`analysis_unit.ion_mode` は Positive 5,035 / Negative 3,647 / Both 717 /
Unknown 456。装置は Thermo Q Exactive 系が最多（Orbitrap 系で約 2,700）。
`review_status` は needs_review 6,837 / unreviewed 3,018 で、**人手レビュー済みの
ユニットは無い**。

### 空のテーブル（0 行）

`metabolite_entity` / `metabolite_observation` / `metabolite_response` /
`analysis_run` / `contrast` / `class_proposal` / `class_assignment` /
`manual_override`。

つまり「どの研究で・どの化合物が・どちら向きに動いたか」という層は、
スキーマとしては用意されているが**まだ空**。これを埋めるには `raw_file` を
MS-DIAL で再処理して `analysis_run` → `metabolite_observation` を作る必要があり、
それは catalog 側パイプラインの仕事で、massbank-context の責務ではない。

### `sample_context` の語彙

| category | 行数 | 異なり数 |
|---|---|---|
| organ | 910,624 | 2,874 |
| species | 419,127 | 1,650 |
| intervention | 233,240 | 5,940 |
| sample_role | 205,620 | 892 |
| sex | 191,303 | 69 |
| condition | 157,884 | 1,994 |
| time | 147,597 | 2,126 |
| genotype | 104,748 | 4,372 |
| disease | 92,382 | 512 |
| aging | 73,893 | 518 |
| cell_type | 12,427 | 374 |

**正規化は粗い。** 実際に別値として並んでいるもの:

- disease: `Pre-MetS` / `MetS` / `CONTROL` / `CASE` / `Healthy Control` / `Normal`
  — 疾患名と対照群ラベルが同じ軸に混在
- organ: `Serum` / `Blood (serum)` / `Plasma` / `Blood (plasma)` が別値
- species: `Homo sapiens` / `Homo sapiens (Human)` が別値。さらに `urine` /
  `blood plasma` / `blood serum` という**検体名が species に混入**

そのまま LLM に選択肢として見せてはならない。出現数の閾値か正規化を挟む必要がある。

## massbank-context との連携

### 使えない: spec §7.2 ツール 7 の定義どおりの実装

`sample_context(inchikeys[])` が返すはずの「catalog 由来の検出実績（accession
付き）」は、`metabolite_observation` が空である限り作れない。
**spec §8 の `upstream_unavailable` を返す形で穴を開けておくのが正しい**
（`kg/` で既に採っている方針と同じ）。

代替として `study_fts` による化合物名の全文検索は可能だが（例: `fatty liver`
で 51 study）、名前ベースで InChIKey に落ちず、MSI レベルも不明なので
「検出実績」とは呼べない。**この代替を検出実績として提示してはならない。**

### 使える 1: `set_objective` の実在語彙

spec §7.3 は「LLM に自由なカテゴリを作らせず、実在語彙から選ばせる」設計で、
その語彙の取得元が §11 の未解決事項（ラボ内 private Virtuoso 待ち）だった。
`sample_context` の disease / organ / species / intervention は実測語彙なので、
**private Virtuoso を待たずに Phase 3 の前提を満たせる経路**になる。
上記のとおり正規化が要る。

### 使える 2: 文脈側からの照会

「文脈 → 研究・試料」の方向は今すぐ引ける。例: organ に `liver` を含む study は
308 件。`raw_file.download_url` があるので、仮説の裏取りに使う公開データの
当たりを付ける用途にも使える。

### 実装方針

`catalog/` を**読み取り専用アダプタ**として追加する（spec §5 の依存方向どおり、
`catalog/` は上位を import しない）。DB 本体は索引・TTL と同様にリポジトリ外へ
置いたままにする（`.gitignore` が `*.sqlite` を見ているのと方針が一致し、
再配布可能なデータを成果物に焼き込まない spec §4 とも整合する）。

## 未確認事項

1. `metabolite_*` を埋める予定と時期（catalog 側パイプラインの状況）
2. `native-smoke` という名前が示すとおり、これがスモークテスト用の生成物なのか、
   本番の成果物なのか
3. metabolights の `completed_with_errors`（1,097 中 806 のみ hydrated）で
   何が落ちたか
4. ライセンス表記が空の 116 study の扱い（再配布可否）
