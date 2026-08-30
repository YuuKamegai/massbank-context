# massbank-context MCP サーバ 設計書

作成日: 2026-08-28

> **この設計書は Lipidmix_with_LLM とは別のシステムを記述する。**
> 本リポジトリの `docs/superpowers/specs/` を設計書置き場として借りているだけで、
> 実装は別リポジトリ（`massbank-context`）に置く。Lipidmix 側からは参照しない。
> 継承するのは規約（payload 規律・`missing_state` 封筒・テスト方針）だけである。

## 1. 目的

ユーザーが実験で得た MS/MS スペクトルを起点に、MassBank-RDF と照合して候補化合物を得、
そこから病態・生理活性・パスウェイのメタデータを取得し、事実と仮説を分離して提示する。
この一連を LLM が MCP 経由で実行できるようにする。

中心的な主張は次の一文である。

> プロダクトイオンのパターンとメタデータのカテゴリには連関がある。
> その連関を**背景分布に対する濃縮度として定量化**し、
> 候補取得の**説明可能なバイアス**として使う。

「説明可能」が要件である。バイアスが結果を何位分動かしたかが常に数値で見えなければ、
この枠組みは確証バイアスで偽同定を生む装置になる。

## 2. スコープ

### やること

- 実験スペクトル（`.msp` / `.mgf` / MS-DIAL `.dcl`）の取り込み
- MassBank-RDF 由来のローカル索引に対する照合と、一致フラグメントの列挙
- InChIKey を鍵としたメタデータ取得（疾患・活性・パスウェイ・生体試料・生物種）
- 目的（自然文）から実在語彙のカテゴリへの写像
- カテゴリとフラグメントの濃縮度統計と、その可視化
- 事実と仮説を分離した提示の骨組み

### やらないこと

索引の分散化 / Web UI / リアルタイム更新 / MassBank への書き戻し /
`msdial-interactive-app` への統合（別 MCP として並立させる）/
MS-DIAL の実行そのもの（`msdial-interactive-app` の領分）

## 3. 前提となる調査結果

設計判断の根拠。推測ではなく実コード・実データから確認した事実のみを置く。

### 3.1 上流の現行契約は MS 証拠を捨てている

`Interpret_mass_rdf/datas/kg_evidence_input.json` の構造:

```
kg_evidence[] : { inchikey, entities: { compounds, diseases, pathways,
                                        biospecimens, organisms, activities } }
```

上流は `InChIKey → メタデータ束` という関数であり、**照合スコアも一致フラグメントも
入口に存在しない**。`Interpret_mass_rdf/CLAUDE.md` の「MS証拠欠落」「同定は正しい前提」は
この契約の帰結である。

本設計はこの境界を広げる。

```
旧: InChIKey                                      → メタデータ束
新: {InChIKey, 一致フラグメント集合, スコア, 未説明ピーク} → メタデータ束 + バイアス説明
```

### 3.2 DB 収載バイアスが偽の説得力を生む実例

同 JSON の 2 件は、どちらも分子式 `C8H10N4O2` の別 InChIKey である。

| InChIKey | 出所 | diseases | pathways | activities |
|---|---|---|---|---|
| `RYYVLZVUVIJVGH-UHFFFAOYSA-N` | KNApSAcK / HMDB | 5 | 23 | 100 |
| `LPHGQDQBBGAPDZ-UHFFFAOYSA-N` | PubChem | 0 | 0 | 0 |

メタデータの厚さだけで見れば前者が圧勝するが、それは生物学的事実ではなく**DB 収載の差**である。
両者を分ける証拠はプロダクトイオンのパターンにしかない。
これが「並び順を決めるのはスペクトル証拠だけ」という §5.3 の契約の根拠である。

### 3.3 MassBank-RDF の TTL スキーマには 2 世代あり、食い違っている

`systemsomicslab/MassBank-RDF` を実査した結果:

| | 世代A: `draft_ttl/`（実体 54,183 件） | 世代B: `SPARQL/` の例が前提とするもの |
|---|---|---|
| ピーク値 | `mbo:has_peak [ mbo:mz 8.60349e+01 ; mbo:intensity ... ]` リテラル直付け | `?n a mbo:mz ; sio:SIO_000300 ?v` — SIO 具象化ノード |
| 中性損失 | 無い（precursor から計算する必要がある） | `mbo:neutral_loss` として格納済み |
| 経路長 | 1〜2 ホップ | Experiment→Spectrum→Peaks→Peak→Annotation→node の 6 ホップ |

**リーダは世代対応にする。** 世代の判定はグラフ内の述語の存在で行い、
どちらでも同一の内部表現に落とす。

### 3.4 世代A TTL が持っているフィールド（実サンプルより）

```
mbo:analytical_methods_and_conditions
    mbo:ac_instrument / instrument_type / ion_mode / ms_type / resolution
    fragmentation_mode / collision_energy / retention_time / column_name / solvent
mbo:compound
    schema:inChIKey / molecularFormula / smiles / inChI / name
    mbo:ch_exact_mass / ch_compound_class
    rdfs:seeAlso  cas: / chemspider: / pubchem:CID
mbo:has_peak              mz / intensity / relative_intensity / peak_sequential_number
mbo:has_peak_annotations  mz / tentative_formula / error
mbo:mass_spectram_meta_data
    precursor_mz_value / precursor_type_value / base_peak
mbo:pk_num_peak / mbo:pk_splash
dcterms:license           （CC0 / CC BY がディレクトリで分離されている）
```

### 3.5 コントリビュータ 10 件の実測

```
contributor          instr_type      resolution  ms_type  precursor  peaks  formula付き
Fac_Eng_Univ_Tokyo   EI-B            —           MS       —          22     0
Eawag                LC-ESI-ITFT     15000       MS2      [M+H]+      2     2
Keio_Univ            LC-ESI-QQ       —           MS2      [M-H]-     20     0
Athens_Univ          LC-ESI-QTOF     35000       MS2      [M+H]+    160   160
Chubu_Univ           LC-ESI-QIT      —           MS2      [M-H]-     42     0
ENTACT_AGILENT       ESI-QTOF        —           MS2      [M+H]+     27     0
LCSB                 LC-ESI-QFT      17500       MS2      [M-H]-      8     8
Antwerp_Univ         LC-ESI-QTOF     7000        MS2      [M+H]+      9     9
BS                   LC-ESI-QTOF     —           MS2      [M-H]-    250     0
JEOL_Ltd             MALDI-TOFTOF    —           MS2      —          49     0
```

レコード分布（上位）: Fac_Eng_Univ_Tokyo 12,378 / Eawag 10,604 / Keio_Univ 4,713 /
Athens_Univ 4,681 / HBM4EU 3,204 / Chubu_Univ 2,544。全 54,183。

ここから導かれる 3 つの帰結を §6 に反映する。

1. `mbo:resolution` は約半数で欠損する。層分けの主根拠にできない。
2. 最大コントリビュータ（全体の 23%）は **GC-EI の MS1 で precursor を持たない**。
   測定原理の異なるスペクトルが同一 TTL 群に混在している。
3. `tentative_formula` の有無は**サブミッター単位で全か無か**に分かれる。

### 3.6 SPARQL では照合が解けないことの実証

`SPARQL/get_inchikey_from_peakmz.sparql` は `VALUES ?targetmz {167.1121}` の
**完全一致**でピーク m/z から InChIKey を引いている。実測 m/z が小数点以下 4 桁まで
一致することはないので、これは実用にならない。トレランス付きの範囲検索と
強度を考慮したスコアリングは SPARQL の不得手な領域である。

したがって**スペクトル照合はローカルの数値エンジンで行い、SPARQL はメタデータ結合に使う**。

## 4. 方式

採用: **ハイブリッド（事前計算されたスペクトル索引 + ライブなメタデータ取得）**

| | A: 全部事前計算 | B: ライブ SPARQL 委譲 | **C: ハイブリッド（採用）** |
|---|---|---|---|
| スペクトル索引 | ローカル事前計算 | 都度 SPARQL | ローカル事前計算 |
| メタデータ | 索引に焼き込み | 都度問い合わせ | 都度問い合わせ |
| 背景分布の算出 | 可 | **事実上不可** | 可 |
| 逆引き速度 | ms | 秒〜分 | ms |
| KG の鮮度 | ビルド時点で凍結 | 最新 | 最新 |
| 再配布・ライセンス | 焼き込みで面倒 | 問題なし | 問題なし |
| 既存 SPARQL 資産 | 使い捨て | 活きる | 活きる |

C を採る理由は、境界が「重くて滅多に変わらないもの」と「軽くて変わる・権利が絡むもの」に
一致するからである。B では**全 MassBank における出現頻度**が算出できず、
本設計の中心的主張（濃縮度としての定量化）が成立しない。

副次的な効用として、**上流実装の所在が未確認のまま Phase 1 を完了できる**（§9）。

## 5. アーキテクチャ

### 5.1 モジュール境界

```
massbank-context/
  core/       FastMCP・設定・セッション状態・パス解決・共通ヘルパ   （leaf）
  spectral/   実験スペクトル取込 → 索引照合 → スコアリング         【事前計算索引】
  kg/         InChIKey → メタデータ束                             【ライブ・差し替え可能】
  bias/       目的 → カテゴリ → 想起拡張 / 濃縮度統計
  catalog/    msdial_repository_catalog の SQLite 読み取り        【読み取り専用】
  report/     事実と仮説を分離した提示の骨組み
  tools/      MCP 公開層
  index/      索引ビルド CLI（MCP からは呼べない。§7.1）
```

依存は上から下への一方向。`bias/` は `spectral/` と `kg/` の両方に依存するが、
**`spectral/` と `kg/` は互いを知らない**。`core/` は leaf であり、
上位モジュールを import してはならない（Lipidmix `mcp_core` と同じ規約）。

### 5.2 kg/ は 3 アダプタで切る

| 実装 | 使う場面 |
|---|---|
| `JsonFixtureAdapter` | 開発・テスト。既存 `kg_evidence_input.json` をそのまま読む |
| `SparqlAdapter` | QLever エンドポイントに問い合わせる |
| `UpstreamAdapter` | 研究室の既存実装をラップする |

**3 つとも返り値を現行契約に揃える**:
`{compounds, diseases, pathways, biospecimens, organisms, activities}`。
これにより `Interpret_mass_rdf` の下流（`organize_kg_data.py` / `interpret_kg.py`）が
そのまま接続でき、かつ上流の所在が未確定のまま設計を完成させられる。

### 5.3 証拠バンドの契約

候補 1 件は必ず次の 4 ブロックを持つ。**スコアを合成してはならない。**

```jsonc
{
  "inchikey": "...",
  "spectral_evidence": {
    "matched_fragments": 7,
    "total_query_peaks": 12,
    "weighted_score": 0.81,
    "cosine": 0.74,
    "unexplained_peak_ratio": 0.25,
    "band": "PASS"
  },
  "compound_context": { },
  "sample_context":   { },
  "bias": {
    "applied": true,
    "category": {"vocabulary": "HMDB:disease", "id": "...", "label": "..."},
    "rank_delta": 14,
    "entered_by": "recall_expansion"
  }
}
```

この型が保証する 3 つの性質:

1. **並び順を決めるのは `spectral_evidence` だけである。**
   `compound_context` が厚いだけで順位は上がらない。§3.2 の実例で言えば、
   メタデータ 0 件の `LPHGQDQBBGAPDZ` もフラグメントが合えば上位に出る。
2. **`weighted_score` は素の一致数ではない。**
   出現頻度の逆数で重み付けする（§6.5）。素の `matched_fragments` も併記し、
   重み付けが何をしたか見えるようにする。
3. **`rank_delta` と `entered_by` は必須である。**
   バイアスが押し上げた順位差が常に数値で見える。これが無ければ査読に耐えない。

`band` の定義:

| 値 | 条件 |
|---|---|
| `PASS` | プロダクトイオンが 2 本以上一致 |
| `WEAK` | プロダクトイオンが 1 本一致 |
| `PRECURSOR_ONLY` | precursor m/z のみ一致。プロダクトイオンの一致なし |

`PRECURSOR_ONLY` と `PASS` を同一視してはならない。
（Lipidmix の `analytical_checks.msms.band` と同じ思想）

## 6. データ設計

### 6.1 索引は acquisition_class で分割する

§3.5 の帰結 2 への対応。**測定原理が異なるスペクトルを同じ背景分布で数えることは無意味である。**
GC-EI のフラグメンテーションと LC-ESI-MS2 の CID プロダクトイオンは別物であり、
EI には precursor という概念自体がない。

`acquisition_class` を `instrument_type` と `ms_type` から導出し、索引を分割する。

| acquisition_class | 導出条件 | 中性損失軸 |
|---|---|---|
| `LC-ESI-MS2` | `instrument_type` が `LC-ESI-*` かつ `ms_type=MS2` | 使える |
| `GC-EI-MS` | `instrument_type` に `EI` を含み `ms_type=MS` | **使えない**（precursor なし） |
| `MALDI-MS2` | `instrument_type` に `MALDI` を含む | precursor 依存 |
| `other` | 上記以外 | 個別判定 |

**背景頻度も濃縮度統計も、同じ `acquisition_class` の中でだけ計算する。**
層を跨いだクエリは `incompatible_index` で拒否する（§8）。

### 6.2 resolution_class は instrument_type から導出する

§3.5 の帰結 1 への対応。`mbo:resolution` は約半数で欠損するため主根拠にできない。

| resolution_class | `instrument_type` の語尾 |
|---|---|
| `high` | `*FT`（ITFT, QFT — Orbitrap 系）、`*TOF`（QTOF, TOFTOF） |
| `low` | `QQ`（三連四重極）、`QIT` / `IT`（イオントラップ） |
| `unknown` | 判定不能 |

`mbo:resolution` が存在する場合は補強に使う。

**既定では `high` のみで索引を構築する。** 低分解能レコードを混ぜると、
それに合わせて全ビンが太り、高分解能側の情報が潰れるため。
`low` は別索引として後から追加できる形にする。

### 6.3 フラグメントは 3 軸で持つ

絶対 m/z だけで索引すると同族体で全てずれる。一方 `precursor − fragment`（中性損失）は
官能基とクラスに直結する（−18 H2O、−141 PE 頭部、−180 hexose など）。
**カテゴリ特異性は絶対 m/z より中性損失に出やすい**ため、中性損失軸を一級市民として扱う。

| 軸 | キー | 捉えるもの | 適用範囲 |
|---|---|---|---|
| `product_ion` | m/z ビン | 部分構造 | 全レコード |
| `neutral_loss` | (precursor − fragment) ビン | 官能基・クラス・頭部 | precursor を持つレコード |
| `formula` | `C9H10NO+` 等 | 装置非依存・解釈可能 | `tentative_formula` を持つレコードのみ |

**`formula` は一次キーにしない。** §3.5 の帰結 3 の通り、有無がサブミッター単位で
全か無かに分かれ、索引に穴が空くため。一次キーは質量ビン、`formula` は
**存在するときだけ付ける解釈用の注釈**である。自前で帰属計算はしない。

濃縮度統計は 3 軸それぞれ独立に回し、結果も軸ごとに分けて出す。

世代B の TTL では `mbo:neutral_loss` が格納済みのため読むだけでよい。
世代A では `precursor_mz_value − mz` で算出する。

### 6.4 ビン幅は実測から決める

`0.01 Da` のような決め打ちには根拠がない。索引ビルド時に、
**同一 InChIKey の複数レコード間で同じフラグメントの m/z がどう散らばるか**を測り、
その分布からビン幅を決める。

- 測定結果を `mass_error_profile` テーブルに `instrument_class` 別で保存する
- 採用したビン幅を `index_meta.bin_width` に記録する
- `index_status()` で読めるようにする

### 6.5 IDF の分母は unique InChIKey 数

MassBank には人気化合物のレコードが装置・CE 違いで多数入っている。
文書頻度をレコード数で数えると、**DB 収載バイアスがそのまま重みに乗る**
（§3.2 で見た罠が重み側に現れる）。

```
idf(f) = log( N_inchikey / df_inchikey(f) )      ※ df は unique InChIKey で数える
```

`mbo:pk_splash` は完全重複レコードの検出に使える（補助）。

### 6.6 濃縮度の背景を引数にする

カテゴリ C とフラグメント f の 2×2 分割表に対する Fisher 正確検定、
オッズ比、Benjamini-Hochberg FDR。標準的である。問題は**背景集団の選び方**にある。

| `background` | 出てくるもの |
|---|---|
| `all` | 「肝疾患カテゴリには脂質フラグメントが多い」= カテゴリの化合物クラス組成の再発見 |
| **`same_class`（既定）** | 「**同じ化合物クラスの中で**このカテゴリに特徴的なフラグメント」 |
| `inchikey[]` | ユーザーが対照集団を明示指定する |

`background="all"` を選んだ場合、返り値に
「クラス組成の再発見である可能性が高い」旨の注記を自動で付ける。

層別に使う化合物クラスは `compound_class(inchikey, source, label)` に持つ。
`mbo:ch_compound_class` は使えるレコード（例: Chubu_Univ の
"Natural Product; Lipid; Fatty acid..."）と使えないレコード（"N/A; Environmental Standard"）が
混在するため、**source 付きで複数持てる形**にし、ClassyFire / LIPID MAPS を
InChIKey 経由で別途取得できるようにする。

### 6.7 スキーマ（SQLite 単一ファイル）

```sql
records(record_id PK, accession, inchikey, precursor_mz, precursor_type,
        polarity, instrument_type, ms_type, acquisition_class, resolution_class,
        resolution_value, collision_energy, splash, license)
peaks(record_id, axis, bin, relative_intensity, formula NULL)
        -- axis: 'product_ion' | 'neutral_loss'
axis_background(acquisition_class, axis, bin, df_inchikey, idf)
compound_class(inchikey, source, label)
compounds(inchikey, molecular_formula, exact_mass, smiles, name)
index_meta(built_at, ttl_source, ttl_generation, n_records, n_inchikey,
           bin_width, resolution_filter, acquisition_classes)
mass_error_profile(instrument_class, percentile, delta_mz)
```

SQLite を選ぶ理由は、`msdial_repository_catalog` が同じ思想
（ローカル・単一ファイル・provenance 付き）で SQLite を採っており、
`catalog/` を読むコードと流儀を揃えられるため。
54,183 レコード × 平均数十ピーク ≒ 数百万行は SQLite の得意領域である。

### 6.8 スコアは IDF 重み付きカバレッジを主とする

```
weighted_score = Σ idf(f) for f in matched / Σ idf(f) for f in query_peaks
```

cosine / entropy similarity も併記するが、主役はカバレッジとする。理由:

1. `unexplained_peak_ratio` と直接対になり、「クエリのピークのうち何割が説明され、
   残りは何か」が一続きの説明になる。cosine の 0.81 は説明できない。
2. cosine は強度比に敏感で、装置間・CE 間の差を拾う。カバレッジは出現の有無なので頑健。

照合手順:

1. `acquisition_class` と polarity が一致する索引を選ぶ（不一致は `incompatible_index`）
2. precursor m/z ± tolerance で候補レコードを絞る
3. 各候補について query peaks と record peaks をビンで結合し、一致集合を得る
4. `weighted_score` を算出、`band` を判定（§5.3）

## 7. MCP ツール面

### 7.1 索引ビルドはツールにしない

TTL 全走査は分〜時間を要する。MCP ツールにすると LLM が会話の途中で起動して長時間固まる。
CLI に出す。

```
python -m massbank_context.index build --ttl-dir <dir> [--acquisition-class LC-ESI-MS2]
```

MCP 側は `index_status()` で状態を読むだけとする。

### 7.2 ツール一覧（10 本）

| # | ツール | 返すもの | 返さないもの |
|---|---|---|---|
| 1 | `index_status()` | 索引の有無・ビルド時刻・レコード数・acquisition_class / polarity のカバー範囲・bin_width | — |
| 2 | `load_spectra(path)` | 件数・polarity 分布・precursor m/z レンジ・id 一覧 | ピーク列（セッション保持） |
| 3 | `inspect_spectrum(id)` | 1 件のピーク列（丸め済み） | — |
| 4 | `set_objective(text)` | 実在語彙にマッチしたカテゴリ候補 | — |
| 5 | `match_spectra(ids?, tol, top_k, apply_bias)` | 候補テーブル（TSV・§5.3 の 4 ブロック） | 一致ピークの座標列 |
| 6 | `compound_context(inchikeys[])` | KG 束（畳み込み済み） | 生の直積 |
| 7 | `sample_context(inchikeys[])` | catalog 由来の検出実績（accession 付き） | — |
| 8 | `fragment_enrichment(category, axis, background)` | 特徴フラグメント上位（TSV・FDR 付き） | 全ビン |
| 9 | `plot_fragment_profile(categories[], axis)` | **PNG 画像** | 座標配列 |
| 10 | `draft_findings()` | 事実と仮説スロットの骨組み | — |

### 7.3 set_objective が「LLM 補助」の実体である

自然文の研究目的を受け取り、**LLM に自由なカテゴリを作らせない**。
`kg/` 側の実在語彙（KNApSAcK category、HMDB disease、pathway ラベル）を検索して候補を返し、
LLM は選ぶだけとする。

```jsonc
{
  "matched": [
    {"vocabulary": "HMDB:disease", "label": "Nonalcoholic fatty liver disease",
     "n_compounds": 312},
    {"vocabulary": "KNApSAcK:category", "label": "Hepatoprotective",
     "n_compounds": 47},
    {"vocabulary": "pathway", "label": "Primary bile acid biosynthesis",
     "n_compounds": 19,
     "warning": "n_compounds < 30 — 濃縮度統計の検出力が低い"}
  ],
  "unmatched_terms": ["胆汁酸"],
  "note": "カテゴリを選んで update_objective で確定してください"
}
```

これにより、目的バイアスが常に (a) 実在する語彙 ID に紐づき、
(b) `n_compounds` という検出力の目安を持ち、(c) 後から再現できる。
自由記述のカテゴリを許すとこの 3 つが全て失われる。

**`unmatched_terms` を分離する理由**: 「胆汁酸」のような化合物クラスは
カテゴリ（バイアス軸）ではなく**絞り込み条件**である。混ぜると `rank_delta` の意味が壊れる。

### 7.4 バイアスの適用方法は想起拡張のみ

バイアスは**候補プールをどこまで広げるかにだけ**効かせる。

- 通常は `weighted_score` の閾値でプールを作る
- バイアス有効時は、カテゴリ適合化合物を閾値下からも拾う（想起拡張）
- **並び順は `weighted_score` のまま変えない**
- 拾われた候補には `bias.entered_by = "recall_expansion"` を付ける
- 全候補に `bias.rank_delta`（バイアス無しなら何位だったか）を付ける

スコアの合成は行わない。ハードフィルタも行わない。

## 8. エラー処理

前提状態が無いときは例外ではなく機械可読な封筒を返す（Lipidmix の `missing_state` を継承）。

```jsonc
// (a) 索引が無い — ツールでは直せない。CLI を案内する
{"error": {"code": "missing_index",
           "required_action": "CLI: python -m massbank_context.index build --ttl-dir <dir>"}}

// (b) スペクトル未ロード — ツールで直せる
{"error": {"code": "missing_state", "state": "spectra",
           "required_tools": ["load_spectra"]}}

// (c) 索引の層が要求と合わない
{"error": {"code": "incompatible_index",
           "query_acquisition_class": "LC-ESI-MS2",
           "available": ["GC-EI-MS"],
           "message": "..."}}

// (d) メタデータ層が応答しない
{"error": {"code": "upstream_unavailable", "adapter": "SparqlAdapter",
           "message": "照合（match_spectra）は索引のみで動作するため継続可能"}}
```

**objective 未設定はエラーにしない。** `match_spectra` は `bias.applied: false` のまま動作し、
警告を 1 行添えるだけとする。理由は、バイアスが任意の追加機能であり、
主機能（スペクトル照合）はそれなしで完全に成立するためである。
ここを `missing_state` にすると「バイアスをかけないと結果が得られない」という
誤ったメッセージを構造で送ることになり、§5.3 の性質 1 と矛盾する。

`(d)` が示す通り、**メタデータ層が落ちても主機能は動き続ける**。方式 C を選んだ効用である。

## 9. 実装フェーズ

各フェーズは単体で使える状態で終わること。

| Phase | 内容 | 終了時点の到達点 | 上流への依存 |
|---|---|---|---|
| 1 | 索引ビルド CLI・世代A/B 対応リーダ・`index_status`・`load_spectra`・`inspect_spectrum`・`match_spectra`（バイアス無し） | MassBank 照合ツールとして単体で有用 | **なし** |
| 2 | `kg/` の `JsonFixtureAdapter`・`compound_context` | 候補にメタデータが付く | 既存 JSON のみ |
| 3 | `set_objective`・`bias/`・`rank_delta`・想起拡張 | 目的バイアスが乗る | 語彙の取得元が必要 |
| 4 | `fragment_enrichment`・`plot_fragment_profile`・`sample_context`・`draft_findings` | 順方向の可視化と仮説提示 | catalog |

Phase 1 が独立して価値を持つことが重要である。
**§11 の確認事項の返答を待たずに Phase 1 に着手できる。**

## 10. テスト戦略

- fixture はテスト自身が作る（Lipidmix の規約を継承）
- **例外**: 実 TTL を数レコード縮小コピーして同梱する（`tests/fixtures/ttl/gen_a/`, `gen_b/`）。
  TTL パースは外部仕様への追従であり、合成 fixture だけではスキーマ変更を検知できないため。
  世代A・世代B の両方を持つ
- 腐敗防止テスト:
  - `test_ttl_schema_contract.py` — 同梱 TTL から期待述語が全て読めるか。世代判定が正しいか
  - `test_statistics.py` — Fisher・IDF・BH-FDR を**手計算した既知解と照合**する。
    数値が静かに狂う場所であり必須
  - `test_acquisition_class.py` — §6.1 の導出表が実 `instrument_type` 文字列で正しく動くか
- ゴールデンテスト: 既知化合物のスペクトルを投げ、正解が top-k に入ることを検証
- `test_bias_invariants.py` — `apply_bias=True/False` で **`weighted_score` が変化しないこと**、
  `rank_delta` が実際の順位差と一致することを検証。§5.3 の性質 1 を機械的に守る

## 11. 未確定事項と、確認すべきこと

実コード調査により、当初の確認事項 4 件のうち 3 件が解決した。
残る 2 件（うち 1 件は新規）を研究室の MassBank-RDF 開発者に確認する。

| # | 確認事項 | 解決するもの | 影響するフェーズ |
|---|---|---|---|
| 1 | 配布する TTL は世代 A / B のどちらか。`draft_ttl/` は廃止されるのか | リーダの主対象。世代対応で実装するため設計は揺れない | Phase 1（運用のみ） |
| 2 | **HMDB / KNApSAcK をどのグラフにどうロードしているか** | `kg_evidence` の diseases・activities の出所。`SPARQL/` の例は PubChem BioAssay しか示していない | Phase 2–3 |

解決済み（確認不要）:

- ピークの述語名 → §3.3 / §3.4 で判明
- `instrument_type` / `resolution` の有無 → §3.5 で実測。ただし `resolution` は
  欠損が多く、§6.2 の導出方式に変更した
- 一致スコアを上流の返り値に足せるか → **問い自体が消滅**。索引を自前で作るため
  上流に追加を依頼する必要がない

その他の未確定事項（自分たちで決められるもの）:

- 層別背景に使う化合物クラスの出所（ClassyFire / LIPID MAPS / ChEBI のいずれか）。
  §6.6 の `compound_class(inchikey, source, label)` は複数 source を許す設計にしてあるため、
  後から追加・差し替えできる
- TTL の取得方法（`git clone` で 54,183 ファイルを一括取得するか、QLever 経由で吸い出すか）。
  Phase 1 の実装時に実測して決める

## 12. 規約（Lipidmix から継承）

戻り値の肥大は LLM の文脈を占め、結論を埋没させる。以下を最初から適用する。

- JSON は共通ヘルパで返す。`json.dumps(..., indent=2)` を書かない
- 行が並ぶ一覧は TSV（列名 1 回）
- float は丸めてから返す
- 全ツールに `structured_output=False` を付ける
  （付けないと content と structuredContent の両方で送られ 2 倍になる）
- 図はサーバで描いて PNG で返す。座標列を LLM に渡さない
  （Lipidmix 実測: 点列 183,578 字 ≒ 数万トークン → PNG 327 画像トークン）
- ピーク列・疎行列などの巨大な中間データはセッションに保持し、payload から外す
- 可変状態の正準は `core/` に置き、モジュール修飾で参照する

## 13. Phase 1 実施後に判明した、Phase 2 以降の課題

Phase 1（実装リポジトリ `C:/Users/yuu18/massbank-context`、commits `<empty-tree>..70d5de0`、185 tests）の
実施と最終レビューで判明した事項。実施中に潰した欠陥は本文に反映済みで、ここに残すのは
**意図的に Phase 1 の範囲外に置いたもの**だけである。

### 13.1 クエリピークの強度閾値（spec レベルの欠落）

最終レビューの指摘: `unexplained_peak_ratio` はクエリの全ピークを分母に数えるため、
デコンボリューション由来の低強度ノイズも含まれる。実際の MS-DIAL `.dcl` スペクトルでは、
**正しい同定であっても 0.8〜0.9 になる** — ノイズはどのライブラリレコードにも含まれないため。
§6.8 はこの指標を「カバレッジの物語の人間可読な半分」として売っているが、
実装のままでは良い一致に対して「ほとんど説明できていない」と読める。

数値自体は正しく、これは解釈の問題である。しかし対処には
`min_relative_intensity`（基準ピーク 999 スケールに対する相対閾値）または上位 N 本制限という
**新しい前処理段とパラメータ意味論**が要り、それは §6.8 に存在しない。
Phase 1 の範囲を実装側の裁量で広げないため保留した。

Phase 2 で扱うときの要件: 閾値は `load_spectra` と `match_spectrum` の間に置き、
**適用した値を payload に反映する**（そうしないと入力源をまたいで数値が比較できず、
後段の濃縮度統計がそのノイズを継承する）。

### 13.2 索引スキーマの版管理

`create_schema` は `CREATE TABLE IF NOT EXISTS` を使うため、既存の索引ファイルに対して
再ビルドすると旧スキーマが温存され、`sqlite3.OperationalError` の生のトレースバックで落ちる
（Phase 1 で `axis_background` に `polarity` 列を足した際に実際に発生する）。
CLI 限定で、復旧は索引ファイルの削除だが、**park した項目の中で実運用で最も踏まれやすい**。
`index_meta` にスキーマ版を持たせ、不一致なら「作り直せ」と明示的に言うのが素直。

### 13.3 質量誤差の同一レコード汚染は狭まったが消えていない

`measure_mass_error` はクラスタに 2 つ以上の**異なるレコード**を要求するようになったが、
spread は依然クラスタ全体の `max - min` で取る。レコード A が `[100.000, 100.030]`、
レコード B が `[100.004]` の場合、「異なる 2 レコード」の条件を満たしつつ 0.030 の spread が出る。
`choose_bin_width` は全クラスの最大を採るため、これが索引全体のビンを広げうる。

### 13.4 `neutral_loss` 軸は索引に入っているが誰も読んでいない

`axis_entries` は両軸を書き、`axis_background` は両軸の背景を持つが、
Phase 1 の照合は `product_ion` しか見ない。Phase 2 以降で
「照合に入れるのか、濃縮度統計専用に留めるのか」を明示的に決めること。

### 13.5 世代B の対応は検出できるようになったが、検証はされていない

§11 の確認事項 1 は未解決のまま。ただし Phase 1 で**検出**は入った:
ピーク 0 本のレコードを弾いて `no_peaks` に数え、TTL 世代ごとの解析件数を報告する。
実際の世代B グラフの入れ子が合成 fixture と違えば
「gen_b: N files, 0 peaks」として構築時に見える。開示だけでなく検出があるので、
実 TTL が手に入るまで安全に運べる。

### 13.6 `--acquisition-class` は未実装（plan の欠落）

§7.1 は挙げているが Phase 1 の plan が落とした。`derive_resolution_class` が
磁場型（`EI-B`）を `low` ではなく `unknown` に落とすため、
**GC-EI だけの索引を作る手段が現状ない**（`--resolution-filter low` では拾えず、
フィルタを無効化すると ion trap が LC-ESI の背景に混ざる）。
既定ビルドがコーパス最大のコントリビュータ（Fac_Eng_Univ_Tokyo、12,378 件 = 23%）を
設計上除外することは README に明記済み。
