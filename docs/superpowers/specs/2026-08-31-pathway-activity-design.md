# パスウェイ活性同定 設計書

**作成**: 2026-08-31
**位置づけ**: [massbank-context MCP サーバ 設計書](2026-08-28-massbank-context-mcp-design.md)（以下「本体 spec」）の
Phase 2〜4 を、目的「群間比較 → 生化学的パスウェイの変化」に向けて具体化したもの。
本体 spec の §5.2（`kg/` の 3 アダプタ）・§8（エラー契約）・§9（フェーズ）・§12（戻り値の規約）を継承する。
節番号が `spec §x` と書かれているときは本体 spec を指す。

**変更対象リポジトリは 2 つある**: `massbank-context` と `Lipidmix_with_LLM`（MCP サーバ名 `ms-data-parser`）。

---

## 1. 目的

2 群の実験データから、**どの生化学的パスウェイに変化が集まっているか**を、
根拠の数値を伴って提示する。

**「活性化した」「不活性化した」と計算で断定しない。** 計算が返すのは
「そのパスウェイの構成メタボライトが、どちらの群で高いか」という記述であり、
活性の解釈は仮説スロットに隔離する（§10）。理由は §3.1 にある。

## 2. スコープ

### やること

- ms-data-parser の群間比較結果を InChIKey 付きで書き出す export ツール
- massbank-context 側でそれを読み、PubChem のパスウェイを引き、濃縮と方向を計算する
- 事実と仮説を分離した提示の骨組みを返す

### やらないこと

- 反応の向き（基質・産物）を使った真の活性推論（§3.1・§13.5）
- 疾患・生理活性メタデータ（HMDB・KNApSAcK。ラボ内エンドポイント待ち・§13.4）
- MS/MS による同定の独立検証を統計に反映すること（今回は印も付けない・§3.5）
- 3 群以上の同時比較（ms-data-parser が 2 群のみ対応。因子トークンで 2 群を切り出す）

---

## 3. 前提となる調査結果

すべて 2026-08-31 に実測した。推測ではない。

### 3.1 PubChem のパスウェイは「参加」しか持たない

`https://rdfportal.org/pubchem/sparql` は認証なしで到達でき、
`FROM <http://rdf.ncbi.nlm.nih.gov/pubchem/pathway>` に対する

```sparql
?pathway a vocab:Pathway ; obo:RO_0000057 ?cid ; dcterms:title ?pathway_label .
```

が実際に返る（InChIKey → CID → Pathway）。

`obo:RO_0000057` は **has participant** であり、その化合物が基質か産物か、
反応の向き、量論を**一切持たない**。したがって代謝物の増減からパスウェイの活性は
直接は導けない。これが §1 の「断定しない」の根拠であり、§10 の `limits` に
毎回出す固定文の中身でもある。

### 3.2 パスウェイの被覆は 14.5%

massbank-context の索引が持つ 8,091 InChIKey から 200 件を無作為抽出（seed 0）し、
50 件ずつ 4 バッチで照会したところ、パスウェイが付いたのは **29 件 = 14.5%**
（標準誤差 約 2.5%）。中心代謝は当たり、脂質はほとんど当たらない。

**この数字は payload に必ず出す。** 出さないと「パスウェイが付かない ＝ 変化が無い」と
読まれる。付かなかった 85% は「調べていない」である。

### 3.3 organism を混ぜると 4 割が非対象生物になる

パスウェイ行には organism が必ず付く（実測 124/124）。ただし内訳は
human 67・**大腸菌 18・緑膿菌 10・シロイヌナズナ 6**・その他であり、
哺乳類試料の解析に非哺乳類のパスウェイが 4 割混ざる勘定になる。

またこの標本では **pathway URI 90 件 = pathway ラベル 90 件**で、
URI とラベルは 1:1 だった。ラベルの名寄せは不要で、URI を集計単位にできる。

### 3.4 mzTab から InChIKey は取れない

MS-DIAL 5.5 の mzTab-M 2.0.0-M 出力（`datasets/2_lipidome_lcms/NEG/`、SML 714 行）の実測:

| 項目 | 実測 |
|---|---|
| `inchi` | **0 / 714**（全行 `null`） |
| `smiles` | 412 / 714 |
| `chemical_name` | 497 / 714 |
| `study_variable` | 4 群（control / LPS / G / ILG）、assay 60 本 |
| `abundance_assay[n]` / `abundance_study_variable[n]` | あり |
| `opt_global_ms2_presence` | True 686 / 714 |
| `opt_global_spectrum_matched` | True 359 / 714 |
| `reliability`（MS-DIAL 独自コード） | 410:265 / 999:217 / 530:138 / 430:89 / 420:5 |

mzTab は**群構造をファイル自体に持つ**という強い利点がある一方、
パスウェイ照会の鍵である InChIKey を持たない。同じ MS-DIAL の `.arf2` は
InChIKey を持つ（Key16）のに、mzTab 出力では落ちている。

**したがって今回の実装は `.arf` / `.arf2` を入力とする。** mzTab への移行方針は §13.1 に記録する。

### 3.5 `.arf` と `.arf2` は MasterAlignmentID で結合できる

ms-data-parser の `lipidmix/arf/tools.py` の `_annotate_with_names` が既に
この結合を実装している。feature 名 `Spot_<id>_<prop>` から ID を取り、
`AlignmentResult_<timestamp>` の語幹が一致する**兄弟 `.arf2` だけ**を許す
（MasterAlignmentID はアラインメント実行ごとに振り直されるため、別バッチの
`.arf2` を引くと ID 対応が黙って崩れる）。

一方 `.arf2` は **MS/MS 取得フラグを持たない**（`arf2_annotate_identities` 自身が
「`has_msms=False`・バンドは `UNKNOWN`」と明記している）。`verify_peak_annotation` は
`.pai2` の**サンプル別ピーク**に対する検証であり、アライメント・スポット単位ではない。
スポット ↔ サンプルピークの対応付けは別実装を要するため、**今回は MS/MS の印を契約に載せない**。
`.mzTab` へ移れば `opt_global_ms2_presence` / `opt_global_spectrum_matched` で自然に埋まる（§3.4）。

### 3.6 ms-data-parser の log2FC は慣習と逆向き

`lipidmix/analysis/differential.py` の現行実装は次のとおりで、**正なら group_a が高い**。

```python
def _log2fc(mean_a, mean_b, pseudo_count):
    return math.log2((max(mean_a, 0) + pseudo_count) / (max(mean_b, 0) + pseudo_count))
# log_transform=True の分岐も fc = lm_a - lm_b
```

一般的な慣習（`log2FC = log2(比較対象 / 基準)`、正 = 上昇）と逆である。
**本設計では ms-data-parser 側を慣習に合わせて修正する**（§6.1 (b)・§12 Phase A）。

---

## 4. 方式（決定事項）

| # | 論点 | 決定 | 理由 |
|---|---|---|---|
| 1 | 3 つの層をどこに置くか | 群統計は ms-data-parser、同定の文脈とパスウェイと提示は **massbank-context** | spec §5.1 の `kg/` `report/` の想定にそのまま収まる |
| 2 | 活性の主張範囲 | **方向集約まで。活性の解釈は仮説スロットへ隔離** | §3.1。参加関係しか無いデータで活性は導けない |
| 3 | 入力ファイル | 今回は `.arf` / `.arf2`。将来 `.mzTab` + `.msp` | §3.4 |
| 4 | サーバ間の受け渡し | ms-data-parser に **export ツールを新設**し、TSV ファイルで渡す | 結合はデータの形を知っている側で行う（§3.5）。会話経由は payload が肥大する（spec §12） |
| 5 | MS/MS 検証の扱い | 今回は使わない（印も付けない）。層別と感度解析は将来 | §3.5。二重実装を避ける |
| 6 | 統計 | **濃縮（Fisher + BH）と方向（記述 + Wilcoxon）を別々に併記**。合成しない | spec §5.3「スコアを合成してはならない」 |
| 7 | パスウェイの持ち方 | **ライブ SPARQL + ローカルキャッシュ** | 解析対象の InChIKey は `.arf2` 由来で、MassBank 索引に無いものを大量に含む。索引に焼き込む方式では引けない |

---

## 5. アーキテクチャ

### 5.1 モジュール境界

massbank-context に 4 モジュールを足す。既存の一方向依存は崩さない。

```
core/        設定・セッション・モデル・直列化・エラー契約        （leaf・既存）
stats/       Fisher / BH-FDR / Wilcoxon 符号順位                （leaf・新規）
experiment/  差次的結果の取り込み（export 契約のリーダ）         （新規）
kg/          InChIKey → パスウェイ（SPARQL + SQLite キャッシュ） （新規・spec §5.2）
report/      事実と仮説を分けた骨組み                           （新規）
tools/       MCP 公開層                                        （既存）
```

依存は `tools/ → report/ → {experiment/, kg/, stats/} → core/`。

- `stats/` は純関数だけを置き、`core/` にも依存しない。§11 の既知解テストがここに集まるため、外部状態を持たせない
- `experiment/` と `kg/` は互いを知らない（spec §5.1 の `spectral/` と `kg/` の関係と同じ）
- 依存は `rdflib` と `mcp` の 2 つから増やさない。Fisher・BH・Wilcoxon は自前で書く

### 5.2 データフロー

```
[ms-data-parser]  .arf（定量・群） + .arf2（InChIKey・Ontology）
      │  arf_export_differential   ← .arf ↔ .arf2 の結合はここ
      ▼
  differential.tsv                 ← 両サーバが知る唯一の契約
      │  ファイルパスを渡す
      ▼
[massbank-context] load_differential(path)     群名・件数・InChIKey 付与率
      ▼            resolve_pathways()          照会 → キャッシュ。被覆と organism 内訳
      ▼            pathway_activity(organism)  濃縮と方向。TSV
      ▼            draft_pathway_findings()    事実と仮説のスロット
```

---

## 6. データ設計

### 6.1 export 契約（`differential.tsv`）

**新ツール**: ms-data-parser の `arf_export_differential(group_a, group_b, output_path)`。
`arf_preprocess` → `arf_differential` の後に呼ぶ。

1 ファイルで完結させる。`#` 始まりのメタ行に来歴を置き、その後に TSV 本体を置く。

```
# contract_version = 1
# exported_at = 2026-08-31T12:00:00Z
# source_arf = <path>
# source_arf2 = <path>
# group_a = gf_AIN    n_a = 5
# group_b = gf_HFD    n_b = 5
# log2fc_sign = positive means group_b is higher
# q_threshold = 0.05  log2fc_threshold = 1.0  log_transform = true
# preprocess = {normalize: ..., blank_removed: ..., qc_excluded: ...}
# n_features_total = 3412   n_with_inchikey = 894   n_unannotated = 2518
spot_id	name	name_source	ontology	inchikey	inchikey_source	msi_level	mz	rt	log2fc	p_value	q_value	mean_a	mean_b	significant
```

契約の要点は 4 つある。

**(a) 有意な行だけでなく、InChIKey が付いた全行を出す。**
Fisher の背景は「検出され、かつパスウェイ照会を実施した化合物」でなければならない。
有意な行だけ渡されると背景が作れない。`significant` はフラグとして持ち、
判定に使った閾値はメタ行に残す。

**(b) 符号の向きを ms-data-parser 側で慣習に合わせる。**
`log2fc = mean_b − mean_a`（log 空間）に修正し、**正 = group_b が高い = 上昇**とする。
`group_a` が基準（対照）、`group_b` が比較対象、という役割を docstring で固定する。
massbank-context 側は反転せずそのまま読む。

これは既存 MCP サーバの破壊的変更であり、過去に出した解析結果の向きが反転する。
payload の caveats に向きを毎回明示し、`differential_contract_version` を持たせて
古い出力と区別できるようにする。修正対象:

- `lipidmix/analysis/differential.py` — `_log2fc` と `two_group_test` の 2 分岐
- `lipidmix/analysis/differential.py` — `summarize_two_group` の `n_up` / `n_down`
- `docs/output_format/arf.md` — 「正=群Aで高い」の記述 3 箇所
- `tests/test_differential.py` / `tests/test_differential_tools.py` / `tests/test_volcano_plot.py`

volcano の `up` / `down` ラベルは符号が直れば意味が正しくなる（現状は
`up` が「群 A で高い」を指しており、ラベルと実体が食い違っている）。

**(c) `n_unannotated` を必ず出す。**
注釈が付かなかった行は「変化が無かった」のではなく「調べていない」。
この数字が payload から消えると、§3.2 の被覆と同じ誤読が起きる。

**(d) 兄弟 `.arf2` が無ければ export は失敗させる。**
InChIKey 空欄の行を黙って出すと、下流で「パスウェイが無い化合物」と区別が付かなくなる。

`msi_level` は `.arf2` 由来の注釈確度（`arf2_annotate_identities` の MSI クラス上限）であり、
**MS/MS の有無ではない**。契約のコメント行にその旨を明記する（§3.5）。

### 6.2 パスウェイキャッシュ

`index/pathway_cache.sqlite`。索引本体とは別ファイルにする（索引の再ビルドで消えない）。

```sql
pathway_hit(inchikey, pathway_id, pathway_label, taxid)
lookup_log(inchikey, endpoint, fetched_at, status)   -- ok / empty / failed
organism(taxid, label)
cache_meta(key, value)                                -- contract_version, endpoint
```

**`lookup_log` が要点である。** 「照会したがヒット 0」と「まだ照会していない」を
区別できないと、被覆の分母が壊れる。`pathway_hit` に行が無いことだけでは両者が同じに見えるため、
照会した事実そのものを別テーブルに記録する。`status='failed'` は次回の再試行対象として残る。

**再取得は明示操作のみ**とする（`refresh=True`）。日数ベースの TTL は付けない。
PubChem のパスウェイは滅多に変わらない一方、日付でキャッシュが勝手に切れると
「先週と同じ解析が違う結果を出す」ことになり、再現性が壊れる。

照会は 50 件/バッチ。バッチが失敗したら**その 50 件を個別に再試行**し、
どの InChIKey で落ちたかを特定できるようにする（`MassbankRdf` の設計を踏襲）。

`organism` テーブルは taxid のラベル解決に使う。パスウェイ照会とは別のクエリで
一度だけ引き、キャッシュする。

### 6.3 集計単位と organism

集計単位は **pathway URI**（§3.3 の実測により、URI とラベルは 1:1）。

**organism は `pathway_activity` の必須引数とする。** 既定で全生物を混ぜるのも、
黙って human に絞るのも避ける。指定が無ければ `incompatible_organism` を返し、
選べる taxid と件数を示して拒否する。層を跨いだ照合を黙ってやらない、という
既存の思想（spec §6.1 の `acquisition_class`）と揃える。

---

## 7. 統計

### 7.1 濃縮（Fisher 片側 + BH-FDR）

背景を次のように取る。

> **背景 = 差次的結果に InChIKey が付き、照会が完了し、指定 organism のパスウェイを 1 つ以上持つ化合物**

パスウェイを 1 つも持たない化合物は背景に含めない。ORA の定石どおりで、
含めると全パスウェイが一律に有意寄りになる。除外件数は `n_no_pathway` として必ず出す。

2×2 は次のとおり。

|  | そのパスウェイのメンバー | 非メンバー |
|---|---|---|
| 差次的（`significant`） | k | n − k |
| 非差次的 | K − k | N − n − K + k |

N = 背景の化合物数、K = 背景のうちそのパスウェイのメンバー数、n = 背景のうち差次的な数。
片側（濃縮方向）で検定し、BH-FDR は**検定したパスウェイ全体**で補正する。

メンバー数が `min_members` 未満のパスウェイは**検定しない**。検出力が無いのに
BH の分母だけ増やして他のパスウェイの q 値を悪化させるためである。
ただし表からは消さず、`tested=false` の行として残す。

### 7.2 方向（記述 + Wilcoxon 符号順位）

方向は、そのパスウェイの**背景メンバー全員**（有意でないものも含む）の `log2fc` から取る。
有意なものだけ集めると、閾値の内側にある一貫した傾向が見えなくなる。
`median_log2fc` も同じ集合（背景メンバー全員）の中央値である。

「有意メンバー」とは契約の `significant == true` の行を指す。その判定は
ms-data-parser 側で行われ（`q_value ≤ q_threshold` かつ `|log2fc| ≥ log2fc_threshold`）、
閾値はメタ行に残る。massbank-context 側でこの判定をやり直さない。

`direction` は恣意的な閾値を避け、次の定義に固定する。**これは検定ではなく記述である。**

| 値 | 条件 |
|---|---|
| `UP_DOMINANT` | 有意メンバーが 2 件以上あり、すべて上昇（group_b が高い） |
| `DOWN_DOMINANT` | 有意メンバーが 2 件以上あり、すべて低下 |
| `MIXED` | 有意メンバーに上昇と低下が両方ある |
| `SINGLE_MEMBER` | 有意メンバーが 1 件だけ |
| `NO_SIGNIFICANT_MEMBER` | 有意メンバーが 0 件 |

Wilcoxon 符号順位検定（中央値 0 に対する）の p 値は別列に置き、`direction` とは合成しない。
メンバー数が閾値未満なら走らせず `wilcoxon_tested=false` とする。

### 7.3 並び順と多重検定

**並び順は濃縮の q 値だけで決める**（同値は `pathway_id` で決定的に）。
方向の p 値では並べ替えない。spec §5.3 の「並び順を決めるのは 1 つだけ」をここでも守る。
BH 補正は濃縮と Wilcoxon で**別々に**掛ける。

出力（TSV）:

```
pathway_id	pathway_label	taxid	n_members	n_significant	n_up	n_down	median_log2fc	direction	fisher_p	fisher_q	tested	wilcoxon_p	wilcoxon_tested
```

payload のヘッダには `n_background` / `n_no_pathway` / `n_pathways_total` /
`n_pathways_tested` / `organism` / `group_a` / `group_b` と、
**符号の向きの説明文**（「上昇 = group_b が高い」）を毎回載せる。

`pathway_coverage` の定義を固定する。**分母は照会が完了した InChIKey 数
（`lookup_log.status` が `ok` または `empty`）、分子はそのうち指定 organism の
パスウェイを 1 つ以上持つ InChIKey 数**である。`failed` は分母に含めない
（§9 のとおり `failed` が残っている間は統計を走らせない）。

`top_k` は表示行数の上限であり、**検定した件数ではない**。打ち切ったときは
`truncated: true` と全件数を必ず payload に出す（spec §12）。
`tested=false` の行も表示対象に含まれるため、`top_k` を小さくすると
検定しなかったパスウェイが見えなくなる — その件数は `n_pathways_total` と
`n_pathways_tested` の差として常に読める。

---

## 8. MCP ツール面

massbank-context に 4 本追加する（既存 4 本と合わせて 8 本）。
すべて `structured_output=False`、戻り値は `core/serialization.py` 経由（spec §12）。

| # | ツール | 返すもの | 返さないもの |
|---|---|---|---|
| 1 | `load_differential(path)` | 件数・群名と n・InChIKey 付与率・`n_unannotated`・符号の向き | 行そのもの（セッション保持） |
| 2 | `resolve_pathways(refresh=False)` | 被覆・キャッシュヒット数・organism 別内訳・取得日時・エンドポイント | パスウェイ行の全量 |
| 3 | `pathway_activity(organism, min_members=3, top_k=30)` | §7.3 の TSV | メンバー化合物の全列挙 |
| 4 | `draft_pathway_findings()` | §10 の骨組み | 解釈そのもの |

ms-data-parser に 1 本追加する。

| ツール | 返すもの |
|---|---|
| `arf_export_differential(group_a, group_b, output_path)` | 書き出したパス・行数・`n_with_inchikey`・`n_unannotated`・符号の向き |

---

## 9. エラー契約

すべて spec §8 の機械可読エンベロープ。**本文 JSON で返す**
（使用中の MCP SDK では `isError=true` の経路が structuredContent を捨てるため）。

| コード | 出る場面 |
|---|---|
| `missing_state` | `load_differential` 未実行 / 照会未完了の InChIKey が残ったまま `pathway_activity` を呼んだ |
| `contract_mismatch` | `contract_version` が未対応、必須列やメタ行が欠ける |
| `incompatible_organism` | organism 未指定、または照会結果に無い organism。選べる taxid と件数を返す |
| `upstream_unavailable` | rdfportal に到達できない（spec §5.2 で名前だけ決まっていたもの） |

**照会が完了していない InChIKey が 1 件でも残っていれば `pathway_activity` は走らせない。**
部分的な母集団で統計を回すと、被覆が見かけ上変わって濃縮の p 値が動く。
`resolve_pathways` をやり直せ、という `missing_state` を返す。

ms-data-parser 側の export は、兄弟 `.arf2` が見つからなければ失敗させる（§6.1 (d)）。

---

## 10. 提示層

`report/` は LLM を呼ばない。**空のスロットを持つ骨組みを返すだけ**で、埋めるのは LLM である。

```jsonc
{
  "facts": {
    "comparison": {"group_a": "gf_AIN", "group_b": "gf_HFD", "n_a": 5, "n_b": 5,
                   "log2fc_sign": "正 = gf_HFD が高い（上昇）"},
    "population": {"n_features_total": 3412, "n_with_inchikey": 894,
                   "n_unannotated": 2518, "n_background": 129,
                   "n_no_pathway": 765, "pathway_coverage": 0.144},
    "provenance": {"endpoint": "https://rdfportal.org/pubchem/sparql",
                   "fetched_at": ["2026-08-31T…", "2026-08-31T…"],
                   "organism": "TAXID10090"},
    "top_pathways": [ /* §7.3 の行。数値のみ、解釈なし */ ]
  },
  "hypotheses": [
    {"slot": "pathway_activity", "pathway_id": null, "claim": null,
     "evidence_rows": [], "required": ["pathway_id", "claim", "evidence_rows"]}
  ],
  "limits": [ /* 下記 5 項目・固定 */ ]
}
```

**「活性化」「不活性化」という語は `hypotheses` の中にしか置けない。**
`facts` 側の語彙は「上昇（group_b が高い）」「低下」で、群名を必ず伴わせる。
方向を語る文が群名なしで独り歩きしないようにするためである。

`limits` は毎回・固定で出す。省略可能にすると、いずれ省略される。

1. パスウェイ被覆は N%（実測値を動的に埋める）。パスウェイが付かなかった化合物は
   「変化が無かった」のではなく「調べていない」
2. **PubChem のパスウェイは `has participant` であり、基質か産物か・反応の向き・量論を持たない。
   代謝物の増減からパスウェイの活性は直接は導けない**
3. 化合物同定は MS-DIAL のライブラリ注釈由来で、MS/MS による独立検証はしていない
4. メンバー数が閾値未満のパスウェイは検定していない（`tested=false`）
5. 異性体は区別していない（MSI Level 2）

2 番目が本質である。§4 の決定 2 を、payload の構造そのもので担保する。

---

## 11. テスト戦略

- **`stats/`** — 教科書の**手計算した既知解**と照合する（Fisher の厳密値、BH の既知例、
  Wilcoxon の小 n 厳密値）。scipy 依存が無いため参照実装との突き合わせができず、
  既知解でしか守れない。spec §10 の `test_statistics.py` の思想を継承する
- **`experiment/`** — 契約リーダ。列欠落・`contract_version` 違い・符号メタ行の欠落を
  **必ず落とす**テストを含める
- **`kg/`** — 既定では実エンドポイントを叩かず、記録済み応答 JSON で回す。
  **ただし腐敗防止テストを 1 本だけ、実エンドポイントに当てる形で別マークで用意する**
  （既定の `pytest` では走らせない）。SPARQL のスキーマ変更は合成 fixture では検知できず、
  spec §10 が実 TTL を同梱したのと同じ理由による
- **`kg/cache.py`** — 「未照会」と「ヒット 0」の区別が壊れないこと、`refresh` の挙動
- **`report/`** — `limits` が必ず 5 項目出ること、`facts` 側に解釈語（「活性化」等）が
  混入しないことを**機械的に**検査する
- **統合** — 小さな合成 `differential.tsv` とキャッシュ済みパスウェイから、
  期待どおりの並び順が出ること
- **ms-data-parser 側** — 符号修正に伴い、`log2fc` の向きを既知解で固定するテストを追加する
  （群 B が明らかに高い合成データで `log2fc > 0` になること）

---

## 12. 実装フェーズ

各段階は単体で使える状態で終わること（spec §9 の方針）。

| Phase | 内容 | 終了時点の到達点 | 触るリポジトリ |
|---|---|---|---|
| A | log2FC の符号修正、`arf_export_differential`、`experiment/`、`load_differential` | InChIKey 付き差次的結果が両サーバ間で渡る | 両方 |
| B | `kg/`（アダプタ + キャッシュ）、`resolve_pathways` | 化合物にパスウェイが付き、被覆が数字で出る | massbank-context |
| C | `stats/`、`pathway_activity` | 濃縮と方向が出る | massbank-context |
| D | `report/`、`draft_pathway_findings` | 事実と仮説が分離された骨組みが出る | massbank-context |

Phase A は破壊的変更（§6.1 (b)）を含むため、単独で完結させ、
ms-data-parser のテストを通してから B に進む。

---

## 13. 将来課題

### 13.1 `.mzTab` + `.msp` への移行

利用者の方針は「将来は `.mzTab` を第一優先の入力、スペクトルは `.msp`」である。
移行時に変わるものを、判明している範囲で記録する。

- **群構造は mzTab 側に移る。** `study_variable[n]` と `assay_refs` が群定義を持ち、
  `abundance_assay[n]` が試料ごとの生強度を持つ（§3.4）。ms-data-parser の
  前処理・検定を経由せず、massbank-context だけで完結させる選択肢が生まれる
- **InChIKey が消える。** mzTab には `inchi` 列があるが実測は全行 `null` で、
  `smiles` が 58% にあるのみ（§3.4）。SMILES → InChIKey 変換（RDKit 依存の追加）か、
  `.msp` の MS/MS を MassBank と照合して InChIKey を得る経路が要る
- **MS/MS の印は逆に手に入る。** `opt_global_ms2_presence` と
  `opt_global_spectrum_matched` が使える（§3.4）
- **未確認**: MS-DIAL のアライメント `.msp` と mzTab の `SML_ID` の対応。
  `SME.spectra_ref` は `ms_run[n]:ms1scanID=… ms2scanID=…` 形式で、
  エクスポートされた `.msp` を直接指してはいない。移行の設計前に実物で確認すること

### 13.2 MS/MS による確度の層別と感度解析

今回は MS/MS 検証を使わない（§4 決定 5）。将来は化合物ごとに
`MSMS_PASS` / `MSMS_WEAK` / `PRECURSOR_ONLY` / `NOT_IN_INDEX` を付け、
主解析（注釈全体）と感度解析（`MSMS_PASS` のみ）を併走させる。
両方で同じパスウェイが出れば主張が強くなる。スコアは合成しない（spec §5.3）。

前提として、アライメント・スポット ↔ サンプル別ピークの対応付けが要る（§3.5）。
massbank-context の照合ヒット率は肝臓リピドミクスで 12.9% であり、
`MSMS_PASS` のみの母集団は小さくなる。感度解析としてのみ意味を持つ。

### 13.3 根拠必須の機械的強制

LLM が `evidence_rows` を空のまま `claim` を書くことをツール側で禁じるには、
書き戻しを受ける `submit_findings` が要る。今回は骨組みを返すだけに留める。

### 13.4 疾患・生理活性

HMDB と KNApSAcK の既定エンドポイントは `http://10.242.145.97:3030/…` の
プライベート IP であり、研究室ネットワーク外からは到達しない。
HMDB は公開 SPARQL 無し、KNApSAcK は公開 API 無し。到達性の確認待ち。
到達できれば `kg/` に同じ形のアダプタを足すだけで済む（spec §5.2 の契約に揃える）。

### 13.5 反応の向きを使った真の活性推論

KEGG / Rhea / Reactome から基質・産物の役割を取得できれば
「基質減・産物増 → 活性化」の推論が可能になる。DB を 1 層足すことになり、
被覆（§3.2）はさらに絞られる。§3.1 の制約を外せる唯一の道であり、
これが実現するまで §1 の「断定しない」は維持する。

---

## 14. 参照

- 本体設計書: [massbank-context MCP サーバ 設計書](2026-08-28-massbank-context-mcp-design.md)
- 到達度調査: [群間比較からパスウェイ変化まで、今どこまで届くか](../../surveys/2026-08-31-massbank-context-capability-assessment.md)
- 上流比較調査: [MassBank-RDF 上流リポジトリ比較調査](../../surveys/2026-08-31-massbank-rdf-upstream-comparison.md)
- Step 1 実施結果: [既知解テスト](../../verification/2026-08-31-selfcheck-step1-results.md)
- クエリの原型: `MassbankRdf/massbank_rdf/services/kg/query_builders/pubchem_query_builder.py`
  の `build_pubchem_pathway_query`（`MotohiroOGAWA/MassbankRdf` `feature/msp-kg` `a7a0233`）
