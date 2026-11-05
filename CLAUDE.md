# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## このリポジトリについて

実験 MS/MS スペクトルを MassBank-RDF と照合し、候補化合物を返す MCP サーバ（Python）。
コード・コメント・docstring・コミットメッセージはすべて**日本語**で書く。

正準の設計書は [docs/superpowers/specs/2026-08-28-massbank-context-mcp-design.md](docs/superpowers/specs/2026-08-28-massbank-context-mcp-design.md)。
コード中の `spec §6.4` のような参照はこの文書の節番号を指す。**設計上の判断を変える前に
必ず該当節を読むこと。** §13 に Phase 2 以降へ意図的に残した課題がある。

現状は Phase 1（索引ビルド + バイアス無し照合）に加え、パスウェイ活性同定の
Phase A〜D（`experiment/` `kg/` `stats/` `report/`、ツール 4 本）まで。
設計は [docs/superpowers/specs/2026-08-31-pathway-activity-design.md](docs/superpowers/specs/2026-08-31-pathway-activity-design.md)。
外来の分子（薬物・代謝物など）の走査（`xenobiotic/`、ツール 3 本）も入っている。
spec §8 の較正は r3 で改訂した期待に合格したが、手順書には入っていない（設計書 §10.0）。
設計は [docs/superpowers/specs/2026-11-04-xenobiotic-screen-design.md](docs/superpowers/specs/2026-11-04-xenobiotic-screen-design.md)。
`bias/` `catalog/` と spec §7.2 のツール 10 本のうち残りは未実装。

## コマンド

Python は `C:/Python314/python.exe`（`python` は PATH に無いことがある）。

```bash
C:/Python314/python.exe -m pytest tests -q
```

単体テスト: `C:/Python314/python.exe -m pytest tests/test_matcher.py::test_name -q`

索引ビルド（TTL 全走査は分〜時間。MCP ツールではなく CLI。spec §7.1）:

```bash
C:/Python314/python.exe -m massbank_context.index build --ttl-dir <TTLのディレクトリ> --output index/massbank.sqlite
```

`--limit N` で先頭 N ファイルだけの動作確認用索引を作れる（`index_meta.limit` に記録される）。

サーバ起動: `C:/Python314/python.exe -m massbank_context.server`。
索引パスは環境変数 `MASSBANK_CONTEXT_INDEX` で上書き可。

lint 設定は無い。依存は `rdflib` と `mcp` のみ。

## アーキテクチャ

依存は上から下への一方向。`core/` は leaf で、上位モジュールを import してはならない。

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
xenobiotic/ 外来の分子の登録と走査（質量関係・MS/MS・検出パターンの 3 つの筋から水準を付ける）
tools/      MCP 公開層（実体はここ、server.py は登録するだけ）
server.py   FastMCP ファサード
```

### 二層の分類が索引の分割単位（spec §6.1 / §6.2）

`core/classify.py` が `instrument_type` / `ms_type` から `acquisition_class`
（`LC-ESI-MS2` / `GC-EI-MS` / `MALDI-MS2` / `other`）と `resolution_class`
（`high` / `low` / `unknown`）を導出する。これは索引の分割単位であり、背景分布の
計算範囲でもある。**測定原理の異なる層を跨いだ照合はしない** — 要求された層が索引に
無ければ `incompatible_index` エラーで明示的に拒否し、黙って別層で代用しない。

既定ビルド（`resolution_class == "high"`）はコーパス最大のコントリビュータ
`Fac_Eng_Univ_Tokyo`（23%）を設計上除外する。GC-EI 磁場型が `unknown` に落ちるため。
取りこぼしではない（README と spec §6.1/§6.2 に理由がある）。

### ビン化の ±1 窓は 2 つの目的で使い分ける

`core/binning.py` が m/z を最近傍ビンに落とし、照会時は `probe_bins()` で ±1 を探る
（境界を挟んだ一致を取りこぼさないため）。一方 `background.build_background()` は
**索引側のピークを ±1 に展開して df を数える**ので、`lookup_idf()` は**中心ビン 1 行だけ**を
引く。ここで再び ±1 に広げると実質 ±2 になり、重み配分が歪んで並び順が変わる。
この非対称は意図的で、両方の docstring に理由が書いてある。

### ビン幅は決め打ちではなく実測（spec §6.4）

`index/mass_error.py` が同一 (instrument_type, inchikey) の**複数レコード間**での
同一フラグメントの m/z 散らばりを測り、ビン幅を決める。1 スペクトル内の近接ピークを
標本に混ぜてはならない（それは別フラグメント）。だからビルドは 2 パス構成になっている。
採用値は `index_meta.bin_width` に残る。

### スコアは IDF 重み付きカバレッジ（spec §6.8）

`weighted_score = Σ idf(matched) / Σ idf(query_peaks)`。並び順を決めるのはこれだけ
（同点は accession で決定的に）。cosine は併記するが並び替えには使わない。
`band` は一致したプロダクトイオン本数由来で、`PRECURSOR_ONLY`（精密質量しか合っていない）を
`PASS` と同一視してはならない。**スコアを合成しないこと**（spec §5.3）。
df は必ず unique InChIKey で数える（レコード数で数えると MassBank の収載バイアスが重みに乗る）。

### TTL には 2 世代のスキーマがある（spec §3.3）

世代A はピーク値をリテラル直付け、世代B は SIO で具象化したノードに `sio:SIO_000300`。
**世代B の TTL 実物は未入手で、fixture は SPARQL パターンからの合成**（spec §11）。
そのため CLI は世代別に `files / records / peaks` を報告し、`gen_b: N files, 0 peaks` が
一目で分かるようにしてある。ピーク 0 本のレコードは `n_no_peaks` で必ず報告する
（黙って弾くとパーサの失敗が成功に見える）。

## 戻り値の規約（spec §12・厳守）

戻り値は LLM の文脈をそのまま占める。

- JSON は `core/serialization.py` の `json_payload()` で返す。`json.dumps(..., indent=2)` を書かない
- 行が並ぶ一覧は `as_tsv()`（列名 1 回）
- float は `round_floats()` で丸めてから返す
- 全ツールに `structured_output=False`（付けないと content と structuredContent で 2 倍送られる）
- ピーク列などの巨大な中間データは `core/session.py` に保持し、payload から外す
- 可変状態（`config.INDEX_PATH` / `session`）は**モジュール修飾で参照**する。
  `from ... import INDEX_PATH` はスナップショット束縛になり、テストの差し替えが効かない
- 打ち切り・スキップ・除外は必ず件数や理由を payload に出す（`truncated` / `n_skipped_empty` /
  `unknown_spectrum_ids`）。黙って減らさない

エラーは `core/errors.py` の機械可読エンベロープを返す（`missing_index` / `missing_state` /
`incompatible_index`）。使用中の MCP SDK では `isError=true` の経路が structuredContent を
捨てるため、本文 JSON に載せている。

## ツールを足すとき

実体を `tools/` の該当モジュールに書き、`__all__` に載せ、`server.py` で
`mcp.tool(structured_output=False)(...)` として登録する。

## 索引成果物

`/index/` と `*.sqlite` は追跡外（.gitignore の先頭 `/` は必須 — アンカーしないと
`massbank_context/index/` まで無視される）。再配布可能な MassBank の中身を成果物に
焼き込まない方針でもある（spec §4）。
