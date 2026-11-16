# massbank-context

実験で得た MS/MS スペクトルを MassBank-RDF と照合し、候補化合物とメタデータを
返す MCP サーバ。設計書は
`docs/superpowers/specs/2026-08-28-massbank-context-mcp-design.md`。

## できること

- **根拠付き2群比較**: `load_evidence_bundle` → `set_analysis_context` →
  `match_feature_spectra` → `resolve_candidate_evidence` → `get_comparison_evidence` →
  `draft_comparison_findings` → `validate_comparison_claims`。既存注釈とMassBank候補を保持し、会話側LLMが解釈する。
  [利用手順](docs/workflows/two-group-interpretation.md)、[入力契約](docs/contracts/evidence-bundle-v2.md)。
  中間成果物は`MASSBANK_CONTEXT_ARTIFACTS`（既定`artifacts/`）へ保存し、再起動後も復元する。
  新exportの上流実装と実MCP経路の検証状況は検証記録を参照する。

- **索引の構築とスペクトル照合**（Phase 1）: MassBank-RDF の TTL からローカル索引を
  作り、実験スペクトルを IDF 重み付きカバレッジで照合する。
- **パスウェイ活性同定**（Phase 2〜3。設計書は
  `docs/superpowers/specs/2026-08-31-pathway-activity-design.md`）:
  ms-data-parser が書き出した差次的結果 TSV を読み、InChIKey から PubChem の
  パスウェイを引いて、パスウェイ単位の濃縮（片側 Fisher + BH-FDR）と方向
  （記述 + Wilcoxon 符号順位）を返す。ツールは `load_differential` /
  `resolve_pathways` / `pathway_activity` / `draft_pathway_findings` の 4 本。

- **外来の分子の走査**（設計書は
  `docs/superpowers/specs/2026-11-04-xenobiotic-screen-design.md`）:
  `register_xenobiotics` で薬物などを登録し、`screen_xenobiotic_features` で全特徴量を
  質量関係・MS/MS・検出パターンの 3 つの筋から走査し、`get_xenobiotic_screen` で行を引く。
  全行は TSV に書き、水準は判定ではなく材料として返す。
  spec §8 の較正は r3 で改訂した期待に合格したが、手順書には入っていない（設計書 §10.0）。

MassBank レコードのメタデータ取得（`catalog/` `bias/`）は未実装。

## セットアップ

リポジトリ直下に仮想環境 `.venv` を作り、以降のコマンドはその `python.exe` で実行する。

```
C:/Python314/python.exe -m venv .venv
.venv/Scripts/python.exe -m pip install -r requirements-dev.txt
```

`requirements.txt` は実行に要るもの（`rdflib` と `mcp`）、`requirements-dev.txt` は
それにテスト用の `pytest` を足したもの。`mcp` は `<2` に固定している（2.x は
`FastMCP` を `MCPServer` に改名しており、このサーバは起動しない）。

## ⚠ 外部への通信

`resolve_pathways` は **rdfportal.org（PubChem SPARQL）へ外向きの HTTP を出す**。
このサーバで外部と通信するのはこの経路だけで、索引ビルドとスペクトル照合は
完全にローカルで完結する。

送るのは差次的結果に含まれる **InChIKey** だけで、50 件ずつのバッチにして
SPARQL クエリの `VALUES` 節に入れる。試料名・群名・強度値・ファイルパスは送らない。

結果はローカルの SQLite にキャッシュされ、2 回目以降は照会しない
（日数ベースの TTL は持たない。日付でキャッシュが切れると、同じ解析が日によって
違う結果を出すため）。引き直したいときは `resolve_pathways(refresh=True)`。

エンドポイントに到達できないときは `upstream_unavailable` を返す。最初のバッチが
全滅した時点で残りを照会せずに打ち切る（打ち切り位置は `n_processed` /
`n_remaining` として返る）。

## 索引を作る

```
.venv/Scripts/python.exe -m massbank_context.index build --ttl-dir <TTLのディレクトリ> --output index/massbank.sqlite
```

既定では `resolution_class == "high"` かつ `acquisition_class != "other"` かつ
InChIKey を持ち、ピークを 1 本以上持つレコードだけを取り込む。ビン幅は同一
InChIKey・同一装置種別のレコード**間**の m/z 散らばりを実測して決まり、
`index_meta.bin_width` に記録される。

### 既定のビルドはコーパス最大のコントリビュータを設計上除外する

`Fac_Eng_Univ_Tokyo` は全 54,183 レコード中 12,378 件（23%）を占める最大の
コントリビュータだが、既定のビルドには**入らない**。同群は GC-EI の磁場型
（`instrument_type` が `EI-B` 等）で、spec §6.2 の導出表では
`resolution_class` が `high` にも `low` にもならず `unknown` に落ちるためである。

これは取りこぼしではなく設計上の選択である（spec §6.1 / §6.2）。GC-EI の
フラグメンテーションは LC-ESI-MS2 の CID プロダクトイオンとは別物で、EI には
precursor という概念自体が無い。測定原理の違うスペクトルを同じ背景分布で
数えることに意味は無いので、索引は `acquisition_class` で分割し、別層として
後から追加できる形にしてある。

ビルド後は `index_status()` の `acquisition_classes` と `n_records`、および
CLI が INFO で出す世代別内訳（`gen_a: N files, M records, K peaks`）で、
実際に何が入ったかを確認できる。

## サーバを起動する

```
.venv/Scripts/python.exe -m massbank_context.server
```

### 環境変数

いずれも任意で、既定値の上書き用。

| 変数 | 既定値 | 用途 |
|---|---|---|
| `MASSBANK_CONTEXT_INDEX` | `index/massbank.sqlite` | スペクトル索引のパス |
| `MASSBANK_CONTEXT_PATHWAY_CACHE` | `index/pathway_cache.sqlite` | パスウェイ照会のキャッシュ。索引とは別ファイルなので、索引を再ビルドしても消えない |
| `MASSBANK_CONTEXT_PUBCHEM_SPARQL` | `https://rdfportal.org/pubchem/sparql` | パスウェイを引く SPARQL エンドポイント。**空にしても通信は止まらない**ので、外部へ出したくないときはサーバを起動しないこと |

キャッシュは `cache_meta.contract_version` を持ち、版が違うファイルは開かずに
`contract_mismatch` を返す。その場合はキャッシュファイルを削除して引き直す。

## テスト

```
.venv/Scripts/python.exe -m pytest tests -q
```
