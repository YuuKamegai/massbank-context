# massbank-context

実験で得た MS/MS スペクトルを MassBank-RDF と照合し、候補化合物とメタデータを
返す MCP サーバ。設計書は
`docs/superpowers/specs/2026-08-28-massbank-context-mcp-design.md`。

## Phase 1 でできること

索引の構築と、バイアス無しのスペクトル照合まで。メタデータ取得（Phase 2 以降）は未実装。

## 索引を作る

```
C:/Python314/python.exe -m massbank_context.index build --ttl-dir <TTLのディレクトリ> --output index/massbank.sqlite
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
C:/Python314/python.exe -m massbank_context.server
```

環境変数 `MASSBANK_CONTEXT_INDEX` で索引のパスを上書きできる。

## テスト

```
C:/Python314/python.exe -m pytest tests -q
```
