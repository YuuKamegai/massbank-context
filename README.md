# massbank-context

実験で得た MS/MS スペクトルを MassBank-RDF と照合し、候補化合物とメタデータを
返す MCP サーバ。設計書は Lipidmix_with_LLM リポジトリの
`docs/superpowers/specs/2026-08-28-massbank-context-mcp-design.md`。

## Phase 1 でできること

索引の構築と、バイアス無しのスペクトル照合まで。メタデータ取得（Phase 2 以降）は未実装。

## 索引を作る

```
C:/Python314/python.exe -m massbank_context.index build --ttl-dir <TTLのディレクトリ> --output index/massbank.sqlite
```

既定では `resolution_class == "high"` かつ `acquisition_class != "other"` かつ
InChIKey を持つレコードだけを取り込む。ビン幅は同一 InChIKey のレコード間の
m/z 散らばりを実測して決まり、`index_meta.bin_width` に記録される。

## サーバを起動する

```
C:/Python314/python.exe -m massbank_context.server
```

環境変数 `MASSBANK_CONTEXT_INDEX` で索引のパスを上書きできる。

## テスト

```
C:/Python314/python.exe -m pytest tests -q
```
