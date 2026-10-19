# 根拠受け渡し・2群比較の検証記録

- 既存回帰: 実装前331 passed、公開ツール追加後355 passed（network除外）。
- 本体stdio: 初期化、14ツール、バンドル読込、別プロセスで再起動後の復元を検証。
- 合成縦断: 明示ID対応のMSP、既存候補、ローカルSPARQL応答を実MCP経由で接続。
  8特徴量の全ページを取得し、確認済み背景、照合候補、RDF根拠、解釈材料を検証した。
- 合成縦断の実行: `C:/Python314/python.exe -m pytest tests/test_evidence_workflow.py -q` → 1 passed。
- 実PubChem接続: この記録時点では未実施。
- ms-data-parserの新export: 実装・検証中。
- 研究者指定の実2群データ: 未指定。
- 実験結果の会話側LLM回答・研究者レビュー: 未実施。

合成テストとサーバ側のMCP成功は、実データの統合完了や解釈の科学的妥当性を意味しない。
