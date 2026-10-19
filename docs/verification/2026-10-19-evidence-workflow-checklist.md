# 根拠受け渡し・2群比較の検証記録

- 既存回帰: 実装前331 passed。最新の本体テスト件数は下記の検証結果記録を参照（network除外）。
- 本体stdio: 初期化、14ツール、バンドル読込、別プロセスで再起動後の復元を検証。
- 合成縦断: 明示ID対応のMSP、既存候補、ローカルSPARQL応答を実MCP経由で接続。
  8特徴量の全ページを取得し、確認済み背景、照合候補、RDF根拠、解釈材料を検証した。
- 合成縦断の実行: `C:/Python314/python.exe -m pytest tests/test_evidence_workflow.py -q` → 1 passed。
- 実PubChem接続: 成功。照会結果と非対応生物種を区別し、適合背景がない経路検定を抑止。
- ms-data-parserの新export: 上流コミット `6cf1545`。ARFとDatasetStateのv2 exportを追加。
  上流全テスト3443 passed / 4 skipped / 1896 subtests passed。既存の数値警告1件あり。
- 研究者指定の実2群データ: 2極性を別々に実MCPで解析し、本体MCPへ受け渡した。
  入力・群割当・前処理・閾値は研究者の確認後に固定。修正後も同じ条件で再検証した。
- 実験結果の会話側LLM回答: ローカル成果物に解釈案と構造化した主張を保存。
  主張のfeature/candidate/evidence/context参照に不整合なし。研究者の科学的レビューは未実施。

実データの入力・背景・数値・回答原本はGit追跡外の`artifacts/kidney-aging-validation/`へ保持する。
研究データをこの検証記録へ複写しない。

合成テストとサーバ側のMCP成功は、実データの統合完了や解釈の科学的妥当性を意味しない。
