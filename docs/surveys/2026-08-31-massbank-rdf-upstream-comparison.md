# MassBank-RDF 上流リポジトリ比較調査

調査日: 2026-08-31  
対象:

- 既存データ基盤: <https://github.com/systemsomicslab/MassBank-RDF>
- 新解析アプリ: <https://github.com/MotohiroOGAWA/MassbankRdf>

## 結論

`MotohiroOGAWA/MassbankRdf` は `systemsomicslab/MassBank-RDF` の fork や
通常の後継ブランチではない。両者に Git の共通祖先はなく、前者は後者の
MassBank-RDF データを利用するために別履歴で作られた Python 解析アプリである。

既定 `main` だけを見ると、新リポジトリの最終更新は 2026-06-17 であり、
旧リポジトリ `main` の 2026-07-24 より古い。新リポジトリ内で実装が最も
進んでいる公開ブランチは `feature/msp-kg` で、先端は 2026-07-29 の
`a7a0233abacc37a1ff569e123a4f7a8bc3243b74` である。このため本調査とローカル
配置では `feature/msp-kg` を「新実装」として扱う。

## 比較基準

| 対象 | ブランチ | 先端コミット | コミット数 | 追跡ファイル数 |
| --- | --- | --- | ---: | ---: |
| 旧 MassBank-RDF | `main` | `086ca6ab0e934642c712b9e517ba78272c058d22` | 210 | 117,674 |
| 新 MassbankRdf | `main` | `21d3d82810631d48a73268e4ffd5af6169d2fbc7` | 17 | 61 |
| 新 MassbankRdf | `feature/msp-kg` | `a7a0233abacc37a1ff569e123a4f7a8bc3243b74` | 41 | 120 |

新 `main` から `feature/msp-kg` のスナップショット差分は、86ファイル、
15,525行追加、526行削除であった。

## 旧リポジトリの役割

旧リポジトリはデータ生成・配布基盤である。追跡ファイルのうち117,619件が
Turtle (`.ttl`) であり、以下を含む。

- MassBank TXTからTTLを生成するRスクリプト
- 複数世代のRDFモデル、prefix、schema
- QLever設定
- SPARQL例
- RDF検証レポート
- 大規模なMassBank-RDFレコード群

このTTLコーパスが `massbank-context` のローカル索引構築元である。

## 新 `main` で追加された機能

- MassBankメタデータTSVとピークJSONをSQLiteへ取り込むPython importer
- SQLAlchemyモデルとMassBank SQLiteアクセス層
- SQLで計算するコサイン類似度スペクトル検索
- FastAPIにGradioアプリをmountするWeb GUI
- 単一スペクトル検索とMassBank候補表示
- PubChem、HMDB、KNApSAcK向けSPARQL生成・照会
- KG結果のJSON化とダウンロード
- Azure OpenAI Structured OutputsによるKG解釈
- 複数スペクトルの共通ピーク抽出と候補検索
- 一時セッション、入力例、カフェインMSP例

## `feature/msp-kg` で追加・変更された機能

### KGデータ基盤

- PubChem、HMDB、KNApSAcK結果の正規化SQLiteスキーマ
- 1NF版KG SQLiteの生成
- バッチSPARQL照会、タイムアウト時の個別再試行、失敗記録、再開
- KG結果をLLM向けに圧縮したevidence JSON
- full InChIKeyに加え、先頭14文字のshort InChIKey検索

### MSPバッチ解析

- 複数 `.msp` ファイルの一括入力
- ファイル単位のサンプルクラス指定
- 各MSPレコードから極性とprecursor m/zを自動取得
- スペクトルごとのMassBank検索と候補InChIKeyの重複排除
- 25スペクトルごとのチェックポイントと再開
- CSV、JSON、SPARQL、設定、要約をZIPで出力

### 候補ランキングと統計

- 類似度しきい値以下の候補を除外
- MassBank類似度順位とKGメタデータ件数順位のrank sum
- サンプルクラス別存在率、enrichment ratio、Fisher正確検定
- 疾患別Fisher正確検定とBenjamini-Hochberg FDR補正

KGメタデータ件数は「よく注釈された化合物」を優先する指標であり、スペクトル
同定の質そのものではない。`massbank-context` の主ランキングを置換する用途には
使わず、同定後の付加情報として扱う。

### LLM・結果対話

- 化合物由来を `endogenous`、`dietary`、`drug`、`exogenous_other` で出力
- 主張のprovenanceと、試料文脈に対する生物学的妥当性を出力
- 生物学的偽陽性候補と化合物横断要約
- 結果内から決定的に絞った根拠だけをAzure OpenAIへ渡す結果チャット
- 広すぎる質問の拒否、直前のInChIKey範囲を使う追質問

## 修正・表示変更

- KG照会へ渡す一意InChIKey数を設定可能にした
- サンプルの上限を3件から5件へ変更した
- SPARQL結果の不要な直積を減らすようクエリを集約した
- 共通ピーク表から内部列を隠した
- 共通ピーク由来のMassBankヒット表からcosine scoreを隠し、match数中心の表示へ変更した
- GUIのボタンID、結果タブ、LLM設定パネルを整理した
- テストコードは84個の `test_` メソッドまで拡大した

## 技術上の注意

### 新実装は旧TTL基盤を置換しない

新リポジトリには約11.7万件のTTL、TXTからTTLへの変換、QLever設定が含まれない。
MassBank SQLite importerはTSVとピークJSONを要求し、旧TTLから直接SQLiteを作る
処理ではない。新アプリを動かすにも、別途作成済みの `massbank.sqlite3` が必要である。

### スペクトル照合は1対1対応になっていない

新実装のSQLは候補を `id, query_peak_index` でpartitionし、各query peakに最も
近いreference peakを選ぶ。同じreference peakが別のquery peakから再利用される
ことを禁止していない。この多対一対応はdot productを水増しし、入力によっては
cosineが1を超える。`massbank-context` で修正済みの1対1対応を維持し、新実装の
照合SQLで置き換えない。

### Common Peak Annotationは厳密なピーク帰属ではない

候補表は選択した全共通ピークと全MassBankヒットを組み合わせる。参照スペクトルの
どのピークへ実際に一致したかを保持しないため、現段階では候補提示であり、個々の
共通ピークの厳密な注釈ではない。

### 配布・実行情報に欠落がある

`feature/msp-kg` のREADMEが参照する `env/requirements.txt`、Dockerfile、
Composeファイル、起動スクリプトは公開ツリーに存在しない。ルートにも
`requirements.txt`、`pyproject.toml`、`setup.py`、`LICENSE` がない。
HMDBとKNApSAcKの既定endpointはプライベートIPであり、研究室ネットワーク外では
設定変更が必要である。

Pythonの構文コンパイルは成功したが、依存定義が欠落しているため、現在のWindows
環境ではSQLAlchemy、SPARQLWrapper、Gradioを再現可能な形で導入できず、上流の
全テスト合格までは確認できなかった。また `tests/run_tests.py` は結果を終了コードへ
反映しないため、テスト失敗時にも終了コード0となる。

## `massbank-context` での採用方針

維持するもの:

- TTL 2世代対応とprovenance
- 測定方式・分解能・極性で層別した索引
- IDF背景とビン幅の母集団制御
- 1対1ピーク対応
- MCPツールと機械可読エラー契約

後段機能として移植を検討できるもの:

- PubChem、HMDB、KNApSAcKクエリ
- short InChIKeyとKG evidence JSON
- MSPバッチ、チェックポイント、ZIP出力
- サンプルクラス・疾患統計
- GUIと根拠を限定した結果チャット

## ローカル配置

外部リポジトリは親Gitへ取り込まず、独立したnested cloneとして次に配置する。

```text
massbank-context/
├── MassBank-RDF/   # systemsomicslab、TTLデータ基盤
└── MassbankRdf/    # MotohiroOGAWA、feature/msp-kg
```

親 `.gitignore` は両ディレクトリをルートアンカー付きで除外する。これにより、
上流の約11.7万ファイルや独立Git履歴を誤って親コミットへ混入させない。

## 配置後の検証結果

2026-08-31 に移動・クローン・パス修正後の再検証を実施した。

- 旧移動元 `C:/Users/yuu18/MassBank-RDF` は存在せず、
  `C:/Users/yuu18/massbank-context/MassBank-RDF` とTTLディレクトリは存在する
- 旧クローンは `main`、作業ツリーclean、ローカルHEADとリモート `main` が
  `086ca6ab0e934642c712b9e517ba78272c058d22` で一致する
- 新クローンは `feature/msp-kg`、作業ツリーclean、ローカルHEADとリモート同ブランチが
  `a7a0233abacc37a1ff569e123a4f7a8bc3243b74` で一致する
- 親Gitでは両nested cloneが `.gitignore` により除外される
- 親リポジトリの実行コード・設定（両nested cloneと、この移動記録を除く）に旧直下パス
  `C:/Users/yuu18/MassBank-RDF` の残存参照は0件
- `python -m pytest tests -q`: **198 passed**
- `python -m compileall -q massbank_context scripts`: 成功
- `python -m massbank_context.index --help`: 成功
- `git diff --check`: エラー0件

移動に伴う修正対象はREADME・実装計画・試行計画の文書参照であり、Python実装と
MCP設定には旧直下パスのハードコードはなかった。

## 参照

- 旧リポジトリ: <https://github.com/systemsomicslab/MassBank-RDF>
- 旧先端コミット: <https://github.com/systemsomicslab/MassBank-RDF/commit/086ca6ab0e934642c712b9e517ba78272c058d22>
- 新 `main` 履歴: <https://github.com/MotohiroOGAWA/MassbankRdf/commits/main/>
- 新 `feature/msp-kg`: <https://github.com/MotohiroOGAWA/MassbankRdf/tree/feature/msp-kg>
- 新開発先端: <https://github.com/MotohiroOGAWA/MassbankRdf/commit/a7a0233abacc37a1ff569e123a4f7a8bc3243b74>
