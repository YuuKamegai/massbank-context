# 根拠受け渡しと2群比較の解釈 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 既存注釈とMassBank照合候補を2群比較に接続し、会話側LLMが出典と不確実性を追って解釈できるようにする。

**Architecture:** ms-data-parserの新しいv2ファイルexportを境界にする。massbank-contextでは不変の解析入力と版付き根拠を保存し、既存照合・統計を再利用する。生物学的評価は会話側LLMで行い、スペクトル順位には反映しない。

**Tech Stack:** Python 3.14、stdlib（dataclasses/csv/json/hashlib/sqlite3/pathlib）、既存rdflib・mcp・pytest。新規の本番依存は追加しない。

**Spec:** [承認済み設計](../specs/2026-10-19-evidence-two-group-interpretation-design.md)


**実行状況（2026-10-19）:** Task 1〜10は実装・検証済み。Task 11は実MCP・実データ・会話側LLMの検証を実施し、研究者による科学的レビューを残す。詳細は[検証記録](../../verification/2026-10-19-two-group-interpretation-results.md)。

## Global Constraints

- 初期開発では、根拠の受け渡しと**2群比較**を一緒に扱う。
- スペクトル順位は保持し、生物学的妥当性は別軸で評価する。LLMやKG情報で並べ替えない。
- 解釈はMCPを利用する会話側LLMが行う。サーバ内でLLM APIを呼ばない。
- `contract_name=massbank-context-evidence`, `contract_version=2`。
- `log2fc_sign` は `positive means group_b is higher` を必須とし、不一致は拒否する。
- 詳細取得は既定50特徴量・最大200、候補照合top_kは既定5・最大20とする。
- 既存の8ツールは維持する。新ツールは内部の型付きサービスを利用し、公開TSVを再パースしない。
- 生の質量データの分析はms-data-parserを使い、本サーバに再実装しない。
- 解析成果物・入力・確認済み背景はGit追跡外とし、研究データを設計書やfixtureに複写しない。
- 外部送信はPubChemへ照会する構造IDに限る。試料背景、群名、強度、ファイルパスは送信しない。
- コード・コメント・docstring・コミットは日本語。既存のモジュール依存方向を維持する。

## Review Focus

1. 同じspot IDを持つ別解析の入力を混同しない。Task 1/2で解析IDとハッシュを検証する。
2. 上位候補だけが同じ構造でも、打切り先に異なる候補がある。Task 3/6で打切りを統計除外する。
3. curationのwrong/redundant/assignで母集団や根拠が消える。Task 10で特徴量を保持し、元注釈と修正を記録する。
4. refresh失敗後に古い成功結果を現在の結果として使う。Task 4/5で失敗版と過去版を区別する。
5. Windowsの日本語パス、UTF-8 BOM、タブを含む名前、再起動・同時更新。Task 1/2/8で往復と版整合性を検証する。

## 実行環境と完了境界

- 本リポジトリ: `C:/Users/yuu18/massbank-context`。Python: `C:/Python314/python.exe`。
- 上流の現在の配置: `C:/Users/yuu18/Metabolomix_with_LLM`、originは `systemsomicslab/Metabolomix_with_LLM`。
- 上流確認箇所: `metabolomix/arf/tools.py::arf_export_differential(output_path, apply_curation=True)`、`metabolomix/arf/identity_join.py::join_identity`、`metabolomix/analysis/dataset_export.py::export_dataset_result`、`metabolomix/analysis/export_contract.py`。
- `.mcp.json`はms-data-parserを上流venvから起動するが、計画作成時の利用可能ツールに同サーバは公開されていない。Task 11の実検証には接続が必要。
- 上流にはARF経路とDatasetState経路がある。新しい共通writerを使い、旧v1の列・版・出力を変更しない。
- 既存データや実稼働DBをfixtureに使わない。作業開始時に両リポジトリの指示・状態を再確認する。
- Task 1〜9の本体とTask 10の上流exportは同じ到達点の依存部分。上流未対応でもfixture開発はできるが、Task 11未完了なら統合完了と報告しない。

## ファイル配置

`core/evidence_models.py`は不変モデル、`core/evidence_ids.py`は識別子生成のみを扱う。
`experiment/bundle.py`は入力契約、`experiment/legacy_bundle.py`はv1の損失を明示する変換。
新しい`evidence/`の`store.py`は永続化、`context.py`は背景、`matching.py`は照合接続、`catalog.py`は索引メタデータを担当する。
`kg/evidence.py`は候補別取得、`kg/evidence_cache.py`は出典別キャッシュを担当する。
`report/comparison_population.py`は統計用射影、`report/comparison.py`は根拠の集約・取得、`report/comparison_findings.py`は解釈材料と引用検証。
`tools/evidence_tools.py`は公開層だけにする。`stats/`にデータモデルやストア依存を持ち込まない。

各Taskは、失敗テスト→最小実装→対象テスト成功→明記したファイルだけをcommitする順に進める。
以下のテストコードで使う共通fixtureはTask 1の`tests/evidence_helpers.py`に定義する。
下流Taskは上流TaskのInterfacesを参照し、公開JSON/TSVを再パースして内部接続しない。

### Task 1: v2契約と型付き入力

**Files:** Create `massbank_context/core/evidence_models.py`, `core/evidence_ids.py`, `experiment/bundle.py`, `experiment/legacy_bundle.py`, `tests/evidence_helpers.py`, `tests/test_evidence_bundle.py`, `docs/contracts/evidence-bundle-v2.md`。

**Interfaces:**
- `load_bundle(path: Path) -> EvidenceBundle`、`import_legacy(path: Path) -> EvidenceBundle`。
- `BundleError(ValueError)`は`code`, `field`, `reason`を持つ。
- `stable_id(kind: str, parts: list[str]) -> str`: UTF-8のJSON配列を区切りなし・キー順固定でSHA-256、`kind:hex`。
- `EvidenceBundle`は`analysis_id`, `manifest`, `features`, `annotations`, `spectrum_links`, `spectra`, `population_complete`を持つ不変dataclass。
- `Feature`は`feature_id`, `spot_id`, `mz`, `rt`, `log2fc`, `p_value`, `q_value`, `mean_a`, `mean_b`, `significant`。欠測数値はNone。
- `Candidate`は`candidate_id`, `feature_id`, `source`, `source_id`, `spectrum_id`, `name`, `inchikey`, `structure_ids`, `rank`, `scores`, `status`, `same_structure_as`, `provenance`, `curation`。statusはactive/rejected/superseded。
- `SpectrumLink`は`feature_id`, `spectrum_file`, `spectrum_id`, `source_analysis_id`, `mapping_evidence`。QuerySpectrumは既存モデルを再利用する。
- `make_bundle(tmp_path, *, case='mixed') -> Path`はmanifestのパスを返す。`mixed`は有意/非有意/欠測/候補競合/名前だけ/MSMSなし/逆符号重複を含む8特徴量。

- [x] **1. 契約を文書化し失敗テストを追加する。** manifestは`source_analysis_id`, `exporter`, `exported_at`, `comparison`, `population`, `files`を必須にする。comparisonは2群とn・方向・上流閾値/方法/前処理、filesはrole/path/sha256。TSV列は上記型に対応し、注釈には`provenance_json`と`curation_json`も持たせる。空の注釈・対応表はヘッダだけで表す。

```python
def test_mixed_bundle_keeps_all_features(tmp_path):
    bundle = load_bundle(make_bundle(tmp_path))
    assert len(bundle.features) == 8
    assert bundle.population_complete is True
    assert any(f.log2fc is None for f in bundle.features)
```

- [x] **2.** `& 'C:/Python314/python.exe' -m pytest tests/test_evidence_bundle.py -q`を実行し、新規モジュール未実装による失敗を確認する。
- [x] **3.** 読込・検証・IDを実装する。TSVは`csv`、UTF-8 BOM許容。manifest正規化内容とファイルハッシュからanalysis_idを生成する。バンドル外の絶対/相対パス・リンク先、NaN/inf、p/q範囲外、重複/欠落ID、別解析の対応、2群でない入力を拒否する。spectrum_idは実ローダ出力との一致を確認し、行順で補わない。
- [x] **4.** 同じコマンドでPASSを確認する。ハッシュ不一致、BOM、日本語パス、引用付きTSV、ゼロ注釈、名前だけ、逆方向、legacy母集団不完全を追加検証する。
- [x] **5.** 上記ファイルをcommit: `feat: 全特徴量と候補を保持する根拠バンドル契約を追加する`。

### Task 2: 再起動可能な解析ストアと背景版

**Files:** Create `massbank_context/evidence/__init__.py`, `evidence/store.py`, `evidence/context.py`, `tests/test_evidence_store.py`, `tests/test_analysis_context.py`; Modify `core/config.py`, `.gitignore`。

**Interfaces:**
- `EvidenceStore(root: Path)`、`register(bundle: EvidenceBundle) -> str`、`load(analysis_id: str) -> EvidenceBundle`。
- `put_snapshot(analysis_id: str, kind: str, payload: dict) -> str`、`get_snapshot(analysis_id: str, revision: str) -> dict`、`latest_revision(analysis_id: str) -> str | None`。
- `set_context(store, analysis_id: str, context: dict, confirmed: bool) -> dict`は`context_revision`, `status`, `context`を返す。

- [x] **1.** 再起動・同一登録・変更後の別ID・書込中断を失敗テストにする。

```python
def test_store_reopens_same_analysis(tmp_path):
    bundle = load_bundle(make_bundle(tmp_path / 'input'))
    root = tmp_path / 'store'
    aid = EvidenceStore(root).register(bundle)
    assert EvidenceStore(root).load(aid).analysis_id == aid
```

- [x] **2.** `python -m pytest tests/test_evidence_store.py tests/test_analysis_context.py -q`を指定Pythonで実行しFAILを確認する。
- [x] **3.** SQLiteトランザクションで不変入力・スナップショット・現在版を保存する。入力のスペクトルファイルもストア内へコピーしハッシュ確認する。設定`MASSBANK_CONTEXT_ARTIFACTS`は既定`<repo>/artifacts`、`/artifacts/`をignore。IDをパスとして解釈しない。背景はspecies/tissue/treatment/groups/measurement/questionの各項目にvalue/source/statusを持ち、推測をconfirmedへ自動昇格しない。
- [x] **4.** 同じ対象テストでPASS。二接続の同時更新、ロールバック、破損ストア、背景更新後の旧確認版保存、未確認・nullを確認する。
- [x] **5.** commit: `feat: 根拠と確認済み背景を版付きで保存する`。

### Task 3: 対応済みスペクトルの照合と候補メタデータ

**Files:** Create `evidence/matching.py`, `evidence/catalog.py`, `tests/test_feature_matching.py`, `tests/test_candidate_catalog.py`; Modify `spectral/matcher.py`。

**Interfaces:**
- 既存`match_spectrum(...) -> list[MatchCandidate]`は維持。
- 新設`match_spectrum_result(conn, spectrum, acquisition_class, precursor_tolerance=0.01, top_k=5) -> MatchResult`。`MatchResult`を`matcher.py`に定義し、`candidates`, `n_total`, `truncated`, `boundary_tie_count`を持たせる。
- `match_features(store, analysis_id: str, feature_ids: list[str] | None, *, top_k=5, acquisition_class='LC-ESI-MS2', precursor_tolerance=0.01) -> dict`。
- `read_candidate_catalog(conn, *, accessions: list[str], inchikeys: list[str]) -> list[dict]`。

- [x] **1.** 新旧照合の同順位・同点・top-k後方の別構造と、背景を変えても順位不変をテストする。

```python
def test_tail_candidates_are_not_unique(conn, query):
    result = match_spectrum_result(conn, query, 'LC-ESI-MS2', top_k=1)
    assert result.truncated is True
    assert result.n_total > len(result.candidates)
```

- [x] **2.** `python -m pytest tests/test_feature_matching.py tests/test_candidate_catalog.py -q`でFAILを確認する。
- [x] **3.** 既存の全候補計算からスライス前情報を取り出す内部関数を共有する。スコアや並び規則を変更しない。順位はスペクトル・出典単位で保持し、複数スペクトル間に合成順位を作らない。索引の内容指紋・meta・条件を根拠版へ記録する。レコード、compounds、compound_classの既存列を使い、取得不能な元TTL参照はnull。候補名・分類・ライセンスを取得する。
- [x] **4.** 対象テストと`tests/test_matcher.py tests/test_tools.py`でPASS。MS/MSなし、空ピーク、索引なし、互換層なし、範囲外top_k、途中の未知feature_idを確認する。
- [x] **5.** commit: `feat: 特徴量と照合候補の根拠を接続する`。

### Task 4: RDF根拠の出典と取得スナップショット

**Files:** Create `kg/evidence.py`, `kg/evidence_cache.py`, `tests/test_candidate_evidence.py`, `tests/test_evidence_cache.py`; Modify `kg/pathway.py`。

**Interfaces:**
- `fetch_pathway_evidence(adapter, inchikeys: list[str]) -> dict[str, dict]`は各keyのstatus/hits/query_version/endpoint/fetched_atを返す。
- `resolve_candidates(store, analysis_id: str, candidate_ids: list[str] | None, *, refresh=False, adapter=None) -> dict`。
- `EvidenceCache(path: Path).lookup(endpoint: str, query_version: str, inchikey: str) -> dict | None`、`.put(snapshot: dict) -> None`。

- [x] **1.** 同じ構造の別候補へ根拠参照を付けること、empty/failed/unsupportedを区別することを失敗テストにする。

```python
def test_same_key_different_endpoints_do_not_share_hits(cache):
    cache.put(snapshot(endpoint='https://one.example/sparql'))
    assert cache.lookup('https://two.example/sparql', '1', KEY_A) is None
```

- [x] **2.** 対象2ファイルを`pytest -q`で実行しFAILを確認する。
- [x] **3.** 既存クエリビルダとtransportを共有し、full pathway URI・predicate URI・taxidを新経路で保持する。既存fetchの戻り値は維持する。新キャッシュはendpoint/query_version/key/取得版をキーにし、旧キャッシュはendpoint一致がmeta/logで確認できる場合だけ来歴不足を付けて移入する。曖昧な既存キャッシュは再照会する。50件バッチ・全面失敗時打切りを再利用し、背景等をHTTPへ入れない。
- [x] **4.** PASS確認。refresh失敗時に旧成功を現在版へ転用しないこと、種不明・別種、index外の既存候補、名前のみ、非公開の入力値がHTTPにないことを確認する。
- [x] **5.** commit: `feat: 候補別のRDF根拠と照会履歴を保持する`。

### Task 5: 版付きの概要・詳細取得

**Files:** Create `report/comparison.py`, `tests/test_comparison_evidence.py`。

**Interfaces:** `get_evidence(store, analysis_id: str, *, feature_ids: list[str] | None=None, cursor: str | None=None, limit: int=50) -> dict`。
返却はanalysis_id/evidence_revision/context_revision/summary/features/candidates/evidence/lookups/truncated/next_cursor。summaryは全件、features等は選択ページ。

- [x] **1.** 上限・ページ合計・全件概要・参照先の同梱を失敗テストにする。

```python
def test_page_does_not_change_population(comparison_store, analysis_id):
    page = get_evidence(comparison_store, analysis_id, limit=1)
    assert page['summary']['n_features_total'] == 8
    assert len(page['features']) == 1
    assert page['next_cursor'] is not None
```

- [x] **2.** `python -m pytest tests/test_comparison_evidence.py -q`でFAILを確認する。
- [x] **3.** feature_id辞書順、明示filterの正規化を固定する。cursorは解析ID・根拠版・背景版・filter hash・offsetを含み検証する。現行版変更後は`stale_cursor`を返し、新旧を混ぜない。被覆は各分母・未照会・失敗件数を分ける。limitは1〜200、候補・RDF詳細はページ内特徴量に限定する。
- [x] **4.** PASS確認。空filterと省略の違い、未知ID、改ざんcursor、途中のrefresh、古い版の根拠ID、欠測、引用付き文字列を確認する。
- [x] **5.** commit: `feat: 比較結果と候補根拠を版付きで取得できるようにする`。

### Task 6: 母集団と保守的なパスウェイ集計

**Files:** Create `report/comparison_population.py`, `tests/test_comparison_population.py`, `tests/test_comparison_pathways.py`。

**Interfaces:** `build_population(bundle: EvidenceBundle, candidates: list[Candidate], match_status: dict) -> PopulationProjection`、`compute_pathways(projection: PopulationProjection, lookups: list[dict], organism: str, min_members: int=3) -> dict`。
`PopulationProjection`を`comparison_population.py`に定義し、compounds/exclusions/representatives/summaryを持たせる。全特徴量の採用・除外理由を返す。

- [x] **1.** 競合/打切り/名前のみ/比較値欠測/逆符号重複、非有意背景をテストする。

```python
def test_opposite_sign_spots_are_excluded(population_case):
    result = build_population(**population_case)
    assert result.exclusions['f-opposite']['reason'] == 'direction_conflict'
    assert all(c['inchikey'] != KEY_CONFLICT for c in result.compounds)
```

- [x] **2.** 対象2ファイルを実行しFAIL確認。
- [x] **3.** 設計§7の規則を実装する。必要な比較値はlog2fc/p/q/mean_a/mean_b、欠測はいずれも除外。候補の未照合状態とMS/MSなしを分け、前者を検定準備未完了とする。採用可能な構造候補は形式検証済full InChIKey、InChI/SMILESのみは自動変換しない。重複代表の最終同点はfeature_id辞書順。ゼロlog2fcは逆符号としない。既存`pathway_activity_rows`へ型付きで渡し、来歴を横持ちする。
- [x] **4.** PASS確認。母集団不完全、採用keyの未照会/失敗、種不一致、経路ゼロ、最小メンバー未満では理由付きで検定を抑止する。旧v1のテストも成功させる。
- [x] **5.** commit: `feat: 候補の曖昧さを保持して2群の経路集計を行う`。

### Task 7: 解釈材料と引用の検証

**Files:** Create `report/comparison_findings.py`, `tests/test_comparison_findings.py`, `docs/workflows/two-group-interpretation.md`。

**Interfaces:** `build_comparison_findings(store, analysis_id: str, organism: str='', include_pathways: bool=False) -> dict`、`validate_claim_references(claims: list[dict], evidence: dict) -> list[dict]`。
claimsはclaim_id/text/kind/feature_ids/candidate_ids/evidence_ids/context_revisionを持つ。

- [x] **1.** 根拠参照、未確認背景、任意経路失敗を失敗テストにする。

```python
def test_unknown_evidence_is_not_a_valid_citation(evidence_page):
    issues = validate_claim_references([claim(evidence_ids=['missing'])], evidence_page)
    assert issues[0]['code'] == 'unknown_evidence'
```

- [x] **2.** `python -m pytest tests/test_comparison_findings.py -q`でFAIL確認。
- [x] **3.** facts/hypotheses/limits/readiness/next_actionsを返す。観測値に活性断定を追加しない。候補依存仮説の空スロットを作り、非有意を同等と書かない。参照が存在しても別feature/candidate/contextなら検証エラー。手順書に会話例、確認方法、model_knowledgeの区分、科学的妥当性は人が評価することを記載する。
- [x] **4.** PASS確認。species未確認でも事実は返し、species依存の経路集計は抑止。RDF失敗でも群間事実を返す。MS/MS未検証と独立検証済を候補ごとに区別する。
- [x] **5.** commit: `feat: 根拠を引用できる2群比較の解釈材料を作る`。

### Task 8: 6ツールの公開とMCP回帰

**Files:** Create `tools/evidence_tools.py`, `tests/test_evidence_tools.py`, `tests/test_evidence_stdio.py`; Modify `server.py`, `tests/test_server_registration.py`, `README.md`。

**Interfaces:** 設計§8の6ツール名を公開する。path/analysis_id等は設計名どおり。getのlimit=50、matchのtop_k=5、refresh=False、include_pathways=False。matchはacquisition_class='LC-ESI-MS2'/precursor_tolerance=0.01も任意引数として公開する。

- [x] **1.** 14ツールの登録・structured_output=False、未知解析ID・エラー復旧案をテストする。

```python
def test_tool_schema_keeps_output_single_copy(registered_tools):
    assert len(registered_tools) == 14
    assert all(t.outputSchema is None for t in registered_tools)
```

- [x] **2.** 新規テストと`tests/test_server_registration.py`を実行しFAIL確認。
- [x] **3.** モジュール修飾のconfigでストアを生成し、Task 1〜7のサービスへ委譲する。docstringへ呼出順・範囲・外部送信を明記し、既存JSONエラー規約へ変換する。invalid入力で既存analysisの現在版を動かさない。
- [x] **4.** `python -m pytest tests/test_evidence_tools.py tests/test_evidence_stdio.py tests/test_server_registration.py -q`でPASS。stdioで初期化・list_tools・load・get・再起動を実際に通す。`python -m pytest tests -q`で全既存回帰を確認する。
- [x] **5.** commit: `feat: 根拠付き2群比較のMCPツールを公開する`。

### Task 9: 固定した合成例でサーバ側の縦断検証

**Files:** Create `tests/test_evidence_workflow.py`, `docs/verification/2026-10-19-evidence-workflow-checklist.md`。

**Interfaces:** Task 8の公開ツールを実MCP stdio経由で利用する。ネットワークだけfixture transportにする。中間JSONを直接改変しない。

- [x] **1.** 背景確認→読込→照合→RDF→ページ→経路→解釈材料のテストを書く。要約8件と全feature ID・candidate ID・evidence IDがつながることをassertする。
- [x] **2.** `python -m pytest tests/test_evidence_workflow.py -q`で不足がある場合の失敗を確認する。既存の組合せだけで成功するなら、その証拠を残し不要な製品コードを追加しない。
- [x] **3.** 縦断検証で見つかった接続部分だけを修正する。確認表にはコード/stdio/実上流/実LLMの証拠を別欄にする。
- [x] **4.** PASSと全テスト成功を確認。input変更・再起動・refresh失敗・背景改訂後に再実行し、旧根拠を現在版と混同しないことを確認する。
- [x] **5.** commit: `test: 根拠と2群比較のMCP縦断経路を検証する`。

### Task 10: ms-data-parserのv2 exporter

**Files（上流リポジトリ相対）:** Create `metabolomix/analysis/evidence_export.py`, `metabolomix/arf/evidence_export.py`, `tests/test_evidence_export.py`, `tests/test_arf_evidence_export.py`; Modify `metabolomix/arf/tools.py`, `metabolomix/tools/dataset_analysis_tools.py`。後者は現在`server.py`からstar importされるため、新ツールを`__all__`にも追加する。

**Interfaces:**
- 共通`write_evidence_bundle(*, comparison: dict, features: list[dict], annotations: list[dict], spectrum_links: list[dict], spectra_files: list[Path], output_dir: Path) -> dict`。
- MCP `arf_export_evidence_bundle(output_dir: str, apply_curation: bool=True) -> str`、`dataset_export_evidence_bundle(output_dir: str, result_id: str, apply_curation: bool=True) -> str`。現行の`session_state.session.dataset`と`ds.results[result_id]`を使い、指定結果がなければ最新結果へ暗黙に置換しない。

- [x] **1.** 実装前に上流AGENTS/CLAUDE・登録箇所・状態・テスト実行Pythonを確認する。ARF兄弟一致とDatasetStateの`assert_current`を再利用する。スペクトルID対応の実フィールドがない経路では空linksと未提供理由を出し、推測接続しない。
- [x] **2.** 全8特徴量を出すfixtureとwrong/suspect/redundant/assignを作り、既存exportでは満たせないテストを実行してFAILを確認する。v1が変わらないテストも併記する。
- [x] **3.** 新exportを追加する。wrongは候補をrejectedにし、assignは元候補をsupersededとして残し新候補を追加、suspectは注記、redundantは特徴量を残して上流の解析対象/除外理由を記録する。rejected/supersededは統計候補集合から除くが詳細で閲覧可能。欠損curation履歴は黙って適用しない。探索専用・古い差次的結果は既存のガードで拒否。原子的にバンドルを保存し、既存v1を書き換えない。
- [x] **4.** 上流の対象テストと既存`test_export_contract.py`, `test_export_differential.py`でPASS。共通writer出力を別プロセスのmassbank-context `load_bundle`へ渡し、全件・群方向・候補来歴が一致することを確認する。
- [x] **5.** 上流だけでcommit: `feat: 全特徴量と候補根拠をv2バンドルへ書き出す`。本体とのコミット対応を確認表へ記録する。

### Task 11: 実データ・実MCP・会話側LLMで検証

**Files:** Modify `docs/verification/2026-10-19-evidence-workflow-checklist.md`; Create `docs/verification/2026-10-19-two-group-interpretation-results.md`（公開可能な要約のみ）。入力・根拠・回答原本はGit追跡外のartifactsへ保存する。

**Interfaces:** ms-data-parserの実MCPで解析/新export、本体の実MCPでTask 8の6ツール。会話側LLMへ同じ根拠版を渡す。

- [x] **1.** 接続できるms-data-parserのツール一覧と実装コミットを確認する。接続なしなら本Taskは未実行と記録し、ローカル関数呼出をMCP成功と報告しない。
- [x] **2.** 公開または研究者指定の2群データ・群割当・前処理・閾値・候補取得条件・評価対象を固定し、ハッシュを保存する。データ未指定なら研究者に対象を確認する。架空の実例や都合のよい候補へ差し替えない。
- [x] **3.** 上流MCPで解析とexport、本体MCPで読込から解釈材料まで実行する。全件数、群方向、ID対応、ページ網羅、RDF被覆、採用/除外を突き合わせる。
- [x] **4.** 会話側LLMへ観測事実/候補根拠/解釈仮説/限界の回答を求め、原文・モデル・取得可能な設定・根拠版を保存する。`validate_claim_references`で参照IDと対応を検証する。
- [ ] **5.** 研究者が実際の主張と科学的妥当性をレビューする。群反転、候補混同、架空引用、失敗を不存在とする説明、未支持の活性断定があれば不合格。修正後は同じ固定条件で再評価する。
- [x] **6.** `python -m pytest tests -q`、`git diff --check`を実行し、実MCP/外部RDF/LLM/研究者レビューの実施・未実施を明記する。成功したテスト件数だけで統合完了にしない。
- [x] **7.** 公開可能な検証記録だけをcommit: `docs: 2群比較の根拠付き解釈を実経路で検証する`。

## 自己レビューと引継ぎ

- 設計§1〜4: Task 7/8/11。§5: Task 1/10。§6: Task 2〜5。§7: Task 6。§8: Task 2/5/8。§9: Task 7/11。§10: Task 1〜5/8。§11: Task 9/11。
- coreのモデル・IDには上位モジュールをimportしない。statsは既存の純関数のまま使う。
- Task 10の上流は別Gitであり、本体へソースやDBをコピーしない。上流の研究者データを外部公開しない。
- Task 10のcuration状態をTask 1のCandidate.statusへ正規化し、Task 6はactive候補だけを採用候補集合にする。
- 設計の範囲追加ではなく、既存コードに存在するcurationとDatasetStateの鮮度ガードを新exportでも保持するための具体化である。
- 研究者の選択に従い、同一セッションで逐次実装した。独立したコードレビューも実施済み。
- チェックは実行記録を反映する。Task 11の科学的妥当性の研究者レビューは未実施。
