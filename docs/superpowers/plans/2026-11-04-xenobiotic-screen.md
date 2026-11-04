# 外来の分子の走査 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 根拠バンドル（v2 / v3）の全特徴について、登録した外来の分子（薬剤・代謝物・溶媒・麻酔薬・併用薬）との関係を 3 つの筋（質量・MS/MS の共有・出方）で示し、水準 A〜D を付けた TSV と要約を返す MCP ツール 3 本を足す。

**Architecture:** 新しいパッケージ `massbank_context/xenobiotic/`（表・質量・登録・MS/MS・出方・走査の 6 モジュール）が `core/`・`stats/`・`spectral/`・`index/`・`experiment/`・`evidence/` を使い、`tools/xenobiotic_tools.py` から公開する。登録と走査の結果は既存の `EvidenceStore`（`config.ARTIFACTS_PATH/evidence.sqlite`）のスナップショットに置く。v3 の `samples` / `intensities` は `experiment/bundle.py` が読み込み時に検証し、走査のときにストアの `sources`（登録時の入力の複製）から読み直す（解析の payload を 28 万行で膨らませないため）。

**Tech Stack:** Python 3.14（`C:/Python314/python.exe`）、SQLite（索引は読み取り専用）、FastMCP。依存は `rdflib` と `mcp` のみ（numpy / scipy は使わない）。

**Spec:** `docs/superpowers/specs/2026-11-04-xenobiotic-screen-design.md`（v3 の書き出し側は `C:\Users\yuu18\Metabolomix_with_LLM\docs\superpowers\specs\2026-11-04-evidence-bundle-v3-samples-design.md`）

**ユーザーの決定（2026-11-04）: 整理してから `main` から始める。** 実装の前に次を順に行う:
(1) fix ブランチの作業ツリーにある未コミットの変更（2026-10-19 からの約 22 ファイル）の中身を確かめ、ユーザーに報告して承認を得てからコミットする。
(2) ユーザーの承認を得てから fix ブランチを `main` にマージする（push はしない）。
(3) `main` から worktree を切る（superpowers:using-git-worktrees）。
いまの MCP サーバはその作業ツリーで動いているので、(1)(2) の後に `index_status` で動作を確かめる。
下の「実行の起点」は、この整理の後に確かめる前提の記録として残す（整理の後は `51e7d2c` が `main` の祖先になっているはず）。

**実行の起点（必ず確かめる）:** ローカルの `main`（70d5de0）は古く、`experiment/` `evidence/` `tools/evidence_tools.py` を持たない。2026-11-04 の時点でこの計画が前提にするコードは
`fix/bin-width-population-and-one-to-one-matching` の先端（51e7d2c）にだけあり、その作業ツリーには他人の未コミットの変更がある（触らない）。
実行者は最初に `git merge-base --is-ancestor 51e7d2c main` を確かめ、真なら `main` から、偽なら fix ブランチの先端から
worktree を切る（superpowers:using-git-worktrees）。この計画が依存するのは 51e7d2c にある次のインターフェースだけ:
`EvidenceStore.read_view / put_snapshot / current / load`、`tools.evidence_tools._invoke`、`core.serialization.json_payload / as_tsv / round_floats`、
`spectral.matcher.match_spectrum_result / DEFAULT_PRECURSOR_TOLERANCE`、`index.background.lookup_idf`、`index.writer.read_meta`、`core.binning.to_bin / probe_bins`。
起点でこれらの名前と引数が変わっていたら、実装の前に報告する。

## Global Constraints

- テストは `C:/Python314/python.exe -m pytest tests -q`。全タスクの終わりで全件が通ること。
- コード・コメント・docstring・コミットメッセージは日本語。
- `core/` は leaf（上位を import しない）。`stats/` も leaf。`xenobiotic/` は `tools/` を import しない。
- 戻り値は `json_payload()`、行の一覧は `as_tsv()`、float は丸める（`_invoke` 経由）。全ツール `structured_output=False`。
- 可変状態は `config.INDEX_PATH` / `config.ARTIFACTS_PATH` のようにモジュール修飾で読む。
- 索引は `sqlite3.connect(config.INDEX_PATH.resolve().as_uri() + '?mode=ro', uri=True)` で読み取り専用に開く。
- 外向きの通信をしない（`resolve_pathways` の例外を増やさない）。
- 名前は照合キーにしない。InChIKey と組成式だけを使う（spec §2 反目標）。
- 特徴を黙って除かない。TSV は全特徴を含む（関係の無い特徴も `level` 空の 1 行）。打ち切り・除外は件数を出す。
- 単独の筋だけの一致は D に留める（spec §2 反目標・§5）。
- 付加体（陽）: `[M+H]+` `[M+Na]+` `[M+NH4]+` `[M+K]+` `[2M+H]+` `[M+2H]2+`。（陰）: `[M-H]-` `[M+Cl]-` `[M+HCOO]-` `[M+CH3COO]-` `[2M-H]-`。
- イオン源の中の損失: `-H2O` `-NH3` `-CO2`（付加体の後に 1 つまで）。
- 代謝の質量差（15 種、組み合わせは 2 つまで）: +O、+2O、−CH2、−C2H4、−2H、+2H、+H2O、+C6H8O6、+SO3、+C10H15N3O6S、+C2H2O、+CH2、+C2H3NO、+C2H5NO2S、−Cl+OH。
- 質量差の表・付加体の表・閾値は版の文字列を持ち、走査の出力に必ず書く。
- MS/MS が無い特徴は「MS/MS なし」（`no_msms`）と明示し、得点 0 と区別する。
- 索引の極性は大文字（`POSITIVE` / `NEGATIVE`）、v3 の manifest の `comparison.polarity` は小文字（`positive` / `negative`）。境界で変換する。
- MS/MS は `acquisition_class = 'LC-ESI-MS2'` だけを使う（層を跨がない。CLAUDE.md）。

## Review Focus

1. **`detected` が全部空（検出状態の取り込み失敗）のバンドル** — 「未検出」と数えると処理群だけの偽の偏りができる。出方の筋は `unavailable`（理由 `detection_unknown`）になり、`detected_test` / `detected_reference` は空。→ Task 6 のテスト `test_unknown_detection_is_unavailable_not_zero`。
2. **m/z が空の特徴**（v2 は数値の欠測を許す） — 落ちずに行が残り、`notes` に `no_mz`。→ Task 7 `test_feature_without_mz_is_kept`。
3. **登録し直した後の古いカーソル** — `get_xenobiotic_screen` が新旧の行を混ぜず `stale_cursor` で止まる。→ Task 7 `test_cursor_is_stale_after_rescreen`。
4. **塩素を持たない分子に −Cl+OH・炭素の足りない分子に −C2H4** — 負の元素数のイオンを作らず、例外も出さない。→ Task 3 `test_shift_never_yields_negative_counts`。
5. **ピーク 0 本のスペクトル・前駆体 m/z の無いスペクトルに結ばれた特徴** — 前者は `empty_spectrum`、後者は product ion だけで照合し、どちらも「得点 0」と区別して `notes` に出る。→ Task 5 `test_empty_and_precursorless_spectra`。

## 設計との食い違い・設計の補足

実装者が迷わないよう、設計書が決めていない点と、コードと食い違う点をここで決める。較正（Task 9）の後、設計書の追記に転記する。

- **artifacts の仕組み**（spec §3.1）: 既存の仕組みは `EvidenceStore`（`artifacts/evidence.sqlite` の `snapshots` / `heads`）。登録は種別 `xenobiotics`、走査は種別 `xenobiotic_screen` のスナップショットにする。副作用として解析の根拠版（`latest_revision`）が変わり、`get_comparison_evidence` の古いカーソルは `stale_cursor` になる（既存の `set_analysis_context` と同じ振る舞い）。
- **v3 の試料の値の置き場所**: `EvidenceBundle` に足さない。ストアは登録時に manifest の全ファイルを `sources` 表に複製しているので、走査はそこから読み直す（`EvidenceStore.read_source` を足す）。`Feature` には任意の `adduct` を足す（`get_comparison_evidence` の features に `adduct` 列が増える。v2 では空）。
- **`intensities_stage` の位置**: ms-data-parser の設計書は「manifest に書く」としか言わない。manifest の最上位キー `intensities_stage` と決め、契約文書に書く（ms-data-parser 側と合わせる。未決の点として報告）。
- **v3 の厳しさ**: `samples` と `intensities` はどちらか片方だけなら誤り。v2 にこの 2 役割があれば誤り。`samples` の `group` が `a` の数は `comparison.n_a`、`b` は `n_b` と一致すること。`intensities` は全特徴 × `group` が a/b の全試料をちょうど 1 回ずつ持つこと（欠けは黙って補わない）。
- **`ppm_bands` の形**（spec §6 は形を決めていない）: ms-data-parser の `band_summary` と同じ `[{"min_mz": float|null, "max_mz": float|null, "ppm": float}]`、半開区間 `[min_mz, max_mz)`、**理論 m/z** で帯を選ぶ。帯の重なりは誤り。どの帯にも入らなければ `ppm` を使う。
- **`relation` の書式**: `<質量差>|<イオン>`。質量差は無ければ `parent`、あれば表の名前を表の順に `,` で連結（例 `demethylation,oxidation|[M+H]+`）。損失はイオンの `]` の前に入れる（`[M+H-H2O]+`）。MS/MS の共有だけの行は `msms_only`。
- **同じ m/z になるイオンの形**: 付加体・損失・質量差の組み合わせには恒等的に同じ m/z になるもの（`[M+NH4]+` から `-NH3` ＝ `[M+H]+`、`+H2O` の質量差と `-H2O` の損失）がある。畳まずに出すと 1 つの特徴に同じ一致が何行も並ぶので、理論 m/z が 1e-6 Da 未満で一致する形は最も単純なもの 1 つに畳む（Task 3）。それでも化学的に不自然な一致（例: ケタミンの `dechlorination_hydroxylation|[M+K-H2O]+` がヒドロキシノルケタミンの [M+H]+ の 0.3 ppm 以内に来る）は残る。これは質量の筋だけでは D に留まる設計で抑える。
- **1 つの (特徴, entry) に複数の質量の一致**: 一致ごとに行を分け（質量差の数 → |ppm| の順）、`level` 以下の組の値は各行に繰り返す。
- **共有の得点**: 分母は参照の上位ビン（product ion と neutral loss それぞれ相対強度の上位 `msms_top_n` 本）の IDF の和、分子は共有したビンのうち上位ビンに入るものの IDF の和（得点を 0〜1 に収める）。`msms_shared_n` と `msms_shared_formulas` は和集合全体で数える。
- **前駆体の残り**: 参照とクエリの両方で、m/z が前駆体 − 0.5 以上の product ion は共有に数えない（前駆体そのものと同位体。同じ質量の内因性の物質が前駆体のビンだけで共有を得るのを防ぐ。`binning._MIN_NEUTRAL_LOSS` と同じ 0.5）。
- **一対一**: 共有の判定は照合と同じく ±1 ビンを探り、中心ビンを優先して参照のビンを一度だけ使う（`spectral/matcher.py` の `_preference_order` と同じ規則。matcher は fix ブランチで変更中なので import せず、同じ規則を `xenobiotic/msms.py` に書く）。
- **`forward_rank`**: 既存の `match_spectrum_result(conn, spectrum, 'LC-ESI-MS2', DEFAULT_PRECURSOR_TOLERANCE, top_k=20)` の候補を**重複を除いた InChIKey の順**に並べたときの entry の InChIKey の順位（1 始まり）。複数のスペクトルが結ばれていれば最小。上位 20 に無ければ空。
- **水準 D**: 「1 つの筋だけ」は質量か MS/MS の共有のどちらか 1 つ。出方だけ（質量も MS/MS も無い）の組は entry との関係の根拠が無いので行にしない（処理群だけの特徴が全 entry に D で付くのを防ぐ）。
- **親との連動の「親」**: entry に `parent_entry_id` があればその entry、無く質量の一致が質量差つきなら entry 自身（本体の水準 A の特徴と比べる）。親の水準 A の特徴が複数あれば ρ が最大のものを取り、`notes` に `parent_feature:<id>`。ρ は test 群の試料のうち、両方の値があり、どちらの `detected` も `false` でない試料で計算する。
- **索引に参照スペクトルの無い代謝物**（ヒドロキシノルケタミン・デヒドロノルケタミンは索引に 0 件。2026-11-04 確認）: 自身の entry の組は MS/MS が `unavailable` なので最高 C。B には親（ケタミン）の entry の組（質量差 2 つ: `demethylation,oxidation` / `demethylation,dehydrogenation`）と親の参照スペクトルの共有で届く。特徴の水準は組の最大なので設計どおり。
- **合格の条件が届かないかもしれない点**（較正で確かめる）: 索引のフルオキセチン（34 件）・ノルフルオキセチン（18 件）は全部 POSITIVE。ST002901 の「処理群だけに出る 3 件」のうち Autophinib と 4-Bromo-2-{…} は陰イオンの特徴なので、陰のバンドルではフルオキセチンの MS/MS の筋が `no_reference_spectra` になり、質量の一致も無ければ C に届かない。この場合も質量の表を広げず、満たせなかったことを設計書に記録する（spec §8）。
- **索引を開けないとき**: 登録も走査も `missing_index` で止める（登録の戻り値に索引のレコード数を出す設計のため）。

## ファイルの構成

| ファイル | 責務 |
|---|---|
| `massbank_context/core/formula.py`（新） | 組成式の解析・単一同位体質量・Hill 表記・元素の加減。純関数 |
| `massbank_context/core/evidence_models.py`（改） | `Feature.adduct`、`Sample` |
| `massbank_context/experiment/bundle.py`（改） | v2 / v3 の読み込み、`parse_samples` / `parse_intensities` |
| `massbank_context/evidence/store.py`（改） | `read_source` を足すだけ |
| `docs/contracts/evidence-bundle-v3.md`（新） | v3 の契約 |
| `massbank_context/stats/correlation.py`（新） | Spearman の ρ |
| `massbank_context/xenobiotic/__init__.py`（新） | 空（docstring だけ） |
| `massbank_context/xenobiotic/tables.py`（新） | 付加体・損失・質量差・閾値の定数と版、イオンの形の列挙 |
| `massbank_context/xenobiotic/mass.py`（新） | 質量の筋 |
| `massbank_context/xenobiotic/registry.py`（新） | 登録の検証・索引の情報・保存と復元、`XenobioticError` |
| `massbank_context/xenobiotic/msms.py`（新） | 索引を開く、参照ビン、共有の得点、`forward_rank` |
| `massbank_context/xenobiotic/pattern.py`（新） | 群の偏り・親との連動 |
| `massbank_context/xenobiotic/screen.py`（新） | 水準の割り当て、行・TSV・要約・スナップショット、行の取得 |
| `massbank_context/tools/xenobiotic_tools.py`（新） | MCP 公開 3 本 |
| `massbank_context/server.py`（改） | 登録 |
| `tests/xenobiotic_helpers.py`（新） | 合成の v3 バンドルと合成の索引 |

---

### Task 1: 組成式の質量（core/formula.py）

**Files:**
- Create: `massbank_context/core/formula.py`
- Test: `tests/test_formula.py`

**Interfaces:**
- Consumes: なし
- Produces:
  - `class FormulaError(ValueError)`
  - `ELEMENT_MASSES: dict[str, float]`（H, C, N, O, P, S, F, Cl, Br, I, Na, K, Si, B, Se の単一同位体質量。C = 12.0、H = 1.00782503207、N = 14.0030740048、O = 15.99491461956、S = 31.97207100、P = 30.97376163、F = 18.99840322、Cl = 34.96885268、Br = 78.9183371、I = 126.904473、Na = 22.9897692809、K = 38.96370668、Si = 27.9769265325、B = 11.0093054、Se = 79.9165213）
  - `ELECTRON_MASS = 0.00054857990946`
  - `parse_formula(text: str) -> dict[str, int]`
  - `monoisotopic_mass(composition: dict[str, int]) -> float`
  - `hill_formula(composition: dict[str, int]) -> str`
  - `shift_composition(base: dict[str, int], add: dict[str, int], remove: dict[str, int]) -> dict[str, int] | None`（どれかの元素が負になれば `None`。0 の元素は落とす）

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_parse_accepts_msdial_order():
    assert parse_formula('C12ClH14NO2') == {'C': 12, 'Cl': 1, 'H': 14, 'N': 1, 'O': 2}

def test_monoisotopic_mass_metformin():
    assert monoisotopic_mass(parse_formula('C4H11N5')) == pytest.approx(129.101445, abs=1e-6)

def test_hill_order():
    assert hill_formula(parse_formula('C12ClH14NO2')) == 'C12H14ClNO2'
    assert hill_formula(parse_formula('ClNa')) == 'ClNa'          # 炭素なしはアルファベット順

@pytest.mark.parametrize('bad', ['', 'C4H11X', 'C4(H2)', 'C4H11N5+', 'c4h11'])
def test_parse_rejects(bad):
    with pytest.raises(FormulaError):
        parse_formula(bad)

def test_shift_composition_refuses_negative():
    base = parse_formula('C4H11N5')
    assert shift_composition(base, {'O': 1, 'H': 1}, {'Cl': 1}) is None
    assert shift_composition(base, {}, {'C': 1, 'H': 2}) == {'C': 3, 'H': 9, 'N': 5}
```

- [ ] **Step 2: 失敗を確かめる** — Run: `C:/Python314/python.exe -m pytest tests/test_formula.py -q` / Expected: FAIL（`ModuleNotFoundError: massbank_context.core.formula`）
- [ ] **Step 3: 実装する** — 正規表現 `([A-Z][a-z]?)(\d*)` で全体を覆えるか確かめ（覆えなければ `FormulaError`）、同じ元素の重複は足し合わせる。stdlib のみ。
- [ ] **Step 4: 通ることを確かめる** — 同じコマンドで PASS。
- [ ] **Step 5: コミット** — `git add massbank_context/core/formula.py tests/test_formula.py` / `git commit -m "feat: 組成式から単一同位体質量を出す純関数を core に置く"`

---

### Task 2: 根拠バンドル v3 の読み込みと契約

**Files:**
- Modify: `massbank_context/core/evidence_models.py`（`Feature` の末尾に `adduct: str | None = None`、新しい `Sample`）
- Modify: `massbank_context/experiment/bundle.py`（`_load_bundle` の版判定・役割の分岐・features の `adduct`、新しい関数）
- Modify: `massbank_context/evidence/store.py`（`read_source` を足す。他は触らない）
- Modify: `massbank_context/tools/evidence_tools.py`（`load_evidence_bundle` の docstring を「v2/v3 manifest」に）
- Create: `docs/contracts/evidence-bundle-v3.md`
- Create: `tests/xenobiotic_helpers.py`（`make_v3_bundle`）
- Test: `tests/test_evidence_bundle_v3.py`

**Interfaces:**
- Consumes: 既存の `load_bundle`, `BundleError`, `EvidenceStore`
- Produces:
  - `@dataclass(frozen=True) class Sample: sample_id: str; group: str; role: str; injection_order: int | None`（`group` は `'a'` / `'b'` / `''`）
  - `SUPPORTED_CONTRACT_VERSIONS = (2, 3)`、`SAMPLE_COLUMNS = ['sample_id', 'group', 'role', 'injection_order']`、`INTENSITY_COLUMNS = ['feature_id', 'sample_id', 'value', 'detected']`
  - `parse_samples(text: str) -> tuple[Sample, ...]`
  - `parse_intensities(text: str, feature_ids: set[str], samples: Sequence[Sample]) -> dict[str, dict[str, tuple[float | None, bool | None]]]`（feature_id → sample_id → (value, detected)。`detected` の空は `None`）
  - `role_path(manifest: dict, role: str) -> str | None`（manifest の `files` から役割のパス）
  - `EvidenceStore.read_source(self, analysis_id: str, path: str) -> bytes`（無ければ `ValueError('入力ファイルが保存されていません: <path>')`）
  - `tests/xenobiotic_helpers.py`: `make_v3_bundle(root: Path, *, version: int = 3, polarity: str | None = 'positive', detection: str = 'known', with_samples: bool = True) -> Path`（manifest のパスを返す。中身は下記）

`make_v3_bundle` の中身（後のタスクのテストはこの値に依存する）:
- comparison: `group_a='Vehicle'`, `group_b='Drug'`, `n_a=6`, `n_b=6`、v3 なら `polarity`、最上位に `intensities_stage='pre_imputation'`。
- samples: `T1`〜`T6`（group `b`）、`R1`〜`R6`（group `a`）、`QC1`（group 空、role `qc`）。
- features（`adduct` 列つき。全部空）と MS/MS（`spectra.msp`。IONMODE Positive。spectrum_id は `spectra:<n>`）:

| feature_id | mz | MS/MS（precursor; peaks） | 検出（test / reference） | 値（test） |
|---|---|---|---|---|
| `k_parent` | 238.0990 | 238.0990; 125.0153, 207.0571, 220.0888, 238.0990 | 6 / 0 | 100,200,…,600 |
| `k_hnk` | 240.0786 | 240.0786; 125.0153, 89.0386, 222.0680 | 6 / 0 | 10,20,…,60 |
| `k_frag` | 220.0885 | なし | 6 / 0 | 5 |
| `x_msms` | 500.3000 | 500.3000; 125.0153, 207.0571 | 4 / 1 | 7 |
| `mass_only` | 254.0942 | なし | 6 / 6 | 9 |
| `endo` | 300.2000 | なし | 6 / 6 | 9 |
| `no_mz` | 空 | なし | 6 / 6 | 9 |

  reference の値は検出が偽なら 1、真なら 9。`detection='unknown'` なら `detected` を全部空にし、`population.detection_unavailable_reason='sample_axis_mismatch'`。`version=2` なら samples / intensities を書かず `polarity` も書かない。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_v3_loads_and_keeps_adduct(tmp_path):
    bundle = load_bundle(make_v3_bundle(tmp_path))
    assert bundle.manifest['comparison']['polarity'] == 'positive'
    assert {f.feature_id: f.adduct for f in bundle.features}['k_parent'] is None

def test_v2_still_loads(tmp_path):
    bundle = load_bundle(make_v3_bundle(tmp_path, version=2))
    assert role_path(bundle.manifest, 'intensities') is None

def test_parse_intensities_round_trip(tmp_path):
    path = make_v3_bundle(tmp_path)
    manifest = json.loads(path.read_text(encoding='utf-8'))
    samples = parse_samples((tmp_path / role_path(manifest, 'samples')).read_text(encoding='utf-8-sig'))
    matrix = parse_intensities((tmp_path / role_path(manifest, 'intensities')).read_text(encoding='utf-8-sig'),
                               {'k_parent', 'k_hnk', 'k_frag', 'x_msms', 'mass_only', 'endo', 'no_mz'}, samples)
    assert matrix['k_parent']['T3'] == (300.0, True)
    assert matrix['k_parent']['R1'] == (1.0, False)
    assert 'QC1' not in matrix['k_parent']

@pytest.mark.parametrize('breakage', [
    'samples_without_intensities', 'v2_with_samples', 'bad_group', 'n_b_mismatch',
    'missing_intensity_row', 'duplicate_intensity_row', 'unknown_sample', 'bad_detected',
    'bad_polarity', 'missing_polarity_v3', 'bad_intensities_stage', 'contract_version_4'])
def test_v3_contract_violations_raise(tmp_path, breakage):
    path = make_v3_bundle(tmp_path)
    break_bundle(path, breakage)          # このテストファイル内のヘルパ。rewrite_table と change_manifest を使う
    with pytest.raises(BundleError):
        load_bundle(path)

def test_store_keeps_intensities_source(tmp_path):
    store = EvidenceStore(tmp_path / 'store')
    bundle = load_bundle(make_v3_bundle(tmp_path / 'in'))
    aid = store.register(bundle)
    raw = store.read_source(aid, role_path(bundle.manifest, 'intensities'))
    assert raw.decode('utf-8-sig').splitlines()[0].split('\t') == INTENSITY_COLUMNS
```

- [ ] **Step 2: 失敗を確かめる** — Run: `C:/Python314/python.exe -m pytest tests/test_evidence_bundle_v3.py -q` / Expected: FAIL（`ImportError: parse_samples`）
- [ ] **Step 3: 実装する** — `_load_bundle` の版判定を `manifest['contract_version'] in SUPPORTED_CONTRACT_VERSIONS`（int 型）にし、誤りの文言を「対応する版は massbank-context-evidence v2 / v3 です」に。v3 では `comparison.polarity` ∈ {`positive`, `negative`} を必須、v2 にあれば同じ検証。役割 `samples` / `intensities` は v3 のときだけ各 1 個まで許し、両方そろうか両方無いか。`_table(path, columns)` の本体を `_table_text(text, name, columns)` に切り出し、`parse_samples` / `parse_intensities` はそれを使う。読み込み時は両方を呼んで検証し、行列は捨てる（保持しない）。`intensities` があれば最上位の `intensities_stage` ∈ {`pre_imputation`, `post_imputation`} を必須に。`value` は空か有限の数（負も可。対数変換の値があるため）、`detected` は `true` / `false` / 空。features の `adduct` は `row.get('adduct') or None`。`read_source` は `SELECT payload FROM sources WHERE analysis_id=? AND path=?`。
- [ ] **Step 4: 契約文書を書く** — `docs/contracts/evidence-bundle-v3.md`。v2 の文書の書きぶりに合わせ、v2 との差分だけを書く: `contract_version` 3、`comparison.polarity`、最上位の `intensities_stage`、`population.detection_unavailable_reason`（任意）、役割 `samples` / `intensities` の列と規則（上の「v3 の厳しさ」）、features の任意列 `adduct`、massbank-context は 2 と 3 の両方を読むこと、見本のバンドルを ms-data-parser の `tests/fixtures/` と共有すること（写しは Task 9 で置く）。
- [ ] **Step 5: 通ることを確かめる** — Run: `C:/Python314/python.exe -m pytest tests -q` / Expected: 全件 PASS（既存の v2 のテストを含む）
- [ ] **Step 6: コミット** — `git commit -m "feat: 根拠バンドル v3（試料と試料ごとの値）を読み、契約を文書にする"`

---

### Task 3: 付加体・質量差の表と質量の筋

**Files:**
- Create: `massbank_context/xenobiotic/__init__.py`, `massbank_context/xenobiotic/tables.py`, `massbank_context/xenobiotic/mass.py`
- Test: `tests/test_xenobiotic_mass.py`

**Interfaces:**
- Consumes: Task 1 の `parse_formula`, `monoisotopic_mass`, `shift_composition`, `ELECTRON_MASS`；Task 2 の `Feature`（`mz`, `adduct`）
- Produces（`tables.py`）:
  - `ADDUCT_TABLE_VERSION = 'xeno-adducts-1'`、`SHIFT_TABLE_VERSION = 'xeno-shifts-1'`、`THRESHOLDS_VERSION = 'xeno-thresholds-provisional'`
  - `THRESHOLDS = {'msms_shared_min': 0.3, 'msms_top_n': 10, 'bias_reference_max': 1, 'rho_min': 0.6, 'rho_min_n': 5}`（仮。Task 9 で決める）
  - `@dataclass(frozen=True) class Adduct: name: str; polarity: str; n_molecules: int; charge: int; add: str; remove: str`。m/z = (n·M + mass(add) − mass(remove) − charge·e) / |charge|。`ADDUCTS: tuple[Adduct, ...]`（Global Constraints の 11 種。`[M+NH4]+` の add は `NH4`、`[M+HCOO]-` は `CHO2`、`[M+CH3COO]-` は `C2H3O2`、`[M-H]-` と `[2M-H]-` は remove `H`、`[M+2H]2+` は add `H2`・charge 2）
  - `SOURCE_LOSSES = ('H2O', 'NH3', 'CO2')`
  - `@dataclass(frozen=True) class Shift: name: str; add: str; remove: str`。`SHIFTS`（この順）: `oxidation`(+O), `dioxidation`(+O2), `demethylation`(−CH2), `deethylation`(−C2H4), `dehydrogenation`(−H2), `hydrogenation`(+H2), `hydrolysis`(+H2O), `glucuronidation`(+C6H8O6), `sulfation`(+SO3), `glutathione_conjugation`(+C10H15N3O6S), `n_acetylation`(+C2H2O), `methylation`(+CH2), `glycine_conjugation`(+C2H3NO), `taurine_conjugation`(+C2H5NO2S), `dechlorination_hydroxylation`(add `OH`, remove `Cl`)
  - `@dataclass(frozen=True) class IonForm: shifts: tuple[str, ...]; ion: str; mz: float`、property `relation -> str`（書式は「設計との食い違い・補足」）
  - `ion_forms(composition: dict[str, int], polarity: str, max_shifts: int) -> list[IonForm]`（polarity は `positive` / `negative`。質量差は重複なしの組み合わせ。正味の質量差が 1e-6 Da 未満の組（−CH2 と +CH2 など）と負の元素数になる組は作らない。**理論 m/z が 1e-6 Da 未満で一致する形は 1 つに畳む**（`[M+NH4-NH3]+` ≡ `[M+H]+`、`hydrolysis` ＋ `-H2O` ≡ 本体など）: 質量差の数が少ない → 損失なし → 付加体の表の順 → 質量差の表の順で最初のものを残す）
  - `ion_label(adduct: str, loss: str | None) -> str`
- Produces（`mass.py`）:
  - `@dataclass(frozen=True) class MassRelation: feature_id: str; entry_id: str; shifts: tuple[str, ...]; ion: str; theoretical_mz: float; ppm: float`、property `relation`
  - `validate_ppm_bands(ppm_bands: list[dict] | None) -> list[dict]`
  - `tolerance_ppm(theoretical_mz: float, ppm: float, bands: list[dict]) -> float`
  - `find_mass_relations(features: Sequence[Feature], entries: Sequence[XenoEntry], polarity: str, ppm: float, ppm_bands: list[dict], max_shifts: int) -> tuple[dict[tuple[str, str], list[MassRelation]], Counter]`（鍵は (feature_id, entry_id)、各リストは (質量差の数, |ppm|) 順。Counter は `adduct_conflict_dropped` / `adduct_unrecognized` / `no_mz`）。`XenoEntry` は Task 4 で定義するが、ここでは `entry_id` と `composition` 属性だけを使う（テストは `SimpleNamespace` で渡す）
  - ppm = (観測 − 理論) / 理論 × 1e6。特徴の m/z を昇順に並べ、`bisect` で窓を引く。

- [ ] **Step 1: 失敗するテストを書く**

```python
def _mz(forms, relation):
    return next(f.mz for f in forms if f.relation == relation)

def test_metformin_positive_ions():
    forms = ion_forms(parse_formula('C4H11N5'), 'positive', 0)
    assert _mz(forms, 'parent|[M+H]+') == pytest.approx(130.108722, abs=1e-5)
    assert _mz(forms, 'parent|[M+Na]+') == pytest.approx(152.090666, abs=1e-5)
    assert _mz(forms, 'parent|[2M+H]+') == pytest.approx(259.210167, abs=1e-5)
    assert _mz(forms, 'parent|[M+2H]2+') == pytest.approx(65.557999, abs=1e-5)
    assert _mz(forms, 'parent|[M+H-NH3]+') == pytest.approx(113.082173, abs=1e-5)

def test_metformin_negative_ions():
    forms = ion_forms(parse_formula('C4H11N5'), 'negative', 0)
    assert _mz(forms, 'parent|[M-H]-') == pytest.approx(128.094169, abs=1e-5)
    assert _mz(forms, 'parent|[M+Cl]-') == pytest.approx(164.070847, abs=1e-5)
    assert _mz(forms, 'parent|[M+HCOO]-') == pytest.approx(174.099648, abs=1e-5)

def test_two_shift_metabolites_of_ketamine():
    forms = ion_forms(parse_formula('C13H16ClNO'), 'positive', 2)
    assert _mz(forms, 'demethylation,oxidation|[M+H]+') == pytest.approx(240.078583, abs=1e-5)
    assert _mz(forms, 'demethylation,dehydrogenation|[M+H]+') == pytest.approx(222.068018, abs=1e-5)
    assert not any(f.shifts == ('demethylation', 'methylation') for f in forms)

def test_shift_never_yields_negative_counts():
    forms = ion_forms(parse_formula('CH4'), 'positive', 2)          # 塩素なし・炭素 1
    assert not any('dechlorination_hydroxylation' in f.shifts for f in forms)
    assert not any('deethylation' in f.shifts for f in forms)

def test_identical_ion_forms_collapse():
    forms = ion_forms(parse_formula('C13H16ClNO'), 'positive', 2)
    at_parent = [f for f in forms if abs(f.mz - 238.099318) < 1e-4]
    assert [f.relation for f in at_parent] == ['parent|[M+H]+']      # [M+NH4-NH3]+ や hydrolysis|[M+H-H2O]+ は畳まれる

def test_max_shifts_limits_combinations():
    assert all(len(f.shifts) <= 1 for f in ion_forms(parse_formula('C4H11N5'), 'positive', 1))

def test_find_relations_orders_by_shift_count_and_ppm():
    entries = [SimpleNamespace(entry_id='X1', composition=parse_formula('C13H16ClNO'))]
    features = [_feature('k_parent', 238.0990), _feature('k_hnk', 240.0786), _feature('no_mz', None)]
    relations, counts = find_mass_relations(features, entries, 'positive', 5.0, [], 2)
    assert relations[('k_parent', 'X1')][0].relation == 'parent|[M+H]+'
    assert relations[('k_parent', 'X1')][0].ppm == pytest.approx(-1.34, abs=0.05)
    hnk = [r.relation for r in relations[('k_hnk', 'X1')]]
    assert 'demethylation,oxidation|[M+H]+' in hnk
    assert [len(r.shifts) for r in relations[('k_hnk', 'X1')]] == sorted(len(r.shifts) for r in relations[('k_hnk', 'X1')])
    assert counts['no_mz'] == 1

def test_adduct_column_drops_conflicting_relations():
    features = [_feature('f', 260.0813, adduct='[M+H]+')]           # ケタミンの [M+Na]+ の位置
    entries = [SimpleNamespace(entry_id='X1', composition=parse_formula('C13H16ClNO'))]
    relations, counts = find_mass_relations(features, entries, 'positive', 10.0, [], 0)
    assert ('f', 'X1') not in relations and counts['adduct_conflict_dropped'] == 1

def test_ppm_bands_choose_by_theoretical_mz():
    bands = validate_ppm_bands([{'min_mz': None, 'max_mz': 200, 'ppm': 10}, {'min_mz': 200, 'max_mz': None, 'ppm': 5}])
    assert tolerance_ppm(150.0, 3.0, bands) == 10 and tolerance_ppm(200.0, 3.0, bands) == 5
    with pytest.raises(ValueError):
        validate_ppm_bands([{'min_mz': None, 'max_mz': 300, 'ppm': 10}, {'min_mz': 200, 'max_mz': None, 'ppm': 5}])
```

（`_feature(fid, mz, adduct=None)` はテスト内で `Feature(fid, fid, mz, 1.0, None, None, None, None, None, False, adduct)` を作る。）

- [ ] **Step 2: 失敗を確かめる** — Run: `C:/Python314/python.exe -m pytest tests/test_xenobiotic_mass.py -q` / Expected: FAIL（モジュールが無い）
- [ ] **Step 3: 実装する** — 付加体の規則: 特徴の `adduct` を空白を除いて比べ、`ion_label` で作れるラベル（全付加体 × 損失なし/3 種、両極性）のどれかに一致し、かつ関係の `ion` と違えば落とす。どれにも一致しなければ残して `adduct_unrecognized` を数える。`validate_ppm_bands` は `ppm` > 0 の有限数、`min_mz < max_mz`、重なり無しを検証し、`ValueError` を出す。
- [ ] **Step 4: 通ることを確かめる** — 同じコマンドで PASS。
- [ ] **Step 5: コミット** — `git commit -m "feat: 外来の分子の付加体・質量差の表と質量の筋を足す"`

---

### Task 4: 登録（register の検証・索引の情報・保存）

**Files:**
- Create: `massbank_context/xenobiotic/registry.py`
- Create: `massbank_context/xenobiotic/msms.py`（このタスクでは `open_index` だけ）
- Modify: `tests/xenobiotic_helpers.py`（`make_index` を足す）
- Test: `tests/test_xenobiotic_registry.py`

**Interfaces:**
- Consumes: Task 1 の formula 関数、`EvidenceStore.read_view / put_snapshot / current`、`index.schema.DDL`、`index.writer.read_meta`、`experiment.bundle.INCHIKEY_PATTERN`
- Produces:
  - `class XenobioticError(ValueError)`: `__init__(self, message: str, code: str = 'invalid_xenobiotic_request')`、属性 `code`（`_invoke` が拾う）
  - `ROLES = ('drug', 'drug_metabolite', 'vehicle', 'anesthetic', 'co_medication', 'other')`、`RECEIVED_BY = ('test', 'reference', 'both')`、`ENTRY_KEYS = {'entry_id', 'role', 'name', 'inchikey', 'formula', 'parent_entry_id', 'received_by'}`
  - `@dataclass(frozen=True) class XenoEntry: entry_id: str; role: str; name: str; inchikey: str | None; formula: str; composition: dict; monoisotopic_mass: float; parent_entry_id: str | None; received_by: str`（`formula` は Hill 表記に正規化）
  - `validate_entries(entries: list[dict]) -> list[XenoEntry]`
  - `describe_entries(conn: sqlite3.Connection, entries: Sequence[XenoEntry]) -> list[dict]`（各 entry に `n_records`（LC-ESI-MS2）, `polarities`（例 `POSITIVE:41`）, `precursor_types`, `warnings`）
  - `register(store: EvidenceStore, analysis_id: str, entries: list[dict]) -> dict`（`registration_id`, `entries`（`as_tsv` の文字列）, `n_entries`）
  - `load_registration(store: EvidenceStore, analysis_id: str) -> tuple[list[XenoEntry], str]`（無ければ `XenobioticError(..., code='missing_state')`）
  - `msms.open_index() -> tuple[sqlite3.Connection, dict]`（索引ファイルが無い・`read_meta` が `None` → `ValueError('missing_index: ...')`、`'LC-ESI-MS2'` が `acquisition_classes` に無い → `ValueError('incompatible_index: ...')`）
  - `tests/xenobiotic_helpers.py`: `make_index(path: Path) -> Path`。`DDL` を流し、`index_meta`（`bin_width=0.0075`, `acquisition_classes=["LC-ESI-MS2"]`, `n_records`, `n_inchikey`, `built_at`）、ケタミン（`YQEZLKZALYSWHR-UHFFFAOYSA-N`, `C13H16ClNO`, exact_mass 237.092）の POSITIVE `[M+H]+` レコード 2 件（precursor 238.0993）、ピーク（product ion: 125.0153 `C7H6Cl+` 強度 999, 207.0571 `C12H12ClO+` 500, 220.0888 `C13H15ClN+` 400, 89.0386 `C7H5+` 100, 238.0993 300。neutral loss は precursor − mz ≥ 0.5 のものを `to_bin` で）、おとりの化合物（`ZZZZZZZZZZZZZZ-ZZZZZZZZZZ-Z`, `C10H10O`, precursor 147.0804）1 件、`axis_background`（`LC-ESI-MS2`, `POSITIVE`）の idf: product ion 125.0153 → 4.0, 207.0571 → 1.0, 220.0888 → 1.0, 89.0386 → 1.0；neutral loss 113.0840 → 1.0, 31.0422 → 0.5, 18.0105 → 0.5, 149.0607 → 1.0（ビンは `to_bin(mz, 0.0075)`）。

登録の警告（`warnings` は `;` 区切り）: `no_inchikey`（InChIKey なし）、`not_in_index`（`compounds` にも `records` にも無い）、`no_lc_esi_ms2_records`、`formula_mismatch_index:<索引の Hill 表記>`、`mass_mismatch_index:<差 Da>`（`compounds.exact_mass` と計算値の差が 0.002 Da を超える）。

- [ ] **Step 1: 失敗するテストを書く**

```python
KET = dict(entry_id='X1', role='drug', name='ketamine', inchikey='YQEZLKZALYSWHR-UHFFFAOYSA-N',
           formula='C13ClH16NO', received_by='test')

def test_register_reports_index_and_normalizes_formula(analysis):     # analysis: v3 バンドル登録済みの (store, aid)、config.INDEX_PATH を make_index に差し替える fixture
    store, aid = analysis
    result = register(store, aid, [KET, dict(entry_id='X2', role='drug_metabolite', name='HNK',
                                             formula='C12H14ClNO2', parent_entry_id='X1', received_by='test')])
    rows = _tsv_rows(result['entries'])
    assert rows['X1']['formula'] == 'C13H16ClNO' and rows['X1']['n_records'] == '2'
    assert rows['X2']['warnings'] == 'no_inchikey'
    entries, registration_id = load_registration(store, aid)
    assert [e.entry_id for e in entries] == ['X1', 'X2'] and registration_id == result['registration_id']

@pytest.mark.parametrize('change', [
    {'entry_id': ''}, {'role': 'poison'}, {'received_by': 'all'}, {'formula': 'C13Xx'},
    {'inchikey': 'not-a-key'}, {'parent_entry_id': 'X9'}, {'parent_entry_id': 'X1'}, {'extra': 1}])
def test_invalid_entries_raise(analysis, change):
    store, aid = analysis
    with pytest.raises(XenobioticError):
        register(store, aid, [dict(KET, **change)])

def test_duplicate_entry_id_and_empty_list_raise(analysis):
    store, aid = analysis
    for entries in ([KET, KET], []):
        with pytest.raises(XenobioticError):
            register(store, aid, entries)

def test_mass_mismatch_warning(analysis):
    store, aid = analysis
    rows = _tsv_rows(register(store, aid, [dict(KET, formula='C13H18ClNO')])['entries'])
    assert 'formula_mismatch_index:C13H16ClNO' in rows['X1']['warnings']
    assert 'mass_mismatch_index' in rows['X1']['warnings']

def test_missing_registration_is_missing_state(analysis):
    store, aid = analysis
    with pytest.raises(XenobioticError) as info:
        load_registration(store, aid)
    assert info.value.code == 'missing_state'

def test_missing_index_raises(analysis, monkeypatch, tmp_path):
    monkeypatch.setattr(config, 'INDEX_PATH', tmp_path / 'none.sqlite')
    store, aid = analysis
    with pytest.raises(ValueError, match='^missing_index'):
        register(store, aid, [KET])
```

（`analysis` fixture と `_tsv_rows` はこのテストファイルに置く。Task 5〜8 のテストも同じ fixture を使うので、`tests/xenobiotic_helpers.py` に `registered_analysis(tmp_path, monkeypatch, *, version=3, detection='known') -> tuple[EvidenceStore, str]` として置いてもよい。）

- [ ] **Step 2: 失敗を確かめる** — Run: `C:/Python314/python.exe -m pytest tests/test_xenobiotic_registry.py -q` / Expected: FAIL
- [ ] **Step 3: 実装する** — `register` は `store.read_view(analysis_id)` で解析の存在と基準の版を取り、`validate_entries` → `open_index` → `describe_entries` → `store.put_snapshot(analysis_id, 'xenobiotics', {'registration_id', 'entries': [asdict(e)]}, expected_revision=base_revision)`。`registration_id = stable_id('xenobiotics', [analysis_id, canonical_json(entries_as_dicts)])`。登録し直しは全体の置き換え。`parent_entry_id` は自分以外の登録済みの entry を指すこと（循環は自分自身の参照だけを拒めば足りる。親の親は辿らない）。
- [ ] **Step 4: 通ることを確かめる** — 同じコマンドで PASS。
- [ ] **Step 5: コミット** — `git commit -m "feat: 外来の分子を検証して解析に登録し、索引の情報と警告を返す"`

---

### Task 5: MS/MS の共有の筋と forward_rank

**Files:**
- Modify: `massbank_context/xenobiotic/msms.py`
- Test: `tests/test_xenobiotic_msms.py`

**Interfaces:**
- Consumes: `core.binning.to_bin / probe_bins / PRODUCT_ION / NEUTRAL_LOSS`、`index.background.lookup_idf`、`spectral.matcher.match_spectrum_result / DEFAULT_PRECURSOR_TOLERANCE`、`core.models.QuerySpectrum`、Task 4 の `open_index`
- Produces:
  - `MSMS_ACQUISITION_CLASS = 'LC-ESI-MS2'`、`PRECURSOR_MARGIN = 0.5`、`FORWARD_TOP_K = 20`
  - `@dataclass(frozen=True) class RefBin: idf: float; intensity: float; formula: str | None`
  - `@dataclass(frozen=True) class ReferenceBins: inchikey: str; n_records: int; bins: dict[tuple[str, int], RefBin]; top: frozenset[tuple[str, int]]`
  - `reference_bins(conn, inchikey: str, index_polarity: str, bin_width: float, top_n: int) -> ReferenceBins | None`（レコードが無ければ `None`）
  - `query_axis_bins(spectrum: QuerySpectrum, bin_width: float) -> list[tuple[str, int]]`
  - `@dataclass(frozen=True) class SharedResult: score: float; n_shared: int; formulas: tuple[str, ...]`
  - `shared_evidence(ref: ReferenceBins, spectra: Sequence[QuerySpectrum], bin_width: float) -> SharedResult | None`（全スペクトルがピーク 0 本なら `None`。複数なら得点が最大のスペクトルの結果）
  - `forward_rank(conn, spectra: Sequence[QuerySpectrum], inchikey: str) -> int | None`

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_reference_bins_exclude_precursor_and_rank_top(index_conn):          # index_conn: make_index を ro で開く fixture
    ref = reference_bins(index_conn, KETAMINE, 'POSITIVE', 0.0075, 10)
    assert ref.n_records == 2
    assert ('product_ion', to_bin(238.0993, 0.0075)) not in ref.bins
    assert ref.bins[('product_ion', to_bin(125.0153, 0.0075))].formula == 'C7H6Cl+'
    assert len(reference_bins(index_conn, KETAMINE, 'POSITIVE', 0.0075, 2).top) == 4   # 軸ごとに 2 本

def test_parent_and_metabolite_scores(index_conn):
    ref = reference_bins(index_conn, KETAMINE, 'POSITIVE', 0.0075, 10)
    parent = shared_evidence(ref, [_spectrum(238.0990, [125.0153, 207.0571, 220.0888, 238.0990])], 0.0075)
    hnk = shared_evidence(ref, [_spectrum(240.0786, [125.0153, 89.0386, 222.0680])], 0.0075)
    assert parent.score == pytest.approx(0.8) and parent.n_shared == 6
    assert hnk.score == pytest.approx(0.55) and hnk.n_shared == 3
    assert hnk.formulas == ('C7H5+', 'C7H6Cl+')

def test_one_reference_bin_is_used_once(index_conn):
    ref = reference_bins(index_conn, KETAMINE, 'POSITIVE', 0.0075, 10)
    twice = shared_evidence(ref, [_spectrum(None, [125.0150, 125.0156])], 0.0075)
    assert twice.n_shared == 1

def test_empty_and_precursorless_spectra(index_conn):
    ref = reference_bins(index_conn, KETAMINE, 'POSITIVE', 0.0075, 10)
    assert shared_evidence(ref, [_spectrum(238.0990, [])], 0.0075) is None
    no_precursor = shared_evidence(ref, [_spectrum(None, [125.0153])], 0.0075)
    assert no_precursor.n_shared == 1 and no_precursor.score == pytest.approx(0.4)

def test_unknown_inchikey_has_no_reference(index_conn):
    assert reference_bins(index_conn, 'AAAAAAAAAAAAAA-AAAAAAAAAA-A', 'POSITIVE', 0.0075, 10) is None

def test_forward_rank(index_conn):
    assert forward_rank(index_conn, [_spectrum(238.0990, [125.0153, 207.0571])], KETAMINE) == 1
    assert forward_rank(index_conn, [_spectrum(500.3000, [125.0153])], KETAMINE) is None
```

（`_spectrum(precursor, mzs)` は `QuerySpectrum('q', precursor, 'POSITIVE', tuple(Peak(mz, 100.0) for mz in mzs))`。得点の期待値は Task 4 の idf から: 上位の分母 10.0、親は 4+1+1+(1+0.5+0.5) = 8.0、HNK は 4+1+0.5 = 5.5。）

- [ ] **Step 2: 失敗を確かめる** — Run: `C:/Python314/python.exe -m pytest tests/test_xenobiotic_msms.py -q` / Expected: FAIL
- [ ] **Step 3: 実装する** — `reference_bins`: `records`（inchikey・`acquisition_class='LC-ESI-MS2'`・polarity）と `peaks` を結合し、product ion で `bin * bin_width >= precursor_mz - PRECURSOR_MARGIN` を除き、(axis, bin) ごとに強度の最大と最初の非空の formula を取る。idf は `lookup_idf(conn, 'LC-ESI-MS2', index_polarity, axis, bin)`（中心ビン 1 行。±1 に広げない — background の docstring）。`top` は軸ごとに強度降順・ビン昇順で `top_n` 本。`query_axis_bins`: product ion は前駆体 − 0.5 未満だけ、neutral loss は前駆体 − mz ≥ 0.5 のものだけ（前駆体が無ければ neutral loss なし）。共有: クエリのビンごとに `probe_bins` の (中心, 下, 上) の順で、同じ軸の未使用の参照ビンを探す。`forward_rank`: スペクトルごとに `match_spectrum_result(conn, s, 'LC-ESI-MS2', DEFAULT_PRECURSOR_TOLERANCE, FORWARD_TOP_K)` の候補の InChIKey を出現順に重複除去し、entry の順位を取る。ピーク 0 本のスペクトルは飛ばす。
- [ ] **Step 4: 通ることを確かめる** — 同じコマンドで PASS。
- [ ] **Step 5: コミット** — `git commit -m "feat: 参照スペクトルとのビンの共有を IDF で重み付けし、順向きの順位を出す"`

---

### Task 6: 出方の筋（群の偏り・親との連動）

**Files:**
- Create: `massbank_context/stats/correlation.py`（`stats/__init__.py` の `__all__` に足す）
- Create: `massbank_context/xenobiotic/pattern.py`
- Test: `tests/test_xenobiotic_pattern.py`

**Interfaces:**
- Consumes: Task 2 の `Sample` と `parse_intensities` の行列の形、Task 3 の `THRESHOLDS`
- Produces:
  - `stats.correlation.spearman_rho(x: Sequence[float], y: Sequence[float]) -> float | None`（同順位は平均順位。n < 2 か分散 0 なら `None`）
  - `@dataclass(frozen=True) class SampleSides: test: tuple[str, ...]; reference: tuple[str, ...]`
  - `sample_sides(samples: Sequence[Sample], test_side: str) -> SampleSides`（`test_side` は `'a'` / `'b'`。group 空の試料は入れない）
  - `detection_counts(row: dict[str, tuple[float | None, bool | None]], sides: SampleSides) -> tuple[int | None, int | None]`（その側の `detected` が全部 `None` なら `None`）
  - `group_bias(counts: tuple[int | None, int | None], sides: SampleSides, received_by: str, thresholds: dict) -> bool | None`（`both` か数が `None` なら `None`。`test`: 検出 test > n_test / 2 かつ reference ≤ `bias_reference_max`。`reference` はその鏡）
  - `parent_rho(child: dict, parent: dict, sides: SampleSides) -> tuple[float | None, int]`（test 側の試料のうち、両方の値があり、どちらの `detected` も `False` でない試料で ρ と n）
  - `parent_linked(rho: float | None, n: int, thresholds: dict) -> bool`（n ≥ `rho_min_n` かつ ρ ≥ `rho_min`）

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_spearman_with_ties():
    assert spearman_rho([1, 2, 3, 4, 5], [5, 6, 7, 8, 7]) == pytest.approx(0.820783, abs=1e-6)
    assert spearman_rho([1, 2, 3], [4, 4, 4]) is None

def test_bias_for_test_only_feature(matrix, samples):                  # make_v3_bundle の行列と試料を読む fixture
    sides = sample_sides(samples, 'b')
    assert sides.test == ('T1', 'T2', 'T3', 'T4', 'T5', 'T6') and 'QC1' not in sides.reference
    assert detection_counts(matrix['x_msms'], sides) == (4, 1)
    assert group_bias((4, 1), sides, 'test', THRESHOLDS) is True
    assert group_bias((6, 6), sides, 'test', THRESHOLDS) is False
    assert group_bias((0, 6), sides, 'reference', THRESHOLDS) is True
    assert group_bias((6, 0), sides, 'both', THRESHOLDS) is None

def test_unknown_detection_is_unavailable_not_zero(unknown_matrix, samples):   # detection='unknown' のバンドル
    sides = sample_sides(samples, 'b')
    assert detection_counts(unknown_matrix['k_parent'], sides) == (None, None)
    assert group_bias((None, None), sides, 'test', THRESHOLDS) is None

def test_parent_rho_skips_imputed(matrix, samples):
    sides = sample_sides(samples, 'b')
    rho, n = parent_rho(matrix['k_hnk'], matrix['k_parent'], sides)
    assert rho == pytest.approx(1.0) and n == 6
    assert parent_linked(rho, n, THRESHOLDS)
    assert not parent_linked(1.0, 4, THRESHOLDS)
```

- [ ] **Step 2: 失敗を確かめる** — Run: `C:/Python314/python.exe -m pytest tests/test_xenobiotic_pattern.py -q` / Expected: FAIL
- [ ] **Step 3: 実装する** — `spearman_rho` は平均順位にしてから Pearson。`stats/` は leaf のまま（core も import しない）。
- [ ] **Step 4: 通ることを確かめる** — 同じコマンドで PASS。
- [ ] **Step 5: コミット** — `git commit -m "feat: 試料ごとの検出で群の偏りと親との連動を見る"`

---

### Task 7: 走査（水準・行・TSV・要約・行の取得）

**Files:**
- Create: `massbank_context/xenobiotic/screen.py`
- Test: `tests/test_xenobiotic_screen.py`

**Interfaces:**
- Consumes: Task 2〜6 のすべて（`role_path`, `parse_samples`, `parse_intensities`, `EvidenceStore.read_source`, `find_mass_relations`, `validate_ppm_bands`, `load_registration`, `open_index`, `reference_bins`, `shared_evidence`, `forward_rank`, `sample_sides`, `detection_counts`, `group_bias`, `parent_rho`, `parent_linked`, `THRESHOLDS*`, 表の版）
- Produces:
  - `SCREEN_COLUMNS = ['feature_id', 'mz', 'rt', 'log2fc', 'q_value', 'significant', 'level', 'entry_id', 'role', 'relation', 'theoretical_mz', 'ppm', 'msms_shared_score', 'msms_shared_n', 'msms_shared_formulas', 'forward_rank', 'detected_test', 'detected_reference', 'rho_parent', 'notes']`（spec §6 の順そのまま）
  - `LEVELS = ('A', 'B', 'C', 'D')`
  - `@dataclass(frozen=True) class PairEvidence: mass_parent: bool; mass_any: bool; msms_hit: bool; forward_rank: int | None; bias: bool | None; parent_link: bool`
  - `assign_level(pair: PairEvidence, received_by: str) -> str | None`
  - `resolve_polarity(manifest: dict, polarity: str | None) -> str`（`XenobioticError` の code: v2 で引数なし → `polarity_required`、v3 と引数の食い違い → `polarity_mismatch`、値が不正 → `invalid_xenobiotic_request`）
  - `resolve_test_side(manifest: dict, test_group: str) -> str`（`'a'` / `'b'`。どちらでもなければ `unknown_test_group`）
  - `run_screen(store: EvidenceStore, analysis_id: str, *, test_group: str, ppm: float, ppm_bands: list[dict] | None = None, polarity: str | None = None, max_shifts: int = 2, output_path: str | None = None) -> dict`（要約）
  - `get_screen_rows(store: EvidenceStore, analysis_id: str, *, feature_ids: list[str] | None = None, level: str | None = None, cursor: str | None = None, limit: int = 50) -> dict`

`assign_level` の規則（spec §5 をこの形に固定する）:
- A: `mass_parent` かつ `forward_rank == 1`
- B: `mass_any` かつ `msms_hit`
- C: (`mass_any` か `msms_hit`) かつ 出方。出方 = `received_by == 'both'` なら `parent_link` だけ、それ以外は `bias is True or parent_link`
- D: `mass_any` か `msms_hit`
- それ以外 `None`

`run_screen` の手順:
1. `store.read_view` → bundle・基準の版。`resolve_polarity`・`resolve_test_side`。`ppm` は正の有限数、`max_shifts` は 0〜2 の int、`validate_ppm_bands`。
2. `load_registration`。`open_index`（`bin_width` は `meta['bin_width']`、索引の極性は `polarity.upper()`）。
3. 質量の筋: `find_mass_relations(bundle.features, entries, polarity, ppm, bands, max_shifts)`。
4. MS/MS: entry ごとに `reference_bins`（InChIKey なし → 理由 `no_inchikey`、`None` → `no_reference_spectra`）。特徴ごとのスペクトルは `bundle.spectrum_links` → `bundle.spectra[file + '::' + id]`。結ばれたスペクトルが無ければ `no_msms`、全部ピーク 0 本なら `empty_spectrum`。全特徴 × 参照のある entry で `shared_evidence`（`msms_hit = score >= THRESHOLDS['msms_shared_min']`）。`forward_rank` は質量の一致か `msms_hit` がある組だけで計算する（全組で回すと照合が特徴 × entry 回になる）。
5. 出方: `intensities` が無ければ理由 `no_intensities`。あれば `read_source` → `parse_samples` / `parse_intensities`、`sample_sides`。特徴ごとに `detection_counts`、組ごとに `group_bias`。
6. 一巡目で A を決め、二巡目で親との連動（「設計との食い違い・補足」の親の決め方）と残りの水準を決める。
7. 行: 組の質量の一致ごとに 1 行（無ければ `relation='msms_only'` の 1 行）。どの組にも水準が無い特徴は `entry_id` 以降が空の 1 行。並びは bundle の特徴の順 → 水準の強い順 → entry_id → 質量差の数 → |ppm|。`notes` は `;` 区切りで `no_mz` / `no_msms` / `empty_spectrum` / `msms_unavailable:<理由>` / `pattern_unavailable:<理由>` / `parent_feature:<id>` / `adduct_unrecognized`。数値は TSV に書く前に 6 桁で丸める（`ppm` と得点は 4 桁）。
8. TSV: `output_path` が無ければ `config.ARTIFACTS_PATH / 'xenobiotic' / f"{analysis_id.split(':', 1)[1][:16]}_{polarity}.tsv"`（`:` は Windows のファイル名に使えない）。親ディレクトリを作り、`csv.writer(delimiter='\t', lineterminator='\n')`、UTF-8（BOM なし）で書く。sha256 を取る。
9. スナップショット: `store.put_snapshot(analysis_id, 'xenobiotic_screen', payload, expected_revision=base_revision)`。payload は `params`, `versions`, `thresholds`, `registration_id`, `summary`, `columns`, `rows`（リストのリスト）, `output_path`, `output_sha256`。戻りの revision を `screen_revision` とする。
10. 要約（戻り値）: `screen_revision`, `registration_id`, `polarity`, `test_group`, `versions`（`adduct_table` / `shift_table` / `thresholds`）, `thresholds`, `params`, `n_features`, `levels`（特徴の水準ごとの数と `none`）, `significant_levels`, `per_entry`（entry_id → `{'A': n, 'B': n}`）, `unavailable`（`msms`: entry_id → 理由、`n_features_without_msms`、`pattern`: 理由か `null`、`detection_unavailable_reason`）, `counts`（`adduct_conflict_dropped` 等）, `n_rows`, `output_path`, `output_sha256`, `next_tools: ['get_xenobiotic_screen']`。行は載せない。

`get_screen_rows`: `store.current(analysis_id, 'xenobiotic_screen')` が空なら `missing_state`。`limit` は 1〜200。未知の `feature_ids` は `unknown_feature_ids` に列挙（止めない）。`level` は A〜D か `None`。カーソルは `base64.urlsafe_b64encode(canonical_json({'rev': screen_revision, 'offset': n}))`。デコード失敗は `ValueError('invalid_cursor')`、`rev` が今の走査と違えば `ValueError('stale_cursor: 走査がやり直されました')`。戻り値は `screen_revision`, `n_matched`, `rows`（`as_tsv(SCREEN_COLUMNS, ...)`）, `next_cursor`, `unknown_feature_ids`。

- [ ] **Step 1: 失敗するテストを書く**

```python
@pytest.mark.parametrize('pair, received_by, expected', [
    (PairEvidence(True, True, True, 1, True, False), 'test', 'A'),
    (PairEvidence(True, True, True, 2, None, False), 'test', 'B'),
    (PairEvidence(False, True, False, None, True, False), 'test', 'C'),
    (PairEvidence(False, False, True, None, True, False), 'test', 'C'),
    (PairEvidence(False, True, False, None, None, False), 'test', 'D'),
    (PairEvidence(False, True, False, None, True, False), 'both', 'D'),
    (PairEvidence(False, True, False, None, None, True), 'both', 'C'),
    (PairEvidence(False, False, False, None, True, True), 'test', None),
])
def test_assign_level(pair, received_by, expected):
    assert assign_level(pair, received_by) == expected

def test_screen_levels_on_fixture(screened):          # 登録 X1=ケタミン(test)、X2=HNK(parent X1, InChIKey なし) の後に run_screen(test_group='Drug', ppm=5)
    summary, rows = screened                           # rows: 書いた TSV を DictReader で読んだもの
    level = _feature_levels(rows)                      # 特徴ごとに最も強い水準
    assert level == {'k_parent': 'A', 'k_hnk': 'B', 'k_frag': 'C', 'x_msms': 'C',
                     'mass_only': 'D', 'endo': '', 'no_mz': ''}
    hnk_x2 = _row(rows, 'k_hnk', 'X2')
    assert hnk_x2['level'] == 'C' and float(hnk_x2['rho_parent']) == pytest.approx(1.0)
    assert 'msms_unavailable:no_inchikey' in hnk_x2['notes']
    assert _row(rows, 'x_msms', 'X1')['relation'] == 'msms_only'
    assert summary['levels'] == {'A': 1, 'B': 1, 'C': 2, 'D': 1, 'none': 2}
    assert summary['versions'] == {'adduct_table': 'xeno-adducts-1', 'shift_table': 'xeno-shifts-1',
                                   'thresholds': THRESHOLDS_VERSION}
    assert 'rows' not in summary and len({r['feature_id'] for r in rows}) == 7

def test_feature_without_mz_is_kept(screened):
    _, rows = screened
    assert 'no_mz' in _row(rows, 'no_mz', '')['notes']

def test_v2_without_polarity_requires_argument(analysis_v2):
    with pytest.raises(XenobioticError) as info:
        run_screen(*analysis_v2, test_group='Drug', ppm=5)
    assert info.value.code == 'polarity_required'

def test_v2_screen_marks_pattern_unavailable(analysis_v2):
    summary = run_screen(*analysis_v2, test_group='Drug', ppm=5, polarity='positive')
    assert summary['unavailable']['pattern'] == 'no_intensities'
    assert summary['levels']['C'] == 0

@pytest.mark.parametrize('kwargs, code', [
    (dict(test_group='Placebo'), 'unknown_test_group'),
    (dict(test_group='Drug', polarity='negative'), 'polarity_mismatch'),
    (dict(test_group='Drug', max_shifts=3), 'invalid_xenobiotic_request'),
    (dict(test_group='Drug', ppm=0), 'invalid_xenobiotic_request'),
])
def test_screen_argument_errors(registered, kwargs, code):
    with pytest.raises(XenobioticError) as info:
        run_screen(*registered, **dict(dict(ppm=5), **kwargs))
    assert info.value.code == code

def test_max_shifts_zero_drops_two_shift_relation(registered):
    run_screen(*registered, test_group='Drug', ppm=5, max_shifts=0)
    rows = _read_tsv(_default_output(registered))
    assert _row(rows, 'k_hnk', 'X1')['relation'] == 'msms_only'

def test_get_rows_pages_and_filters(screened_store):
    store, aid = screened_store
    first = get_screen_rows(store, aid, level='A', limit=1)
    assert first['n_matched'] == 1 and first['rows'].splitlines()[0].split('\t') == SCREEN_COLUMNS
    page = get_screen_rows(store, aid, limit=2)
    assert page['next_cursor'] and get_screen_rows(store, aid, cursor=page['next_cursor'], limit=2)['rows']
    assert get_screen_rows(store, aid, feature_ids=['nope'])['unknown_feature_ids'] == ['nope']

def test_cursor_is_stale_after_rescreen(screened_store):
    store, aid = screened_store
    cursor = get_screen_rows(store, aid, limit=1)['next_cursor']
    run_screen(store, aid, test_group='Drug', ppm=10)
    with pytest.raises(ValueError, match='^stale_cursor'):
        get_screen_rows(store, aid, cursor=cursor)
```

（fixture `registered` / `screened` / `screened_store` / `analysis_v2` はこのファイルで Task 4 の `registered_analysis` と `register` から作る。`_default_output` は既定の出力パスの規則を写したもの。）

- [ ] **Step 2: 失敗を確かめる** — Run: `C:/Python314/python.exe -m pytest tests/test_xenobiotic_screen.py -q` / Expected: FAIL
- [ ] **Step 3: 実装する** — 上の手順どおり。索引の接続は `finally` で閉じる。
- [ ] **Step 4: 通ることを確かめる** — Run: `C:/Python314/python.exe -m pytest tests -q` / Expected: 全件 PASS
- [ ] **Step 5: コミット** — `git commit -m "feat: 3 つの筋から水準を付け、全特徴の TSV と要約を書く"`

---

### Task 8: MCP 公開と登録・文書

**Files:**
- Create: `massbank_context/tools/xenobiotic_tools.py`
- Modify: `massbank_context/server.py`（`evidence_tools` と同じループで `xenobiotic_tools.__all__` を登録）
- Modify: `tests/test_server_registration.py`（`EXPECTED_TOOLS` に 3 本）
- Modify: `CLAUDE.md`（アーキテクチャの図に `xenobiotic/` の 1 行、現状の段落に走査の 1 文と設計書へのリンク）、`README.md`（ツール一覧があればそこに 3 本）
- Test: `tests/test_xenobiotic_tools.py`

**Interfaces:**
- Consumes: `tools.evidence_tools._invoke`（`EvidenceStore(config.ARTIFACTS_PATH)` を渡し、`ValueError` の `code` 属性と `missing_index` / `incompatible_index` / `stale_cursor` / `invalid_cursor` の接頭辞を機械可読のエラーに変える）、Task 4 の `register`、Task 7 の `run_screen` / `get_screen_rows`
- Produces（`__all__` の順）:
  - `register_xenobiotics(analysis_id: str, entries: list[dict]) -> str`
  - `screen_xenobiotic_features(analysis_id: str, test_group: str, ppm: float, ppm_bands: list[dict] | None = None, polarity: str | None = None, max_shifts: int = 2, output_path: str | None = None) -> str`
  - `get_xenobiotic_screen(analysis_id: str, feature_ids: list[str] | None = None, level: str | None = None, cursor: str | None = None, limit: int = 50) -> str`

docstring は LLM が読む説明そのもの。各 1〜2 文で: 登録は「entries の欄、formula は必須、名前は照合に使わない」、走査は「test_group は comparison の group_a/group_b の名前、v2 は polarity 必須、全行は TSV に書き戻り値は要約だけ、水準は判定ではなく材料」、取得は「最大 200 行、next_cursor を cursor へ、走査し直すと古いカーソルは stale_cursor」。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_tools_round_trip(tmp_path, monkeypatch):
    aid = _load_v3(tmp_path, monkeypatch)                 # ARTIFACTS_PATH・INDEX_PATH を差し替えて load_evidence_bundle を呼ぶ
    reg = json.loads(register_xenobiotics(aid, [KET]))
    assert reg['registration_id'].startswith('xenobiotics:')
    summary = json.loads(screen_xenobiotic_features(aid, 'Drug', 5.0))
    assert summary['levels']['A'] == 1 and Path(summary['output_path']).is_file()
    rows = json.loads(get_xenobiotic_screen(aid, level='A'))
    assert rows['rows'].count('\n') == 1

def test_errors_are_machine_readable(tmp_path, monkeypatch):
    aid = _load_v3(tmp_path, monkeypatch)
    error = json.loads(screen_xenobiotic_features(aid, 'Drug', 5.0))['error']
    assert error['code'] == 'missing_state'
    error = json.loads(screen_xenobiotic_features(aid, 'Nobody', 5.0))['error']
    assert error['code'] == 'unknown_test_group'          # 引数の検証は登録の読み出しより先（Task 7 の手順 1→2）
```

`tests/test_server_registration.py` の `EXPECTED_TOOLS` に `"register_xenobiotics"`, `"screen_xenobiotic_features"`, `"get_xenobiotic_screen"` を足す（起点で `validate_comparison_claims` が既に載っていればそれも残す）。

- [ ] **Step 2: 失敗を確かめる** — Run: `C:/Python314/python.exe -m pytest tests/test_xenobiotic_tools.py tests/test_server_registration.py -q` / Expected: FAIL
- [ ] **Step 3: 実装する** — 3 本とも `_invoke(lambda store: ...)`。`test_no_structured_output` が通ること。
- [ ] **Step 4: 通ることを確かめる** — Run: `C:/Python314/python.exe -m pytest tests -q` / Expected: 全件 PASS
- [ ] **Step 5: サーバが立つことを確かめる** — Run: `C:/Python314/python.exe -c "import asyncio; from massbank_context import server; print(len(asyncio.run(server.mcp.list_tools())))"` / Expected: 起点のツール数 + 3
- [ ] **Step 6: コミット** — `git commit -m "feat: 外来の分子の登録・走査・行の取得を MCP に公開する"`

---

### Task 9: 実データでの較正と合格の確認（ゲートつき）

**ゲート:** ms-data-parser（`C:\Users\yuu18\Metabolomix_with_LLM`）の `main` に v3 の書き出しが入っていること。確かめ方: `metabolomix/analysis/evidence_export.py` に `contract_version=3` と `evidence-v3` があり、`tests/fixtures/` に v3 の見本のバンドルがある。2026-11-04 の時点では未実装（`contract_version=2`）。**無ければこのタスクを始めず、Task 8 までで止めて報告する。**

**Files:**
- Create: `tests/fixtures/evidence_v3/`（ms-data-parser の見本の写し）
- Create: `tests/fixtures/xenobiotic/calibration_expectations.json`
- Create: `tests/test_xenobiotic_calibration.py`
- Modify: `massbank_context/xenobiotic/tables.py`（`THRESHOLDS` と `THRESHOLDS_VERSION`）
- Modify: `pytest.ini`（marker `realdata` を足し、`addopts = -m "not network and not realdata"`）
- Modify: `docs/superpowers/specs/2026-11-04-xenobiotic-screen-design.md`（末尾に「## 10. 較正の結果」）、`docs/contracts/evidence-bundle-v3.md`（見本の共有の節を実際のパスで更新）

**Interfaces:**
- Consumes: Task 8 のツール（MCP 経由でも関数直接でもよい）、ms-data-parser の `dataset_export_evidence_bundle`
- Produces: 決まった `THRESHOLDS` / `THRESHOLDS_VERSION = 'xeno-thresholds-<較正した日付 YYYY-MM-DD>'`

- [ ] **Step 1: 見本の契約テスト** — ms-data-parser の v3 見本を `tests/fixtures/evidence_v3/` に写し、`tests/test_evidence_bundle_v3.py` に `test_shared_fixture_from_ms_data_parser_loads`（`load_bundle` が通り、`parse_intensities` の行数 = 特徴数 × a/b の試料数）を足す。通らなければ契約文書との食い違いを直す（どちらが正かは契約文書。ms-data-parser の設計書 §5）。`intensities_stage` の位置が食い違えばここで合わせる。
- [ ] **Step 2: v3 のバンドルを作る** — 各 run の主の比較（割れの組をまとめる前のスポット単位）を、reanalysis-study の `runs/<accession>/decisions.jsonl` と `runs/<accession>/results/differential_*.tsv` の頭のメタ（group_a / group_b）から再現し、`dataset_export_evidence_bundle` で `C:\Users\yuu18\datasets\reanalysis\<accession>\xenobiotic_calibration\<pos|neg>\bundle\` に書く。対象: ST003841 pos、ST002901 pos・neg、MTBLS9365 pos・neg。Console の出力は `C:\Users\yuu18\datasets\reanalysis\<accession>\pipeline\<pos|neg>\`。
- [ ] **Step 3: 登録して走らせる** — `ppm=10`（帯が run の `results/mass_accuracy_<pol>.json` にあれば `ppm_bands` で渡す）。entries:
  - ST003841（test = メトホルミン群）: `X1` drug metformin `C4H11N5` `XZWYZXLIPXDOLR-UHFFFAOYSA-N` test。
  - ST002901（test = フルオキセチン群）: `X1` drug fluoxetine `C17H18F3NO` `RTHCYVBBDHJXIQ-UHFFFAOYSA-N` test；`X2` drug_metabolite norfluoxetine `C16H16F3NO` `WIQRCHMSJFFONW-UHFFFAOYSA-N` parent X1 test；`X3` anesthetic ketamine `C13H16ClNO` `YQEZLKZALYSWHR-UHFFFAOYSA-N` both；`X4` drug_metabolite norketamine `C12H14ClNO` `BEQZHFIKTBVCAU-UHFFFAOYSA-N` parent X3 both；`X5` drug_metabolite hydroxynorketamine `C12H14ClNO2`（InChIKey なし。索引に無い）parent X3 both；`X6` drug_metabolite dehydronorketamine `C12H12ClNO`（同）parent X3 both；`X7` anesthetic xylazine `C12H16N2S`（InChIKey なし。索引に無い）both。
  - MTBLS9365（おとり。test = MCD 群）: `X1` drug amlodipine `C20H25ClN2O5` `HTIQEAQVCYTUBX-UHFFFAOYSA-N` test；`X2` drug fluoxetine（上と同じ）test。
- [ ] **Step 4: 特徴 ID を固定する** — 各 run の `runs/<accession>/report.md` の表の m/z・RT・組成式から v3 のバンドルの特徴を引き直し、`calibration_expectations.json` に `{accession, polarity, bundle_dir, test_group, entries, expect: [{feature_id, label, min_level}]}` で書く。min_level: 本体（メトホルミン・フルオキセチン・ケタミン）の特徴は `A`、既知の代謝物（ノルフルオキセチン・ノルケタミン・ヒドロキシノルケタミン・デヒドロノルケタミン）は `B`、レポートが「薬剤由来と推定」とした特徴（ST003841 の 1-Methylbiguanide・N-Isopropylimidodicarbonimidic diamide・N-Cyano-N',N'-dimethylguanidine・6-Propyl-1,3,5-triazine-2,4-diamine の 4 件、ST002901 の Autophinib（陰）・2-Ethyl-2-(p-tolyl)malonamide（陽）・4-Bromo-2-{…}（陰）の 3 件）は `C`。おとりは `{max_A: 0, max_B: 0}` と C の件数の上限。検出されていない本体・代謝物（バンドルに該当の特徴が無い）は `expect` に入れず、理由を `skipped` に書く。
- [ ] **Step 5: 閾値を決める** — `msms_shared_min` と `msms_top_n` を動かし（出方の既定値は設計どおりから始める）、陽性の期待をすべて満たし、おとりで A・B が 0 になる値を取る。満たす値が無ければ**質量の表を広げず**、満たせなかった特徴と理由を記録して Step 7 へ（手順書には入れない。spec §8）。
- [ ] **Step 6: 実データのテストを書いて通す** — `tests/test_xenobiotic_calibration.py` に `@pytest.mark.realdata` の `test_calibration_expectations`（期待の JSON を読み、バンドルが無ければ `pytest.skip`、各 run を走らせて `min_level` / `max_A` / `max_B` を確かめる）と、マーク無しの `test_thresholds_match_design_record`（`THRESHOLDS` と `THRESHOLDS_VERSION` が決めた値と等しい）。Run: `C:/Python314/python.exe -m pytest tests -q -m realdata` / Expected: PASS。Run: `C:/Python314/python.exe -m pytest tests -q` / Expected: 全件 PASS（realdata は選ばれない）。
- [ ] **Step 7: 設計書に追記する** — `## 10. 較正の結果（YYYY-MM-DD）`: run ごとの水準の件数（全特徴・有意な特徴）、陽性の対照の各特徴の水準、おとりの C・D の件数（偽陽性の目安）、決めた閾値と版、合否、満たせなかった点。あわせて、この計画の「設計との食い違い・設計の補足」の各項を「実装で決めたこと」として転記する。
- [ ] **Step 8: コミット** — `git commit -m "test: 陽性とおとりの対照で外来の分子の走査の閾値を較正する"`

---

## Self-Review の記録

- spec §3.1 登録 → Task 4 / §3.2 v3 → Task 2 / §4.1 → Task 3 / §4.2 → Task 5 / §4.3 → Task 6 / §5 → Task 7 / §6 出力 → Task 7・8 / §7 置き場所・版・`no_msms` → Task 3・5・7 / §8 → Task 9（単体テストの 5 項目は Task 1・3・5・6・7・2）/ §9 は対象外。
- 型の一貫: `XenoEntry`（Task 4）を Task 3 は属性 `entry_id` / `composition` だけで使う。`role_path`・`parse_samples`・`parse_intensities`・`read_source`（Task 2）を Task 7 が使う。`THRESHOLDS` のキーは Task 3 で定義し Task 5〜7・9 が同じ名前で読む。
