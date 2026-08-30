# massbank-context Phase 1 実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** MassBank-RDF の TTL からローカル SQLite 索引を構築する CLI と、実験スペクトルを照合して候補化合物を返す MCP ツール 4 本を作る。

**Architecture:** TTL（世代A/B の両対応）を rdflib で読み、`acquisition_class` と `resolution_class` で層別しながら SQLite に格納する。ビン幅は同一 InChIKey 内の m/z 散らばりを実測して決める。照合は IDF 重み付きカバレッジで行い、バイアス層は Phase 3 まで載せない。

**Tech Stack:** Python 3.14.4 (`C:/Python314/python.exe`) / rdflib（新規導入）/ 標準ライブラリの sqlite3 / FastMCP / pytest

**Spec:** `C:/Users/yuu18/Lipidmix_with_LLM/docs/superpowers/specs/2026-08-28-massbank-context-mcp-design.md`

**実装先リポジトリ:** `C:/Users/yuu18/massbank-context`（新規。Lipidmix とは別リポジトリ）

## Global Constraints

- Python は `C:/Python314/python.exe` を使う。`.venv` は作らない（Lipidmix と同じ流儀）
- テストは常にリポジトリルートから `C:/Python314/python.exe -m pytest tests -q`
- **全 MCP ツールに `structured_output=False` を付ける**。付けないと FastMCP が outputSchema を導出し、同じ内容が content と structuredContent の両方で送られる（＝2倍）
- **JSON は `json_payload()` で返す。`json.dumps(..., indent=2)` を書かない**
- **float は返す前に丸める**（既定 6 桁、m/z は 4 桁）
- **行が並ぶ一覧は TSV（列名 1 回）で返す**
- `massbank_context/core/` は依存グラフの leaf。上位モジュールを import してはならない
- **スコアを合成してはならない。並び順を決めるのは `spectral_evidence` だけ**
- `band` の定義: プロダクトイオン一致 2 本以上 = `PASS` / 1 本 = `WEAK` / 0 本（precursor のみ一致）= `PRECURSOR_ONLY`
- 索引の既定は `resolution_class == "high"` のみ
- **`acquisition_class` を跨ぐ照合は行わない**。跨いだ要求は `incompatible_index` で拒否する

---

## File Structure

```
massbank-context/
  massbank_context/
    __init__.py
    core/
      __init__.py
      models.py          SpectrumRecord / Peak / QuerySpectrum のデータクラス（leaf）
      classify.py        acquisition_class / resolution_class の導出（leaf）
      serialization.py   json_payload / round_floats（leaf・Lipidmix から移植）
      errors.py          missing_index / missing_state / incompatible_index 封筒（leaf）
      binning.py         m/z ↔ ビン番号の変換（leaf）
      session.py         セッション状態（読み込んだクエリスペクトル）
    ttl/
      __init__.py
      generation.py      世代A/B の判定
      reader_gen_a.py    世代A TTL → SpectrumRecord
      reader_gen_b.py    世代B TTL → SpectrumRecord
      reader.py          世代を判定して適切なリーダへ振り分ける入口
    index/
      __init__.py
      schema.py          SQLite の DDL
      writer.py          SpectrumRecord → SQLite
      mass_error.py      質量誤差プロファイルの計測とビン幅の決定
      background.py      axis_background（IDF）の算出
      __main__.py        索引ビルド CLI
    query/
      __init__.py
      msp.py             .msp / .mgf リーダ
      dcl.py             MS-DIAL .dcl リーダ（Lipidmix から移植）
      loader.py          拡張子で振り分ける入口
    spectral/
      __init__.py
      matcher.py         照合エンジン（IDF 重み付きカバレッジ・band 判定）
    tools/
      __init__.py
      index_tools.py     index_status
      spectra_tools.py   load_spectra / inspect_spectrum
      match_tools.py     match_spectra
    server.py            薄いファサード。import 副作用でツールを登録する
  tests/
    fixtures/ttl/gen_a/  実 TTL を縮小コピーしたもの
    fixtures/ttl/gen_b/  SPARQL 例から再構成した合成 TTL
    fixtures/query/      .msp / .mgf のサンプル
    test_classify.py
    test_ttl_generation.py
    test_reader_gen_a.py
    test_reader_gen_b.py
    test_binning.py
    test_index_writer.py
    test_mass_error.py
    test_background.py
    test_query_loader.py
    test_matcher.py
    test_tools.py
  requirements.txt
  pytest.ini
  README.md
```

---

### Task 1: リポジトリ骨格と層の導出ロジック

`acquisition_class` / `resolution_class` の導出は spec §6.1 / §6.2 の中心であり、
かつ純粋関数で I/O を持たない。ここを最初に固める。

**Files:**
- Create: `C:/Users/yuu18/massbank-context/requirements.txt`
- Create: `C:/Users/yuu18/massbank-context/pytest.ini`
- Create: `C:/Users/yuu18/massbank-context/massbank_context/__init__.py`
- Create: `C:/Users/yuu18/massbank-context/massbank_context/core/__init__.py`
- Create: `C:/Users/yuu18/massbank-context/massbank_context/core/classify.py`
- Test: `C:/Users/yuu18/massbank-context/tests/test_classify.py`

**Interfaces:**
- Consumes: なし
- Produces:
  - `derive_acquisition_class(instrument_type: str | None, ms_type: str | None) -> str`
    返り値は `"LC-ESI-MS2"` / `"GC-EI-MS"` / `"MALDI-MS2"` / `"other"` のいずれか
  - `derive_resolution_class(instrument_type: str | None, resolution_value: float | None = None) -> str`
    返り値は `"high"` / `"low"` / `"unknown"` のいずれか
  - `NEUTRAL_LOSS_CAPABLE: frozenset[str]` — 中性損失軸を使える acquisition_class の集合

- [ ] **Step 1: リポジトリを作り、依存を書く**

```bash
mkdir -p C:/Users/yuu18/massbank-context/massbank_context/core
mkdir -p C:/Users/yuu18/massbank-context/tests
cd C:/Users/yuu18/massbank-context
git init
```

`requirements.txt`:

```
rdflib>=7.0
mcp>=1.0
```

`pytest.ini`:

```
[pytest]
testpaths = tests
```

`massbank_context/__init__.py` と `massbank_context/core/__init__.py` は空ファイル。

- [ ] **Step 2: 失敗するテストを書く**

テストケースは spec §3.5 で実測した `instrument_type` 文字列をそのまま使う。

`tests/test_classify.py`:

```python
import pytest

from massbank_context.core.classify import (
    NEUTRAL_LOSS_CAPABLE,
    derive_acquisition_class,
    derive_resolution_class,
)


@pytest.mark.parametrize("instrument_type,ms_type,expected", [
    # spec §3.5 の実測値
    ("EI-B", "MS", "GC-EI-MS"),
    ("LC-ESI-ITFT", "MS2", "LC-ESI-MS2"),
    ("LC-ESI-QQ", "MS2", "LC-ESI-MS2"),
    ("LC-ESI-QTOF", "MS2", "LC-ESI-MS2"),
    ("LC-ESI-QIT", "MS2", "LC-ESI-MS2"),
    ("ESI-QTOF", "MS2", "LC-ESI-MS2"),      # LC 接頭辞が無い版
    ("LC-ESI-QFT", "MS2", "LC-ESI-MS2"),
    ("MALDI-TOFTOF", "MS2", "MALDI-MS2"),
    # 判定不能
    ("LC-ESI-QTOF", "MS3", "other"),
    ("EI-B", "MS2", "other"),               # EI の MS2 は想定外
    (None, "MS2", "other"),
    ("LC-ESI-QTOF", None, "other"),
])
def test_derive_acquisition_class(instrument_type, ms_type, expected):
    assert derive_acquisition_class(instrument_type, ms_type) == expected


def test_esi_is_not_mistaken_for_ei():
    """'ESI' の部分文字列一致で 'EI' を拾ってはならない。"""
    assert derive_acquisition_class("LC-ESI-ITFT", "MS2") == "LC-ESI-MS2"


@pytest.mark.parametrize("instrument_type,expected", [
    ("LC-ESI-ITFT", "high"),     # Orbitrap 系
    ("LC-ESI-QFT", "high"),
    ("LC-ESI-QTOF", "high"),
    ("MALDI-TOFTOF", "high"),
    ("ESI-QTOF", "high"),
    ("LC-ESI-QQ", "low"),        # 三連四重極
    ("LC-ESI-QIT", "low"),       # イオントラップ
    ("EI-B", "unknown"),         # 磁場型。MassBank の EI は公称質量が多い
    (None, "unknown"),
])
def test_derive_resolution_class(instrument_type, expected):
    assert derive_resolution_class(instrument_type) == expected


def test_resolution_value_does_not_override_class():
    """mbo:resolution は補強にのみ使い、instrument_type の判定を覆さない。

    spec §6.2: resolution は約半数で欠損するため主根拠にできない。
    """
    assert derive_resolution_class("LC-ESI-QTOF", 7000.0) == "high"
    assert derive_resolution_class("LC-ESI-QQ", 50000.0) == "low"


def test_neutral_loss_capable_excludes_gc_ei():
    """GC-EI-MS は precursor を持たないため中性損失軸を使えない（spec §6.1）。"""
    assert "LC-ESI-MS2" in NEUTRAL_LOSS_CAPABLE
    assert "GC-EI-MS" not in NEUTRAL_LOSS_CAPABLE
```

- [ ] **Step 3: テストが落ちることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_classify.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'massbank_context.core.classify'`

- [ ] **Step 4: 実装する**

`massbank_context/core/classify.py`:

```python
"""instrument_type / ms_type から索引の層を導出する。

spec §6.1 / §6.2。この 2 つの層は索引の分割単位であり、背景分布と濃縮度統計の
計算範囲を決める。測定原理の異なるスペクトルを同じ背景で数えることは無意味なので、
判定不能なものは "other" / "unknown" に落として明示的に扱う。

依存グラフの leaf（stdlib のみ）。
"""
from __future__ import annotations

#: 中性損失軸（precursor - fragment）を使える acquisition_class。
#: GC-EI-MS は precursor という概念自体を持たないため除外する。
NEUTRAL_LOSS_CAPABLE = frozenset({"LC-ESI-MS2", "MALDI-MS2"})

_ACQUISITION_TABLE = {
    ("LC-ESI", "MS2"): "LC-ESI-MS2",
    ("GC-EI", "MS"): "GC-EI-MS",
    ("MALDI", "MS2"): "MALDI-MS2",
}

_HIGH_ANALYZER_SUFFIXES = ("FT", "TOF")
_LOW_ANALYZERS = frozenset({"QQ", "QIT", "IT", "Q"})


def _tokens(instrument_type: str | None) -> list[str]:
    """'LC-ESI-ITFT' -> ['LC', 'ESI', 'ITFT']。

    部分文字列一致だと 'ESI' から 'EI' を誤って拾う経路が生まれるため、
    必ずトークンに割ってから判定する。
    """
    if not instrument_type:
        return []
    return [t for t in instrument_type.upper().replace("_", "-").split("-") if t]


def _ionization_family(tokens: list[str]) -> str | None:
    token_set = set(tokens)
    if "MALDI" in token_set:
        return "MALDI"
    if "ESI" in token_set:
        return "LC-ESI"
    if "EI" in token_set:
        return "GC-EI"
    return None


def derive_acquisition_class(instrument_type: str | None,
                             ms_type: str | None) -> str:
    """索引の分割単位を返す（spec §6.1）。

    判定できないものは "other" を返す。呼び出し側は "other" を索引対象外として
    扱ってよいが、黙って別の層に混ぜてはならない。
    """
    family = _ionization_family(_tokens(instrument_type))
    if family is None or not ms_type:
        return "other"
    return _ACQUISITION_TABLE.get((family, ms_type.strip().upper()), "other")


def derive_resolution_class(instrument_type: str | None,
                            resolution_value: float | None = None) -> str:
    """質量分解能の層を返す（spec §6.2）。

    mbo:resolution は実測で約半数が欠損していたため、判定は instrument_type の
    アナライザ部を主根拠とする。resolution_value は記録・補強のためだけに受け取り、
    判定を覆さない。
    """
    tokens = _tokens(instrument_type)
    if not tokens:
        return "unknown"
    analyzer = tokens[-1]
    if analyzer in _LOW_ANALYZERS:
        return "low"
    if any(analyzer.endswith(suffix) for suffix in _HIGH_ANALYZER_SUFFIXES):
        return "high"
    return "unknown"
```

- [ ] **Step 5: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_classify.py -q`
Expected: PASS（16 件）

- [ ] **Step 6: コミット**

```bash
cd C:/Users/yuu18/massbank-context && git add -A && git commit -m "feat: acquisition_class / resolution_class の導出とリポジトリ骨格"
```

---

### Task 2: データモデルと共通ヘルパ

**Files:**
- Create: `massbank_context/core/models.py`
- Create: `massbank_context/core/serialization.py`
- Create: `massbank_context/core/errors.py`
- Test: `tests/test_models.py`

**Interfaces:**
- Consumes: Task 1 の `derive_acquisition_class` / `derive_resolution_class`
- Produces:
  - `Peak(mz: float, relative_intensity: float, formula: str | None, neutral_loss: float | None)`
  - `SpectrumRecord` — 下記フィールドと `acquisition_class` / `resolution_class` プロパティ
  - `QuerySpectrum(spectrum_id: str, precursor_mz: float | None, polarity: str | None, peaks: tuple[Peak, ...], source: str)`
  - `json_payload(obj) -> str` / `round_floats(obj, digits=6)`
  - `missing_index(required_action: str) -> str`
  - `missing_state(state: str, required_tools: list[str], message: str) -> str`
  - `incompatible_index(query_class: str, available: list[str], message: str) -> str`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_models.py`:

```python
import json

from massbank_context.core.errors import (
    incompatible_index,
    missing_index,
    missing_state,
)
from massbank_context.core.models import Peak, QuerySpectrum, SpectrumRecord
from massbank_context.core.serialization import json_payload, round_floats


def _record(**overrides) -> SpectrumRecord:
    base = dict(
        accession="MSBNK-TEST-0001",
        inchikey="XSDQTOBWRPYKKA-UHFFFAOYSA-N",
        precursor_mz=230.0552,
        precursor_type="[M+H]+",
        polarity="POSITIVE",
        instrument_type="LC-ESI-ITFT",
        ms_type="MS2",
        resolution_value=15000.0,
        collision_energy="35% (nominal)",
        splash="splash10-001i-0290000000-657972e0af96e9da0162",
        license="CC0",
        molecular_formula="C6H8ClN7O",
        exact_mass=229.0479,
        compound_names=("Amiloride",),
        ch_compound_class="N/A; Environmental Standard",
        peaks=(Peak(mz=86.0349, relative_intensity=11.0,
                    formula="C2H4N3O+", neutral_loss=None),),
    )
    base.update(overrides)
    return SpectrumRecord(**base)


def test_record_derives_layers():
    rec = _record()
    assert rec.acquisition_class == "LC-ESI-MS2"
    assert rec.resolution_class == "high"


def test_record_without_precursor_is_gc_ei():
    rec = _record(instrument_type="EI-B", ms_type="MS", precursor_mz=None,
                  precursor_type=None)
    assert rec.acquisition_class == "GC-EI-MS"


def test_json_payload_has_no_whitespace():
    """戻り値の空白は LLM の文脈をそのまま食う（spec §12）。"""
    out = json_payload({"a": 1, "b": [1, 2]})
    assert out == '{"a":1,"b":[1,2]}'
    assert " " not in out


def test_round_floats_preserves_bool():
    assert round_floats({"x": 1.23456789, "flag": True}) == {
        "x": 1.234568, "flag": True}


def test_missing_index_envelope():
    payload = json.loads(missing_index("CLI: python -m massbank_context.index build"))
    assert payload["error"]["code"] == "missing_index"
    assert "required_action" in payload["error"]


def test_missing_state_envelope():
    payload = json.loads(missing_state("spectra", ["load_spectra"], "説明"))
    assert payload["error"]["code"] == "missing_state"
    assert payload["error"]["required_tools"] == ["load_spectra"]


def test_incompatible_index_envelope():
    payload = json.loads(
        incompatible_index("LC-ESI-MS2", ["GC-EI-MS"], "説明"))
    assert payload["error"]["code"] == "incompatible_index"
    assert payload["error"]["query_acquisition_class"] == "LC-ESI-MS2"
    assert payload["error"]["available"] == ["GC-EI-MS"]


def test_query_spectrum_holds_source():
    q = QuerySpectrum(spectrum_id="q1", precursor_mz=230.0552,
                      polarity="POSITIVE",
                      peaks=(Peak(86.0349, 11.0, None, None),),
                      source="sample.msp")
    assert q.source == "sample.msp"
    assert len(q.peaks) == 1
```

- [ ] **Step 2: テストが落ちることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_models.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'massbank_context.core.models'`

- [ ] **Step 3: 実装する**

`massbank_context/core/serialization.py`（Lipidmix `lipidmix/core/serialization.py` と同じ挙動）:

```python
"""LLM へ返す文字列の直列化。依存グラフの leaf（stdlib のみ）。

戻り値は LLM の文脈をそのまま占める。json.dumps(..., indent=2) は 1 ネストにつき
改行と空白を全キーに掛けるため、Lipidmix の実測では戻り値の 15〜57% が空白だった。
人間向けの整形はクライアント側の責務なので、サーバは最小形で返す。
"""
from __future__ import annotations

import json

DEFAULT_FLOAT_DIGITS = 6
MZ_DIGITS = 4


def json_payload(obj) -> str:
    """ツールの戻り値を最小トークンの JSON 文字列にする。"""
    return json.dumps(obj, ensure_ascii=False, separators=(",", ":"))


def round_floats(obj, digits: int = DEFAULT_FLOAT_DIGITS):
    """入れ子の dict/list を辿って float を丸める（非破壊・型は保つ）。

    bool は int の派生だが数値として丸めてはいけないので素通しする。
    """
    if isinstance(obj, bool):
        return obj
    if isinstance(obj, float):
        return round(obj, digits)
    if isinstance(obj, list):
        return [round_floats(item, digits) for item in obj]
    if isinstance(obj, tuple):
        return tuple(round_floats(item, digits) for item in obj)
    if isinstance(obj, dict):
        return {key: round_floats(value, digits) for key, value in obj.items()}
    return obj


def as_tsv(columns: list[str], rows: list[list]) -> str:
    """行が並ぶ一覧を TSV で返す（列名は 1 回だけ）。

    JSON の配列オブジェクトは行ごとにキー名を繰り返すため、行数に比例して無駄が増える。
    """
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join("" if v is None else str(v) for v in row))
    return "\n".join(lines)
```

`massbank_context/core/errors.py`:

```python
"""クライアント非依存のエラーエンベロープ（spec §8）。

MCP クライアントが「次に何をすればよいか」をサーバ固有の日本語文面を解釈せずに
機械的に判断できるようにする。本文の JSON に載せるのは、使用中の MCP SDK では
isError=true の経路が structuredContent を捨てるため。

依存グラフの leaf（serialization のみ）。
"""
from __future__ import annotations

from massbank_context.core.serialization import json_payload


def missing_index(required_action: str) -> str:
    """索引が無い。ツールでは直せないので CLI を案内する。"""
    if not required_action:
        raise ValueError("required_action は必須です。")
    return json_payload({
        "error": {
            "code": "missing_index",
            "required_action": required_action,
        }
    })


def missing_state(state: str, required_tools: list[str], message: str) -> str:
    """前提となるセッション状態が無い。ツールで復旧できる。"""
    if not state or not required_tools or not message:
        raise ValueError("state / required_tools / message はいずれも必須です。")
    return json_payload({
        "error": {
            "code": "missing_state",
            "state": state,
            "required_tools": list(required_tools),
            "message": message,
        }
    })


def incompatible_index(query_acquisition_class: str,
                       available: list[str],
                       message: str) -> str:
    """索引の層が要求と合わない（spec §6.1）。

    測定原理の異なるスペクトルを混ぜて照合することは無意味なので、黙って
    別の層で代用せず明示的に拒否する。
    """
    if not query_acquisition_class or not message:
        raise ValueError("query_acquisition_class / message は必須です。")
    return json_payload({
        "error": {
            "code": "incompatible_index",
            "query_acquisition_class": query_acquisition_class,
            "available": list(available),
            "message": message,
        }
    })
```

`massbank_context/core/models.py`:

```python
"""索引とクエリのデータモデル。依存グラフの leaf（classify のみ）。"""
from __future__ import annotations

from dataclasses import dataclass, field

from massbank_context.core.classify import (
    derive_acquisition_class,
    derive_resolution_class,
)


@dataclass(frozen=True)
class Peak:
    """1 本のピーク。

    neutral_loss は世代B TTL では mbo:neutral_loss から読み、世代A では
    precursor_mz - mz で後から算出する（spec §6.3）。算出できない場合は None。
    """
    mz: float
    relative_intensity: float
    formula: str | None = None
    neutral_loss: float | None = None


@dataclass(frozen=True)
class SpectrumRecord:
    """MassBank レコード 1 件（spec §6.7 の records + peaks）。"""
    accession: str
    inchikey: str | None
    precursor_mz: float | None
    precursor_type: str | None
    polarity: str | None
    instrument_type: str | None
    ms_type: str | None
    resolution_value: float | None
    collision_energy: str | None
    splash: str | None
    license: str | None
    molecular_formula: str | None
    exact_mass: float | None
    compound_names: tuple[str, ...] = ()
    ch_compound_class: str | None = None
    peaks: tuple[Peak, ...] = ()

    @property
    def acquisition_class(self) -> str:
        return derive_acquisition_class(self.instrument_type, self.ms_type)

    @property
    def resolution_class(self) -> str:
        return derive_resolution_class(self.instrument_type,
                                       self.resolution_value)


@dataclass(frozen=True)
class QuerySpectrum:
    """ユーザーが持ち込んだ実験スペクトル 1 本。"""
    spectrum_id: str
    precursor_mz: float | None
    polarity: str | None
    peaks: tuple[Peak, ...] = ()
    source: str = ""
    retention_time: float | None = None
```

- [ ] **Step 4: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_models.py -q`
Expected: PASS（8 件）

- [ ] **Step 5: コミット**

```bash
cd C:/Users/yuu18/massbank-context && git add -A && git commit -m "feat: データモデル・直列化ヘルパ・エラー封筒"
```

---

### Task 3: 世代判定と世代A TTL リーダ

**Files:**
- Create: `massbank_context/ttl/__init__.py`
- Create: `massbank_context/ttl/generation.py`
- Create: `massbank_context/ttl/reader_gen_a.py`
- Create: `tests/fixtures/ttl/gen_a/MSBNK-TEST-A0001.ttl`
- Test: `tests/test_reader_gen_a.py`

**Interfaces:**
- Consumes: `SpectrumRecord` / `Peak`（Task 2）
- Produces:
  - `detect_generation(graph: rdflib.Graph) -> str` — `"gen_a"` / `"gen_b"` / `"unknown"`
  - `parse_gen_a(graph: rdflib.Graph) -> list[SpectrumRecord]`
  - `MBO: rdflib.Namespace` — `http://www.massbank.jp/ontology/`

- [ ] **Step 1: 依存を入れる**

```bash
C:/Python314/python.exe -m pip install "rdflib>=7.0"
```

- [ ] **Step 2: fixture を作る**

`tests/fixtures/ttl/gen_a/MSBNK-TEST-A0001.ttl` — 実 TTL（`MSBNK-HBM4EU-HB000104`）を
ピーク 3 本に縮小したもの。**構造は実物のまま変えない**。

```turtle
@prefix dcterms: <http://purl.org/dc/terms/> .
@prefix inchikey: <http://identifiers.org/inchikey/> .
@prefix massbank: <https://identifiers.org/massbank/> .
@prefix mbo: <http://www.massbank.jp/ontology/> .
@prefix obo: <http://purl.obolibrary.org/obo/> .
@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
@prefix schema: <https://schema.org/> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .

massbank:MSBNK-TEST-A0001 a mbo:MassSpectrum ;
    dcterms:identifier "MSBNK-TEST-A0001" ;
    dcterms:license <https://creativecommons.org/publicdomain/zero/1.0/> ;
    mbo:analytical_methods_and_conditions [ a mbo:AnalyticalMethods ;
            mbo:ac_instrument "LTQ Orbitrap XL Thermo Scientific" ;
            mbo:collision_energy "35% (nominal)" ;
            mbo:fragmentation_mode "CID" ;
            mbo:instrument_type "LC-ESI-ITFT" ;
            mbo:ion_mode "POSITIVE" ;
            mbo:ms_type "MS2" ;
            mbo:resolution "15000" ] ;
    mbo:compound [ a schema:chemicalsubstance ;
            rdfs:label "Amiloride" ;
            mbo:ch_compound_class "N/A; Environmental Standard" ;
            mbo:ch_exact_mass 2.290479e+02 ;
            schema:inChIKey inchikey:XSDQTOBWRPYKKA-UHFFFAOYSA-N ;
            schema:molecularFormula "C6H8ClN7O" ;
            schema:name "Amiloride" ] ;
    mbo:has_peak [ mbo:intensity 5.92827e+04 ;
            mbo:mz 8.60349e+01 ;
            mbo:peak_sequential_number 1 ;
            mbo:relative_intensity 11 ],
        [ mbo:intensity 1.271685e+06 ;
            mbo:mz 1.710069e+02 ;
            mbo:peak_sequential_number 2 ;
            mbo:relative_intensity 226 ],
        [ mbo:intensity 5.61971e+06 ;
            mbo:mz 2.300552e+02 ;
            mbo:peak_sequential_number 3 ;
            mbo:relative_intensity 999 ] ;
    mbo:has_peak_annotations [ mbo:mz 8.60349e+01 ;
            mbo:tentative_formula "C2H4N3O+" ],
        [ mbo:mz 1.710069e+02 ;
            mbo:tentative_formula "C5H4ClN4O+" ] ;
    mbo:mass_spectram_meta_data [ a mbo:MassSpectramMetaData ;
            mbo:base_peak 2.300551e+02 ;
            mbo:precursor_mz_value "230.0552"^^xsd:float ;
            mbo:precursor_type_value "[M+H]+" ] ;
    mbo:pk_num_peak 3 ;
    mbo:pk_splash "splash10-001i-0290000000-657972e0af96e9da0162" .
```

- [ ] **Step 3: 失敗するテストを書く**

`tests/test_reader_gen_a.py`:

```python
from pathlib import Path

import rdflib

from massbank_context.ttl.generation import detect_generation
from massbank_context.ttl.reader_gen_a import parse_gen_a

FIXTURE = Path(__file__).parent / "fixtures" / "ttl" / "gen_a" / "MSBNK-TEST-A0001.ttl"


def _graph() -> rdflib.Graph:
    g = rdflib.Graph()
    g.parse(FIXTURE, format="turtle")
    return g


def test_detects_generation_a():
    assert detect_generation(_graph()) == "gen_a"


def test_parses_single_record():
    records = parse_gen_a(_graph())
    assert len(records) == 1


def test_record_metadata():
    rec = parse_gen_a(_graph())[0]
    assert rec.accession == "MSBNK-TEST-A0001"
    assert rec.inchikey == "XSDQTOBWRPYKKA-UHFFFAOYSA-N"
    assert rec.precursor_mz == 230.0552
    assert rec.precursor_type == "[M+H]+"
    assert rec.polarity == "POSITIVE"
    assert rec.instrument_type == "LC-ESI-ITFT"
    assert rec.ms_type == "MS2"
    assert rec.resolution_value == 15000.0
    assert rec.molecular_formula == "C6H8ClN7O"
    assert rec.splash.startswith("splash10-")


def test_peaks_are_sorted_by_mz():
    rec = parse_gen_a(_graph())[0]
    mzs = [p.mz for p in rec.peaks]
    assert mzs == sorted(mzs)
    assert len(mzs) == 3


def test_tentative_formula_is_joined_by_mz():
    """has_peak と has_peak_annotations は別ブロックで、m/z で突き合わせる。"""
    rec = parse_gen_a(_graph())[0]
    by_mz = {round(p.mz, 4): p.formula for p in rec.peaks}
    assert by_mz[86.0349] == "C2H4N3O+"
    assert by_mz[171.0069] == "C5H4ClN4O+"
    assert by_mz[230.0552] is None      # アノテーションが無いピーク


def test_neutral_loss_is_none_before_derivation():
    """世代A では mbo:neutral_loss が無い。算出は Task 5 が行う。"""
    rec = parse_gen_a(_graph())[0]
    assert all(p.neutral_loss is None for p in rec.peaks)


def test_layers_are_derived():
    rec = parse_gen_a(_graph())[0]
    assert rec.acquisition_class == "LC-ESI-MS2"
    assert rec.resolution_class == "high"
```

- [ ] **Step 4: テストが落ちることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_reader_gen_a.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'massbank_context.ttl'`

- [ ] **Step 5: 実装する**

`massbank_context/ttl/generation.py`:

```python
"""TTL のスキーマ世代を判定する（spec §3.3）。

MassBank-RDF には 2 世代のスキーマが存在する。世代A（draft_ttl/ の実体）は
ピーク値をリテラルで直付けし、世代B（SPARQL/ の例が前提とするもの）は
SIO で具象化したノードに sio:SIO_000300 で値を付ける。判定は述語の存在で行う。
"""
from __future__ import annotations

import rdflib

MBO = rdflib.Namespace("http://www.massbank.jp/ontology/")
SIO = rdflib.Namespace("http://semanticscience.org/resource/")
SCHEMA = rdflib.Namespace("https://schema.org/")

#: 世代B で値を保持する述語。
SIO_HAS_VALUE = SIO["SIO_000300"]

GEN_A = "gen_a"
GEN_B = "gen_b"
UNKNOWN = "unknown"


def detect_generation(graph: rdflib.Graph) -> str:
    """グラフの世代を返す。

    世代B は mz を rdf:type mbo:mz のノードとして持つ。世代A は
    mbo:mz をリテラル述語として使う。両者は排他的なので type 文の有無で分かれる。
    """
    if (None, rdflib.RDF.type, MBO["mz"]) in graph:
        return GEN_B
    for _ in graph.subject_objects(MBO["has_peak"]):
        return GEN_A
    return UNKNOWN
```

`massbank_context/ttl/reader_gen_a.py`:

```python
"""世代A の TTL を SpectrumRecord に変換する（spec §3.4）。

世代A の構造:
  <accession> mbo:has_peak [ mbo:mz ; mbo:intensity ; mbo:relative_intensity ]
              mbo:has_peak_annotations [ mbo:mz ; mbo:tentative_formula ]
              mbo:analytical_methods_and_conditions [ ... ]
              mbo:compound [ schema:inChIKey ; ... ]
              mbo:mass_spectram_meta_data [ mbo:precursor_mz_value ; ... ]

has_peak と has_peak_annotations は別ブロックなので m/z で突き合わせる。
アノテーションは全レコードにあるわけではない（サブミッター単位で全か無か）。
"""
from __future__ import annotations

import rdflib

from massbank_context.core.models import Peak, SpectrumRecord
from massbank_context.ttl.generation import MBO, SCHEMA

#: m/z を突き合わせる際の丸め桁。TTL 内の同一ピークは同一リテラルなので
#: 完全一致でよいが、指数表記の揺れを吸収するために丸めて鍵にする。
_MZ_KEY_DIGITS = 6


def _one(graph: rdflib.Graph, subject, predicate):
    for value in graph.objects(subject, predicate):
        return value
    return None


def _text(value) -> str | None:
    return None if value is None else str(value)


def _number(value) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _inchikey(value) -> str | None:
    """inchikey:XSDQ... という URI から末尾のキーだけを取り出す。"""
    if value is None:
        return None
    text = str(value)
    return text.rsplit("/", 1)[-1] if "/" in text else text


def parse_gen_a(graph: rdflib.Graph) -> list[SpectrumRecord]:
    """グラフに含まれる全レコードを返す（通常は 1 ファイル 1 レコード）。"""
    records: list[SpectrumRecord] = []
    for subject in graph.subjects(rdflib.RDF.type, MBO["MassSpectrum"]):
        records.append(_parse_one(graph, subject))
    records.sort(key=lambda r: r.accession)
    return records


def _parse_one(graph: rdflib.Graph, subject) -> SpectrumRecord:
    conditions = _one(graph, subject, MBO["analytical_methods_and_conditions"])
    compound = _one(graph, subject, MBO["compound"])
    meta = _one(graph, subject, MBO["mass_spectram_meta_data"])

    names = tuple(sorted(
        str(v) for v in graph.objects(compound, SCHEMA["name"])
    )) if compound is not None else ()

    formula_by_mz: dict[float, str] = {}
    for annotation in graph.objects(subject, MBO["has_peak_annotations"]):
        mz = _number(_one(graph, annotation, MBO["mz"]))
        formula = _text(_one(graph, annotation, MBO["tentative_formula"]))
        if mz is not None and formula:
            formula_by_mz[round(mz, _MZ_KEY_DIGITS)] = formula

    peaks: list[Peak] = []
    for node in graph.objects(subject, MBO["has_peak"]):
        mz = _number(_one(graph, node, MBO["mz"]))
        if mz is None:
            continue
        peaks.append(Peak(
            mz=mz,
            relative_intensity=_number(
                _one(graph, node, MBO["relative_intensity"])) or 0.0,
            formula=formula_by_mz.get(round(mz, _MZ_KEY_DIGITS)),
            neutral_loss=None,
        ))
    peaks.sort(key=lambda p: p.mz)

    license_value = _one(graph, subject, rdflib.URIRef(
        "http://purl.org/dc/terms/license"))

    return SpectrumRecord(
        accession=_text(_one(graph, subject, rdflib.URIRef(
            "http://purl.org/dc/terms/identifier"))) or str(subject).rsplit("/", 1)[-1],
        inchikey=_inchikey(_one(graph, compound, SCHEMA["inChIKey"]))
        if compound is not None else None,
        precursor_mz=_number(_one(graph, meta, MBO["precursor_mz_value"]))
        if meta is not None else None,
        precursor_type=_text(_one(graph, meta, MBO["precursor_type_value"]))
        if meta is not None else None,
        polarity=_text(_one(graph, conditions, MBO["ion_mode"]))
        if conditions is not None else None,
        instrument_type=_text(_one(graph, conditions, MBO["instrument_type"]))
        if conditions is not None else None,
        ms_type=_text(_one(graph, conditions, MBO["ms_type"]))
        if conditions is not None else None,
        resolution_value=_number(_one(graph, conditions, MBO["resolution"]))
        if conditions is not None else None,
        collision_energy=_text(_one(graph, conditions, MBO["collision_energy"]))
        if conditions is not None else None,
        splash=_text(_one(graph, subject, MBO["pk_splash"])),
        license=_text(license_value),
        molecular_formula=_text(_one(graph, compound, SCHEMA["molecularFormula"]))
        if compound is not None else None,
        exact_mass=_number(_one(graph, compound, MBO["ch_exact_mass"]))
        if compound is not None else None,
        compound_names=names,
        ch_compound_class=_text(_one(graph, compound, MBO["ch_compound_class"]))
        if compound is not None else None,
        peaks=tuple(peaks),
    )
```

`massbank_context/ttl/__init__.py` は空ファイル。

- [ ] **Step 6: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_reader_gen_a.py -q`
Expected: PASS（7 件）

- [ ] **Step 7: コミット**

```bash
cd C:/Users/yuu18/massbank-context && git add -A && git commit -m "feat: 世代判定と世代A TTL リーダ"
```

---

### Task 4: 世代B TTL リーダ

**世代B の実 TTL は入手できていない。** `SPARQL/get_neutralloss_from_inchikey.sparql` と
`get_inchikey_from_peakmz.sparql` のトリプルパターンから構造を再構成した**合成 fixture**を使う。
spec §11 の確認事項 1 が解決して実物が手に入ったら、この fixture を差し替えること。

**Files:**
- Create: `massbank_context/ttl/reader_gen_b.py`
- Create: `massbank_context/ttl/reader.py`
- Create: `tests/fixtures/ttl/gen_b/MSBNK-TEST-B0001.ttl`
- Test: `tests/test_reader_gen_b.py`

**Interfaces:**
- Consumes: `detect_generation` / `parse_gen_a`（Task 3）
- Produces:
  - `parse_gen_b(graph: rdflib.Graph) -> list[SpectrumRecord]`
  - `read_ttl_file(path: str | Path) -> list[SpectrumRecord]` — 世代を判定して振り分ける入口

- [ ] **Step 1: 合成 fixture を作る**

`tests/fixtures/ttl/gen_b/MSBNK-TEST-B0001.ttl`:

```turtle
@prefix dcterms: <http://purl.org/dc/terms/> .
@prefix massbank: <https://identifiers.org/massbank/> .
@prefix mbo: <http://www.massbank.jp/ontology/> .
@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
@prefix schema: <https://schema.org/> .
@prefix sio: <http://semanticscience.org/resource/> .

# SPARQL/get_neutralloss_from_inchikey.sparql のパターンから再構成した合成データ。
# 実物が入手できたら差し替えること（spec §11 確認事項 1）。
massbank:MSBNK-TEST-B0001 a mbo:MassSpectrum ;
    dcterms:identifier "MSBNK-TEST-B0001" ;
    mbo:analytical_methods_and_conditions [ a mbo:AnalyticalMethods ;
            mbo:instrument_type "LC-ESI-QTOF" ;
            mbo:ion_mode "NEGATIVE" ;
            mbo:ms_type "MS2" ;
            mbo:resolution "35000" ] ;
    mbo:compound [ a schema:chemicalsubstance ;
            schema:inChIKey [ a sio:CHEMINF_000059 ;
                    sio:SIO_000300 "XHKOSGCDDIGCOM-HIAJHCFCSA-N" ] ;
            schema:molecularFormula "C10H12N4O5" ;
            schema:name "TestCompoundB" ] ;
    mbo:mass_spectram_meta_data [ a mbo:MassSpectramMetaData ;
            mbo:precursor_mz_value "267.0855"^^<http://www.w3.org/2001/XMLSchema#float> ;
            mbo:precursor_type_value "[M-H]-" ] ;
    mbo:has_peak [ sio:SIO_000255 [
                sio:SIO_000255 [ a mbo:mz ; sio:SIO_000300 1.350434e+02 ],
                                [ a mbo:relative_intensity ; sio:SIO_000300 999 ],
                                [ a mbo:neutral_loss ; sio:SIO_000300 1.320421e+02 ] ] ],
        [ sio:SIO_000255 [
                sio:SIO_000255 [ a mbo:mz ; sio:SIO_000300 1.100349e+02 ],
                                [ a mbo:relative_intensity ; sio:SIO_000300 250 ],
                                [ a mbo:neutral_loss ; sio:SIO_000300 1.570506e+02 ] ] ] ;
    mbo:pk_num_peak 2 .
```

- [ ] **Step 2: 失敗するテストを書く**

`tests/test_reader_gen_b.py`:

```python
from pathlib import Path

import rdflib

from massbank_context.ttl.generation import detect_generation
from massbank_context.ttl.reader import read_ttl_file
from massbank_context.ttl.reader_gen_b import parse_gen_b

FIXTURES = Path(__file__).parent / "fixtures" / "ttl"
GEN_B = FIXTURES / "gen_b" / "MSBNK-TEST-B0001.ttl"
GEN_A = FIXTURES / "gen_a" / "MSBNK-TEST-A0001.ttl"


def _graph(path) -> rdflib.Graph:
    g = rdflib.Graph()
    g.parse(path, format="turtle")
    return g


def test_detects_generation_b():
    assert detect_generation(_graph(GEN_B)) == "gen_b"


def test_parses_metadata():
    rec = parse_gen_b(_graph(GEN_B))[0]
    assert rec.accession == "MSBNK-TEST-B0001"
    assert rec.inchikey == "XHKOSGCDDIGCOM-HIAJHCFCSA-N"
    assert rec.precursor_mz == 267.0855
    assert rec.polarity == "NEGATIVE"
    assert rec.acquisition_class == "LC-ESI-MS2"


def test_reads_neutral_loss_directly():
    """世代B は mbo:neutral_loss を持つので算出不要（spec §6.3）。"""
    rec = parse_gen_b(_graph(GEN_B))[0]
    by_mz = {round(p.mz, 4): p.neutral_loss for p in rec.peaks}
    assert by_mz[135.0434] == 132.0421
    assert by_mz[110.0349] == 157.0506


def test_peaks_sorted_and_complete():
    rec = parse_gen_b(_graph(GEN_B))[0]
    mzs = [p.mz for p in rec.peaks]
    assert mzs == sorted(mzs)
    assert len(mzs) == 2


def test_read_ttl_file_dispatches_gen_a():
    records = read_ttl_file(GEN_A)
    assert records[0].accession == "MSBNK-TEST-A0001"


def test_read_ttl_file_dispatches_gen_b():
    records = read_ttl_file(GEN_B)
    assert records[0].accession == "MSBNK-TEST-B0001"


def test_read_ttl_file_returns_empty_for_unknown(tmp_path):
    path = tmp_path / "empty.ttl"
    path.write_text("@prefix ex: <http://example.org/> .\n", encoding="utf-8")
    assert read_ttl_file(path) == []
```

- [ ] **Step 3: テストが落ちることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_reader_gen_b.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'massbank_context.ttl.reader_gen_b'`

- [ ] **Step 4: 実装する**

`massbank_context/ttl/reader_gen_b.py`:

```python
"""世代B の TTL を SpectrumRecord に変換する（spec §3.3）。

世代B は値を SIO で具象化する:
  <accession> mbo:has_peak ?peak
  ?peak  sio:SIO_000255 ?annotation
  ?annotation sio:SIO_000255 ?content
  ?content a mbo:mz | mbo:relative_intensity | mbo:intensity | mbo:neutral_loss
  ?content sio:SIO_000300 ?value

世代A と違い、中性損失が既に入っているので算出しない。
"""
from __future__ import annotations

import rdflib

from massbank_context.core.models import Peak, SpectrumRecord
from massbank_context.ttl.generation import MBO, SCHEMA, SIO, SIO_HAS_VALUE

SIO_HAS_ATTRIBUTE = SIO["SIO_000255"]
CHEMINF_INCHIKEY = SIO["CHEMINF_000059"]


def _one(graph: rdflib.Graph, subject, predicate):
    for value in graph.objects(subject, predicate):
        return value
    return None


def _text(value) -> str | None:
    return None if value is None else str(value)


def _number(value) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _typed_value(graph: rdflib.Graph, annotation, type_uri):
    """?annotation sio:SIO_000255 [ a <type_uri> ; sio:SIO_000300 ?v ] の ?v。"""
    for content in graph.objects(annotation, SIO_HAS_ATTRIBUTE):
        if (content, rdflib.RDF.type, type_uri) in graph:
            return _one(graph, content, SIO_HAS_VALUE)
    return None


def _inchikey(graph: rdflib.Graph, compound) -> str | None:
    """世代B では inChIKey も具象化ノードになりうる。両形を受ける。"""
    node = _one(graph, compound, SCHEMA["inChIKey"])
    if node is None:
        return None
    if isinstance(node, rdflib.Literal):
        return str(node)
    value = _one(graph, node, SIO_HAS_VALUE)
    if value is not None:
        return str(value)
    text = str(node)
    return text.rsplit("/", 1)[-1] if "/" in text else text


def parse_gen_b(graph: rdflib.Graph) -> list[SpectrumRecord]:
    records: list[SpectrumRecord] = []
    for subject in graph.subjects(rdflib.RDF.type, MBO["MassSpectrum"]):
        records.append(_parse_one(graph, subject))
    records.sort(key=lambda r: r.accession)
    return records


def _parse_one(graph: rdflib.Graph, subject) -> SpectrumRecord:
    conditions = _one(graph, subject, MBO["analytical_methods_and_conditions"])
    compound = _one(graph, subject, MBO["compound"])
    meta = _one(graph, subject, MBO["mass_spectram_meta_data"])

    peaks: list[Peak] = []
    for peak_node in graph.objects(subject, MBO["has_peak"]):
        for annotation in graph.objects(peak_node, SIO_HAS_ATTRIBUTE):
            mz = _number(_typed_value(graph, annotation, MBO["mz"]))
            if mz is None:
                continue
            peaks.append(Peak(
                mz=mz,
                relative_intensity=_number(_typed_value(
                    graph, annotation, MBO["relative_intensity"])) or 0.0,
                formula=_text(_typed_value(
                    graph, annotation, MBO["tentative_formula"])),
                neutral_loss=_number(_typed_value(
                    graph, annotation, MBO["neutral_loss"])),
            ))
    peaks.sort(key=lambda p: p.mz)

    return SpectrumRecord(
        accession=_text(_one(graph, subject, rdflib.URIRef(
            "http://purl.org/dc/terms/identifier"))) or str(subject).rsplit("/", 1)[-1],
        inchikey=_inchikey(graph, compound) if compound is not None else None,
        precursor_mz=_number(_one(graph, meta, MBO["precursor_mz_value"]))
        if meta is not None else None,
        precursor_type=_text(_one(graph, meta, MBO["precursor_type_value"]))
        if meta is not None else None,
        polarity=_text(_one(graph, conditions, MBO["ion_mode"]))
        if conditions is not None else None,
        instrument_type=_text(_one(graph, conditions, MBO["instrument_type"]))
        if conditions is not None else None,
        ms_type=_text(_one(graph, conditions, MBO["ms_type"]))
        if conditions is not None else None,
        resolution_value=_number(_one(graph, conditions, MBO["resolution"]))
        if conditions is not None else None,
        collision_energy=_text(_one(graph, conditions, MBO["collision_energy"]))
        if conditions is not None else None,
        splash=_text(_one(graph, subject, MBO["pk_splash"])),
        license=_text(_one(graph, subject, rdflib.URIRef(
            "http://purl.org/dc/terms/license"))),
        molecular_formula=_text(_one(graph, compound, SCHEMA["molecularFormula"]))
        if compound is not None else None,
        exact_mass=_number(_one(graph, compound, MBO["ch_exact_mass"]))
        if compound is not None else None,
        compound_names=tuple(sorted(
            str(v) for v in graph.objects(compound, SCHEMA["name"])
        )) if compound is not None else (),
        ch_compound_class=_text(_one(graph, compound, MBO["ch_compound_class"]))
        if compound is not None else None,
        peaks=tuple(peaks),
    )
```

`massbank_context/ttl/reader.py`:

```python
"""TTL ファイルを読む入口。世代を判定して適切なリーダへ振り分ける。"""
from __future__ import annotations

from pathlib import Path

import rdflib

from massbank_context.core.models import SpectrumRecord
from massbank_context.ttl.generation import GEN_A, GEN_B, detect_generation
from massbank_context.ttl.reader_gen_a import parse_gen_a
from massbank_context.ttl.reader_gen_b import parse_gen_b


def read_ttl_file(path: str | Path) -> list[SpectrumRecord]:
    """1 ファイルを読む。世代が判定できない場合は空リストを返す。

    索引ビルドは数万ファイルを走査するため、1 件の異常でビルド全体を止めない。
    呼び出し側（CLI）が件数を集計して報告する。
    """
    graph = rdflib.Graph()
    graph.parse(str(path), format="turtle")
    generation = detect_generation(graph)
    if generation == GEN_A:
        return parse_gen_a(graph)
    if generation == GEN_B:
        return parse_gen_b(graph)
    return []
```

- [ ] **Step 5: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_reader_gen_b.py -q`
Expected: PASS（7 件）

- [ ] **Step 6: コミット**

```bash
cd C:/Users/yuu18/massbank-context && git add -A && git commit -m "feat: 世代B TTL リーダと世代振り分けの入口"
```

---

### Task 5: ビン化と中性損失の算出

ビン境界でピークが分断される問題を、**`round(mz/w)` で索引し、照会時に ±1 ビンを探る**
ことで回避する。境界損失が無く、索引を二重化する必要もない。

**Files:**
- Create: `massbank_context/core/binning.py`
- Test: `tests/test_binning.py`

**Interfaces:**
- Consumes: `Peak` / `SpectrumRecord`（Task 2）、`NEUTRAL_LOSS_CAPABLE`（Task 1）
- Produces:
  - `to_bin(mz: float, bin_width: float) -> int`
  - `probe_bins(mz: float, bin_width: float) -> tuple[int, int, int]` — 照会時に探る 3 ビン
  - `derive_neutral_losses(record: SpectrumRecord) -> tuple[Peak, ...]`
  - `axis_entries(record: SpectrumRecord, bin_width: float) -> list[tuple[str, int, float, str | None]]`
    — `(axis, bin, relative_intensity, formula)` の列

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_binning.py`:

```python
import pytest

from massbank_context.core.binning import (
    axis_entries,
    derive_neutral_losses,
    probe_bins,
    to_bin,
)
from massbank_context.core.models import Peak, SpectrumRecord


def _record(**overrides) -> SpectrumRecord:
    base = dict(
        accession="A", inchikey="IK", precursor_mz=230.0552,
        precursor_type="[M+H]+", polarity="POSITIVE",
        instrument_type="LC-ESI-ITFT", ms_type="MS2", resolution_value=15000.0,
        collision_energy=None, splash=None, license=None,
        molecular_formula=None, exact_mass=None,
        peaks=(Peak(86.0349, 11.0, "C2H4N3O+", None),
               Peak(171.0069, 226.0, None, None)),
    )
    base.update(overrides)
    return SpectrumRecord(**base)


def test_to_bin_rounds_to_nearest():
    assert to_bin(86.0349, 0.01) == 8603
    assert to_bin(86.0351, 0.01) == 8604


def test_probe_bins_covers_neighbours():
    assert probe_bins(86.0349, 0.01) == (8602, 8603, 8604)


def test_bin_boundary_is_recoverable():
    """境界を挟んだ 2 つの m/z が、互いの探索範囲に入ること。"""
    a, b = 86.0349, 86.0351          # 別ビンに落ちる近接ペア
    assert to_bin(a, 0.01) != to_bin(b, 0.01)
    assert to_bin(b, 0.01) in probe_bins(a, 0.01)
    assert to_bin(a, 0.01) in probe_bins(b, 0.01)


def test_derive_neutral_losses_for_lc_esi():
    peaks = derive_neutral_losses(_record())
    losses = [round(p.neutral_loss, 4) for p in peaks]
    assert losses == [144.0203, 59.0483]


def test_derive_neutral_losses_skips_gc_ei():
    """GC-EI-MS は precursor を持たないので中性損失は None のまま（spec §6.1）。"""
    rec = _record(instrument_type="EI-B", ms_type="MS",
                  precursor_mz=None, precursor_type=None)
    peaks = derive_neutral_losses(rec)
    assert all(p.neutral_loss is None for p in peaks)


def test_derive_neutral_losses_keeps_existing_values():
    """世代B が読み込んだ値を上書きしない。"""
    rec = _record(peaks=(Peak(135.0434, 999.0, None, 132.0421),))
    peaks = derive_neutral_losses(rec)
    assert peaks[0].neutral_loss == 132.0421


def test_axis_entries_emits_both_axes():
    entries = axis_entries(_record(), bin_width=0.01)
    axes = {axis for axis, _bin, _intensity, _formula in entries}
    assert axes == {"product_ion", "neutral_loss"}
    assert len(entries) == 4      # ピーク 2 本 × 2 軸


def test_axis_entries_omits_neutral_loss_for_gc_ei():
    rec = _record(instrument_type="EI-B", ms_type="MS",
                  precursor_mz=None, precursor_type=None)
    entries = axis_entries(rec, bin_width=0.01)
    axes = {axis for axis, _bin, _intensity, _formula in entries}
    assert axes == {"product_ion"}


def test_negative_neutral_loss_is_dropped():
    """precursor より重いピーク（同位体・多量体）は中性損失軸に載せない。"""
    rec = _record(peaks=(Peak(300.0, 10.0, None, None),))
    entries = axis_entries(rec, bin_width=0.01)
    axes = {axis for axis, _bin, _intensity, _formula in entries}
    assert axes == {"product_ion"}


@pytest.mark.parametrize("bin_width", [0.0, -0.01])
def test_invalid_bin_width_raises(bin_width):
    with pytest.raises(ValueError):
        to_bin(100.0, bin_width)
```

- [ ] **Step 2: テストが落ちることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_binning.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'massbank_context.core.binning'`

- [ ] **Step 3: 実装する**

`massbank_context/core/binning.py`:

```python
"""m/z のビン化と中性損失の算出（spec §6.3 / §6.4）。

固定幅ビンには境界問題がある。近接した 2 つの m/z が境界を挟むと別ビンに落ち、
同一フラグメントが一致しなくなる。ここでは round() で最近傍ビンに落とし、
照会時に ±1 ビンを探ることで境界損失を無くす。索引を二重化する必要がない。

依存グラフの leaf（classify / models のみ）。
"""
from __future__ import annotations

from massbank_context.core.classify import NEUTRAL_LOSS_CAPABLE
from massbank_context.core.models import Peak, SpectrumRecord

PRODUCT_ION = "product_ion"
NEUTRAL_LOSS = "neutral_loss"

#: 中性損失としてこれ未満は扱わない。precursor 自身のピーク（loss≒0）と
#: 同位体ピークを軸に載せないため。
_MIN_NEUTRAL_LOSS = 0.5


def to_bin(mz: float, bin_width: float) -> int:
    """m/z を最近傍のビン番号にする。"""
    if bin_width <= 0:
        raise ValueError(f"bin_width は正の数でなければならない: {bin_width}")
    return round(mz / bin_width)


def probe_bins(mz: float, bin_width: float) -> tuple[int, int, int]:
    """照会時に探るビン番号。境界を挟んだ相手を取りこぼさないための ±1。"""
    center = to_bin(mz, bin_width)
    return (center - 1, center, center + 1)


def derive_neutral_losses(record: SpectrumRecord) -> tuple[Peak, ...]:
    """中性損失が未設定のピークについて precursor_mz - mz を埋める。

    世代B が読み込んだ値は上書きしない。precursor を持たない acquisition_class
    （GC-EI-MS）では何もしない。
    """
    if (record.acquisition_class not in NEUTRAL_LOSS_CAPABLE
            or record.precursor_mz is None):
        return record.peaks

    filled: list[Peak] = []
    for peak in record.peaks:
        if peak.neutral_loss is not None:
            filled.append(peak)
            continue
        loss = record.precursor_mz - peak.mz
        filled.append(Peak(
            mz=peak.mz,
            relative_intensity=peak.relative_intensity,
            formula=peak.formula,
            neutral_loss=loss if loss >= _MIN_NEUTRAL_LOSS else None,
        ))
    return tuple(filled)


def axis_entries(record: SpectrumRecord,
                 bin_width: float) -> list[tuple[str, int, float, str | None]]:
    """レコードを (axis, bin, relative_intensity, formula) の列に展開する。"""
    peaks = derive_neutral_losses(record)
    entries: list[tuple[str, int, float, str | None]] = []
    for peak in peaks:
        entries.append((PRODUCT_ION, to_bin(peak.mz, bin_width),
                        peak.relative_intensity, peak.formula))
        if peak.neutral_loss is not None:
            entries.append((NEUTRAL_LOSS, to_bin(peak.neutral_loss, bin_width),
                            peak.relative_intensity, peak.formula))
    return entries
```

- [ ] **Step 4: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_binning.py -q`
Expected: PASS（11 件）

- [ ] **Step 5: コミット**

```bash
cd C:/Users/yuu18/massbank-context && git add -A && git commit -m "feat: ビン化（最近傍＋±1 探索）と中性損失の算出"
```

---

### Task 6: 質量誤差プロファイルとビン幅の決定

spec §6.4。`0.01 Da` のような決め打ちを避け、**同一 InChIKey の複数レコード間で
同じフラグメントの m/z がどう散らばるか**を実測してビン幅を決める。

**Files:**
- Create: `massbank_context/index/__init__.py`
- Create: `massbank_context/index/mass_error.py`
- Test: `tests/test_mass_error.py`

**Interfaces:**
- Consumes: `SpectrumRecord`（Task 2）
- Produces:
  - `measure_mass_error(records: list[SpectrumRecord], cluster_tolerance: float = 0.05) -> dict[str, list[tuple[float, float]]]`
    — `{instrument_class: [(percentile, delta_mz), ...]}`。percentile は 50/90/95/99
  - `choose_bin_width(profile: dict[str, list[tuple[float, float]]], percentile: float = 95.0, fallback: float = 0.01) -> float`
  - `MASS_ERROR_PERCENTILES: tuple[float, ...]`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_mass_error.py`:

```python
from massbank_context.core.models import Peak, SpectrumRecord
from massbank_context.index.mass_error import (
    MASS_ERROR_PERCENTILES,
    choose_bin_width,
    measure_mass_error,
)


def _record(accession: str, inchikey: str, mzs: list[float],
            instrument_type: str = "LC-ESI-QTOF") -> SpectrumRecord:
    return SpectrumRecord(
        accession=accession, inchikey=inchikey, precursor_mz=300.0,
        precursor_type="[M+H]+", polarity="POSITIVE",
        instrument_type=instrument_type, ms_type="MS2", resolution_value=None,
        collision_energy=None, splash=None, license=None,
        molecular_formula=None, exact_mass=None,
        peaks=tuple(Peak(mz, 100.0, None, None) for mz in mzs),
    )


def test_measures_dispersion_within_same_inchikey():
    """同じ InChIKey の 2 レコードで同じフラグメントが 0.004 ずれている場合。"""
    records = [
        _record("A1", "IK1", [100.000, 200.000]),
        _record("A2", "IK1", [100.004, 200.004]),
    ]
    profile = measure_mass_error(records)
    assert "LC-ESI-QTOF" in profile
    deltas = dict(profile["LC-ESI-QTOF"])
    assert abs(deltas[95.0] - 0.004) < 1e-6


def test_ignores_singleton_inchikeys():
    """レコードが 1 件しかない InChIKey からは散らばりを測れない。"""
    profile = measure_mass_error([_record("A1", "IK1", [100.0])])
    assert profile == {}


def test_separates_instrument_classes():
    records = [
        _record("A1", "IK1", [100.000], "LC-ESI-QTOF"),
        _record("A2", "IK1", [100.002], "LC-ESI-QTOF"),
        _record("B1", "IK2", [100.000], "LC-ESI-QQ"),
        _record("B2", "IK2", [100.400], "LC-ESI-QQ"),
    ]
    profile = measure_mass_error(records)
    assert set(profile) == {"LC-ESI-QTOF", "LC-ESI-QQ"}
    assert dict(profile["LC-ESI-QTOF"])[95.0] < dict(profile["LC-ESI-QQ"])[95.0]


def test_distant_peaks_are_not_clustered():
    """0.05 Da 以上離れたピークは同一フラグメントとみなさない。"""
    records = [
        _record("A1", "IK1", [100.000]),
        _record("A2", "IK1", [180.000]),
    ]
    assert measure_mass_error(records) == {}


def test_choose_bin_width_uses_max_across_classes():
    """索引は 1 つのビン幅で作るので、最も散らばる装置に合わせる。"""
    profile = {
        "LC-ESI-QTOF": [(95.0, 0.004)],
        "LC-ESI-ITFT": [(95.0, 0.002)],
    }
    assert choose_bin_width(profile) == 0.004


def test_choose_bin_width_falls_back_when_empty():
    assert choose_bin_width({}) == 0.01


def test_percentiles_are_declared():
    assert MASS_ERROR_PERCENTILES == (50.0, 90.0, 95.0, 99.0)
```

- [ ] **Step 2: テストが落ちることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_mass_error.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'massbank_context.index'`

- [ ] **Step 3: 実装する**

`massbank_context/index/__init__.py` は空ファイル。

`massbank_context/index/mass_error.py`:

```python
"""質量誤差の実測とビン幅の決定（spec §6.4）。

決め打ちのビン幅には根拠がない。MassBank は装置が混在しており、必要なトレランスは
レコード集団によって違う。同一 InChIKey の複数レコード（装置・CE 違い）で同じ
フラグメントの m/z がどう散らばるかを測れば、必要な幅が実測で出る。
"""
from __future__ import annotations

import statistics
from collections import defaultdict

from massbank_context.core.models import SpectrumRecord

MASS_ERROR_PERCENTILES = (50.0, 90.0, 95.0, 99.0)

#: これ以上離れたピークは別フラグメントとみなす。装置間差より十分大きく、
#: 隣接フラグメント（多くは 1 Da 以上離れる）より十分小さい値。
DEFAULT_CLUSTER_TOLERANCE = 0.05

DEFAULT_BIN_WIDTH = 0.01


def _percentile(values: list[float], percentile: float) -> float:
    """線形補間なしの単純な順位統計。値が少ないときも壊れない。"""
    if not values:
        return 0.0
    ordered = sorted(values)
    index = min(len(ordered) - 1,
                max(0, int(round(percentile / 100.0 * (len(ordered) - 1)))))
    return ordered[index]


def _instrument_class(record: SpectrumRecord) -> str:
    return record.instrument_type or "unknown"


def measure_mass_error(
    records: list[SpectrumRecord],
    cluster_tolerance: float = DEFAULT_CLUSTER_TOLERANCE,
) -> dict[str, list[tuple[float, float]]]:
    """装置種別ごとの m/z 散らばりを返す。

    同一 (instrument_type, inchikey) のレコード群を集め、各ピークを
    cluster_tolerance 内で貪欲にクラスタリングして、2 件以上のクラスタの
    最大幅（max - min）を散らばりの標本とする。
    """
    grouped: dict[tuple[str, str], list[list[float]]] = defaultdict(list)
    for record in records:
        if not record.inchikey:
            continue
        grouped[(_instrument_class(record), record.inchikey)].append(
            sorted(peak.mz for peak in record.peaks))

    spreads: dict[str, list[float]] = defaultdict(list)
    for (instrument, _inchikey), spectra in grouped.items():
        if len(spectra) < 2:
            continue
        all_mz = sorted(mz for spectrum in spectra for mz in spectrum)
        cluster: list[float] = []
        for mz in all_mz:
            if cluster and mz - cluster[0] > cluster_tolerance:
                if len(cluster) >= 2:
                    spreads[instrument].append(cluster[-1] - cluster[0])
                cluster = []
            cluster.append(mz)
        if len(cluster) >= 2:
            spreads[instrument].append(cluster[-1] - cluster[0])

    profile: dict[str, list[tuple[float, float]]] = {}
    for instrument, values in spreads.items():
        if not values:
            continue
        profile[instrument] = [
            (percentile, _percentile(values, percentile))
            for percentile in MASS_ERROR_PERCENTILES
        ]
    return profile


def choose_bin_width(profile: dict[str, list[tuple[float, float]]],
                     percentile: float = 95.0,
                     fallback: float = DEFAULT_BIN_WIDTH) -> float:
    """索引全体で使う 1 つのビン幅を決める。

    索引は単一のビン幅で作るため、最も散らばる装置に合わせる。
    測定できなかった場合は fallback を返し、index_meta に記録する。
    """
    candidates = [
        value
        for entries in profile.values()
        for pct, value in entries
        if pct == percentile and value > 0
    ]
    if not candidates:
        return fallback
    return max(candidates)
```

- [ ] **Step 4: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_mass_error.py -q`
Expected: PASS（7 件）

- [ ] **Step 5: コミット**

```bash
cd C:/Users/yuu18/massbank-context && git add -A && git commit -m "feat: 質量誤差プロファイルの実測とビン幅の決定"
```

---

### Task 7: SQLite スキーマと書き込み

**Files:**
- Create: `massbank_context/index/schema.py`
- Create: `massbank_context/index/writer.py`
- Test: `tests/test_index_writer.py`

**Interfaces:**
- Consumes: `SpectrumRecord`（Task 2）、`axis_entries`（Task 5）
- Produces:
  - `create_schema(conn: sqlite3.Connection) -> None`
  - `IndexWriter(conn, bin_width)` — メソッド `add_record(record) -> bool` / `finalize(meta: dict) -> None`
  - `open_index(path: str | Path, create: bool = False) -> sqlite3.Connection`
  - `read_meta(conn) -> dict | None`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_index_writer.py`:

```python
import sqlite3

import pytest

from massbank_context.core.models import Peak, SpectrumRecord
from massbank_context.index.writer import IndexWriter, create_schema, read_meta


def _record(accession="A1", inchikey="IK1", instrument_type="LC-ESI-QTOF",
            ms_type="MS2", precursor_mz=230.0552) -> SpectrumRecord:
    return SpectrumRecord(
        accession=accession, inchikey=inchikey, precursor_mz=precursor_mz,
        precursor_type="[M+H]+", polarity="POSITIVE",
        instrument_type=instrument_type, ms_type=ms_type, resolution_value=35000.0,
        collision_energy="35", splash="splash10-xxx", license="CC0",
        molecular_formula="C6H8ClN7O", exact_mass=229.0479,
        compound_names=("Amiloride",), ch_compound_class="Natural Product; Lipid",
        peaks=(Peak(86.0349, 11.0, "C2H4N3O+", None),
               Peak(171.0069, 226.0, None, None)),
    )


@pytest.fixture()
def conn():
    connection = sqlite3.connect(":memory:")
    create_schema(connection)
    yield connection
    connection.close()


def test_schema_creates_all_tables(conn):
    names = {row[0] for row in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    assert {"records", "peaks", "axis_background", "compound_class",
            "compounds", "index_meta", "mass_error_profile"} <= names


def test_add_record_stores_layers(conn):
    writer = IndexWriter(conn, bin_width=0.01)
    assert writer.add_record(_record()) is True
    row = conn.execute(
        "SELECT acquisition_class, resolution_class FROM records").fetchone()
    assert row == ("LC-ESI-MS2", "high")


def test_add_record_writes_both_axes(conn):
    writer = IndexWriter(conn, bin_width=0.01)
    writer.add_record(_record())
    axes = {row[0] for row in conn.execute("SELECT DISTINCT axis FROM peaks")}
    assert axes == {"product_ion", "neutral_loss"}


def test_low_resolution_is_rejected_by_default(conn):
    """既定は resolution_class == 'high' のみ（spec §6.2）。"""
    writer = IndexWriter(conn, bin_width=0.01)
    assert writer.add_record(_record(instrument_type="LC-ESI-QQ")) is False
    assert conn.execute("SELECT COUNT(*) FROM records").fetchone()[0] == 0


def test_other_acquisition_class_is_rejected(conn):
    writer = IndexWriter(conn, bin_width=0.01)
    assert writer.add_record(_record(ms_type="MS3")) is False


def test_record_without_inchikey_is_rejected(conn):
    """InChIKey が無いレコードは IDF の分母にも候補にもならない。"""
    writer = IndexWriter(conn, bin_width=0.01)
    assert writer.add_record(_record(inchikey=None)) is False


def test_compound_class_is_split_on_semicolon(conn):
    writer = IndexWriter(conn, bin_width=0.01)
    writer.add_record(_record())
    labels = {row[0] for row in conn.execute(
        "SELECT label FROM compound_class")}
    assert labels == {"Natural Product", "Lipid"}


def test_na_compound_class_is_not_stored(conn):
    writer = IndexWriter(conn, bin_width=0.01)
    record = _record()
    writer.add_record(SpectrumRecord(**{
        **record.__dict__, "ch_compound_class": "N/A; Environmental Standard"}))
    labels = {row[0] for row in conn.execute("SELECT label FROM compound_class")}
    assert "N/A" not in labels


def test_finalize_writes_meta(conn):
    writer = IndexWriter(conn, bin_width=0.004)
    writer.add_record(_record())
    writer.finalize({"ttl_source": "/tmp/ttl", "ttl_generation": "gen_a",
                     "resolution_filter": "high"})
    meta = read_meta(conn)
    assert meta["bin_width"] == 0.004
    assert meta["n_records"] == 1
    assert meta["n_inchikey"] == 1
    assert meta["acquisition_classes"] == ["LC-ESI-MS2"]
    assert meta["built_at"]


def test_read_meta_returns_none_when_absent(conn):
    assert read_meta(conn) is None
```

- [ ] **Step 2: テストが落ちることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_index_writer.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'massbank_context.index.writer'`

- [ ] **Step 3: スキーマを書く**

`massbank_context/index/schema.py`:

```python
"""SQLite 索引の DDL（spec §6.7）。

SQLite を選ぶのは、msdial_repository_catalog が同じ思想（ローカル・単一ファイル・
provenance 付き）で SQLite を採っており、catalog を読むコードと流儀を揃えられるため。
54,183 レコード × 平均数十ピーク ≒ 数百万行は SQLite の得意領域。
"""

DDL = """
CREATE TABLE IF NOT EXISTS records (
    record_id         INTEGER PRIMARY KEY,
    accession         TEXT NOT NULL UNIQUE,
    inchikey          TEXT NOT NULL,
    precursor_mz      REAL,
    precursor_type    TEXT,
    polarity          TEXT,
    instrument_type   TEXT,
    ms_type           TEXT,
    acquisition_class TEXT NOT NULL,
    resolution_class  TEXT NOT NULL,
    resolution_value  REAL,
    collision_energy  TEXT,
    splash            TEXT,
    license           TEXT
);

CREATE TABLE IF NOT EXISTS peaks (
    record_id          INTEGER NOT NULL,
    axis               TEXT NOT NULL,
    bin                INTEGER NOT NULL,
    relative_intensity REAL,
    formula            TEXT
);

CREATE TABLE IF NOT EXISTS axis_background (
    acquisition_class TEXT NOT NULL,
    axis              TEXT NOT NULL,
    bin               INTEGER NOT NULL,
    df_inchikey       INTEGER NOT NULL,
    idf               REAL NOT NULL,
    PRIMARY KEY (acquisition_class, axis, bin)
);

CREATE TABLE IF NOT EXISTS compound_class (
    inchikey TEXT NOT NULL,
    source   TEXT NOT NULL,
    label    TEXT NOT NULL,
    PRIMARY KEY (inchikey, source, label)
);

CREATE TABLE IF NOT EXISTS compounds (
    inchikey          TEXT PRIMARY KEY,
    molecular_formula TEXT,
    exact_mass        REAL,
    name              TEXT
);

CREATE TABLE IF NOT EXISTS index_meta (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS mass_error_profile (
    instrument_class TEXT NOT NULL,
    percentile       REAL NOT NULL,
    delta_mz         REAL NOT NULL,
    PRIMARY KEY (instrument_class, percentile)
);

-- 照合は precursor m/z の範囲絞り込みから始まる（spec §6.8 手順 2）
CREATE INDEX IF NOT EXISTS idx_records_precursor
    ON records (acquisition_class, polarity, precursor_mz);
-- 一致集合の算出はビンからレコードを引く
CREATE INDEX IF NOT EXISTS idx_peaks_bin
    ON peaks (axis, bin, record_id);
CREATE INDEX IF NOT EXISTS idx_peaks_record
    ON peaks (record_id);
"""
```

- [ ] **Step 4: writer を書く**

`massbank_context/index/writer.py`:

```python
"""SpectrumRecord を SQLite 索引へ書き込む。

既定では resolution_class == "high" かつ acquisition_class != "other" かつ
InChIKey を持つレコードだけを受け入れる。弾いた件数は呼び出し側が集計する。
"""
from __future__ import annotations

import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from massbank_context.core.binning import axis_entries
from massbank_context.core.models import SpectrumRecord
from massbank_context.index.schema import DDL

#: 化合物クラスとして意味をなさないラベル。実データに頻出する（spec §6.6）。
_USELESS_CLASS_LABELS = frozenset({"N/A", "NA", "", "Environmental Standard"})

_META_KEYS_JSON = frozenset({"acquisition_classes"})


def create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(DDL)
    conn.commit()


def open_index(path: str | Path, create: bool = False) -> sqlite3.Connection:
    """索引を開く。create=True なら無ければ作る。"""
    path = Path(path)
    if not create and not path.exists():
        raise FileNotFoundError(path)
    conn = sqlite3.connect(str(path))
    if create:
        create_schema(conn)
    return conn


def read_meta(conn: sqlite3.Connection) -> dict | None:
    """index_meta を dict にして返す。空なら None。"""
    try:
        rows = conn.execute("SELECT key, value FROM index_meta").fetchall()
    except sqlite3.OperationalError:
        return None
    if not rows:
        return None
    meta: dict = {}
    for key, value in rows:
        if key in _META_KEYS_JSON:
            meta[key] = json.loads(value)
        elif key in {"n_records", "n_inchikey"}:
            meta[key] = int(value)
        elif key == "bin_width":
            meta[key] = float(value)
        else:
            meta[key] = value
    return meta


class IndexWriter:
    """索引への書き込みを 1 箇所に集約する。"""

    def __init__(self, conn: sqlite3.Connection, bin_width: float,
                 resolution_filter: str = "high") -> None:
        self._conn = conn
        self._bin_width = bin_width
        self._resolution_filter = resolution_filter
        self._acquisition_classes: set[str] = set()
        self._inchikeys: set[str] = set()
        self._n_records = 0

    def add_record(self, record: SpectrumRecord) -> bool:
        """受け入れたら True、層の条件で弾いたら False を返す。"""
        if not record.inchikey:
            return False
        if record.acquisition_class == "other":
            return False
        if (self._resolution_filter
                and record.resolution_class != self._resolution_filter):
            return False

        cursor = self._conn.execute(
            "INSERT OR IGNORE INTO records ("
            " accession, inchikey, precursor_mz, precursor_type, polarity,"
            " instrument_type, ms_type, acquisition_class, resolution_class,"
            " resolution_value, collision_energy, splash, license)"
            " VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (record.accession, record.inchikey, record.precursor_mz,
             record.precursor_type, record.polarity, record.instrument_type,
             record.ms_type, record.acquisition_class, record.resolution_class,
             record.resolution_value, record.collision_energy, record.splash,
             record.license))
        if cursor.rowcount == 0:      # 同一 accession が既にある
            return False
        record_id = cursor.lastrowid

        self._conn.executemany(
            "INSERT INTO peaks (record_id, axis, bin, relative_intensity, formula)"
            " VALUES (?,?,?,?,?)",
            [(record_id, axis, bin_no, intensity, formula)
             for axis, bin_no, intensity, formula
             in axis_entries(record, self._bin_width)])

        self._conn.execute(
            "INSERT OR IGNORE INTO compounds"
            " (inchikey, molecular_formula, exact_mass, name) VALUES (?,?,?,?)",
            (record.inchikey, record.molecular_formula, record.exact_mass,
             record.compound_names[0] if record.compound_names else None))

        for label in self._class_labels(record.ch_compound_class):
            self._conn.execute(
                "INSERT OR IGNORE INTO compound_class (inchikey, source, label)"
                " VALUES (?,?,?)", (record.inchikey, "massbank", label))

        self._acquisition_classes.add(record.acquisition_class)
        self._inchikeys.add(record.inchikey)
        self._n_records += 1
        return True

    @staticmethod
    def _class_labels(raw: str | None) -> list[str]:
        """'Natural Product; Lipid; Fatty acid' を個別ラベルに割る。

        'N/A' 等の無意味なラベルは落とす。層別背景の質に直結する（spec §6.6）。
        """
        if not raw:
            return []
        labels = []
        for part in raw.split(";"):
            label = part.strip()
            if label and label not in _USELESS_CLASS_LABELS:
                labels.append(label)
        return labels

    def write_mass_error_profile(
            self, profile: dict[str, list[tuple[float, float]]]) -> None:
        for instrument_class, entries in profile.items():
            for percentile, delta_mz in entries:
                self._conn.execute(
                    "INSERT OR REPLACE INTO mass_error_profile"
                    " (instrument_class, percentile, delta_mz) VALUES (?,?,?)",
                    (instrument_class, percentile, delta_mz))

    def finalize(self, meta: dict) -> None:
        """index_meta を書いてコミットする。"""
        values = {
            "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "bin_width": repr(self._bin_width),
            "n_records": str(self._n_records),
            "n_inchikey": str(len(self._inchikeys)),
            "acquisition_classes": json.dumps(
                sorted(self._acquisition_classes)),
            "resolution_filter": self._resolution_filter,
        }
        values.update({k: str(v) for k, v in meta.items()})
        for key, value in values.items():
            self._conn.execute(
                "INSERT OR REPLACE INTO index_meta (key, value) VALUES (?,?)",
                (key, value))
        self._conn.commit()
```

- [ ] **Step 5: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_index_writer.py -q`
Expected: PASS（10 件）

- [ ] **Step 6: コミット**

```bash
cd C:/Users/yuu18/massbank-context && git add -A && git commit -m "feat: SQLite スキーマと索引ライタ"
```

---

### Task 8: 背景分布（IDF）の算出

spec §6.5。**df はレコード数ではなく unique InChIKey 数**で数える。
また照会時に ±1 ビンを探るため、**背景側も ±1 に展開した集合で df を数える**
（そうしないと照会の意味と df の意味がずれる）。

**Files:**
- Create: `massbank_context/index/background.py`
- Test: `tests/test_background.py`

**Interfaces:**
- Consumes: `open_index`（Task 7）
- Produces:
  - `build_background(conn: sqlite3.Connection) -> int` — 書き込んだ行数を返す
  - `lookup_idf(conn, acquisition_class: str, axis: str, bins: tuple[int, ...]) -> float`
    — 探索ビン群のうち最大の idf を返す（無ければ 0.0）

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_background.py`:

```python
import math
import sqlite3

import pytest

from massbank_context.index.background import build_background, lookup_idf
from massbank_context.index.writer import create_schema


@pytest.fixture()
def conn():
    connection = sqlite3.connect(":memory:")
    create_schema(connection)
    yield connection
    connection.close()


def _insert(conn, accession, inchikey, bins, axis="product_ion",
            acquisition_class="LC-ESI-MS2"):
    cursor = conn.execute(
        "INSERT INTO records (accession, inchikey, acquisition_class,"
        " resolution_class, polarity) VALUES (?,?,?,'high','POSITIVE')",
        (accession, inchikey, acquisition_class))
    for bin_no in bins:
        conn.execute(
            "INSERT INTO peaks (record_id, axis, bin, relative_intensity)"
            " VALUES (?,?,?,100.0)", (cursor.lastrowid, axis, bin_no))


def test_df_counts_unique_inchikeys_not_records(conn):
    """同じ化合物の 3 レコードは df=1（spec §6.5）。"""
    _insert(conn, "A1", "IK1", [1000])
    _insert(conn, "A2", "IK1", [1000])
    _insert(conn, "A3", "IK1", [1000])
    _insert(conn, "B1", "IK2", [2000])
    build_background(conn)
    row = conn.execute(
        "SELECT df_inchikey FROM axis_background WHERE bin=1000").fetchone()
    assert row[0] == 1


def test_idf_is_log_of_n_over_df(conn):
    _insert(conn, "A1", "IK1", [1000])
    _insert(conn, "B1", "IK2", [1000])
    _insert(conn, "C1", "IK3", [3000])
    _insert(conn, "D1", "IK4", [4000])
    build_background(conn)
    row = conn.execute(
        "SELECT idf FROM axis_background WHERE bin=1000").fetchone()
    assert abs(row[0] - math.log(4 / 2)) < 1e-9


def test_neighbouring_bins_are_counted(conn):
    """照会が ±1 を探るので、背景側も ±1 に展開して数える。"""
    _insert(conn, "A1", "IK1", [1000])
    build_background(conn)
    bins = {row[0] for row in conn.execute("SELECT bin FROM axis_background")}
    assert bins == {999, 1000, 1001}


def test_axes_are_separate(conn):
    _insert(conn, "A1", "IK1", [1000], axis="product_ion")
    _insert(conn, "B1", "IK2", [1000], axis="neutral_loss")
    build_background(conn)
    axes = {row[0] for row in conn.execute(
        "SELECT DISTINCT axis FROM axis_background")}
    assert axes == {"product_ion", "neutral_loss"}


def test_acquisition_classes_are_separate(conn):
    """測定原理の違うものを同じ背景で数えない（spec §6.1）。"""
    _insert(conn, "A1", "IK1", [1000], acquisition_class="LC-ESI-MS2")
    _insert(conn, "B1", "IK2", [1000], acquisition_class="GC-EI-MS")
    build_background(conn)
    classes = {row[0] for row in conn.execute(
        "SELECT DISTINCT acquisition_class FROM axis_background")}
    assert classes == {"LC-ESI-MS2", "GC-EI-MS"}


def test_lookup_idf_returns_max_of_probe_bins(conn):
    _insert(conn, "A1", "IK1", [1000])
    _insert(conn, "B1", "IK2", [1000])
    _insert(conn, "C1", "IK3", [5000])
    build_background(conn)
    value = lookup_idf(conn, "LC-ESI-MS2", "product_ion", (999, 1000, 1001))
    assert value > 0


def test_lookup_idf_returns_zero_for_unknown_bin(conn):
    _insert(conn, "A1", "IK1", [1000])
    build_background(conn)
    assert lookup_idf(conn, "LC-ESI-MS2", "product_ion", (90000,)) == 0.0
```

- [ ] **Step 2: テストが落ちることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_background.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'massbank_context.index.background'`

- [ ] **Step 3: 実装する**

`massbank_context/index/background.py`:

```python
"""背景頻度（IDF）の算出（spec §6.5）。

df をレコード数で数えると、人気化合物のレコードが装置・CE 違いで何十件も入っている
MassBank の収載バイアスがそのまま重みに乗る。df は必ず unique InChIKey で数える。

照会は probe_bins() で ±1 ビンを探るため、背景側も各ピークを ±1 に展開した集合で
df を数える。そうしないと「照会でヒットするのに df に数えられていない」ビンが生じ、
idf が過大になる。
"""
from __future__ import annotations

import math
import sqlite3
from collections import defaultdict


def build_background(conn: sqlite3.Connection) -> int:
    """axis_background を作り直し、書き込んだ行数を返す。"""
    conn.execute("DELETE FROM axis_background")

    # (acquisition_class, axis) ごとの母数 = その層に現れる unique InChIKey 数
    totals = {
        (acquisition_class, axis): n
        for acquisition_class, axis, n in conn.execute(
            "SELECT r.acquisition_class, p.axis, COUNT(DISTINCT r.inchikey)"
            " FROM peaks p JOIN records r ON r.record_id = p.record_id"
            " GROUP BY r.acquisition_class, p.axis")
    }

    # 各 (層, 軸, ビン±1) に現れる InChIKey の集合
    seen: dict[tuple[str, str, int], set[str]] = defaultdict(set)
    for acquisition_class, axis, bin_no, inchikey in conn.execute(
            "SELECT DISTINCT r.acquisition_class, p.axis, p.bin, r.inchikey"
            " FROM peaks p JOIN records r ON r.record_id = p.record_id"):
        for neighbour in (bin_no - 1, bin_no, bin_no + 1):
            seen[(acquisition_class, axis, neighbour)].add(inchikey)

    rows = []
    for (acquisition_class, axis, bin_no), inchikeys in seen.items():
        total = totals.get((acquisition_class, axis), 0)
        df = len(inchikeys)
        if total <= 0 or df <= 0:
            continue
        rows.append((acquisition_class, axis, bin_no, df,
                     math.log(total / df)))

    conn.executemany(
        "INSERT INTO axis_background"
        " (acquisition_class, axis, bin, df_inchikey, idf) VALUES (?,?,?,?,?)",
        rows)
    conn.commit()
    return len(rows)


def lookup_idf(conn: sqlite3.Connection, acquisition_class: str, axis: str,
               bins: tuple[int, ...]) -> float:
    """探索ビン群の idf のうち最大値を返す。該当が無ければ 0.0。

    最大を採るのは、隣接ビンのうち最も稀な（＝情報量の大きい）ビンで
    ヒットしうるため。合計にすると同じフラグメントを重複計上する。
    """
    if not bins:
        return 0.0
    placeholders = ",".join("?" for _ in bins)
    row = conn.execute(
        "SELECT MAX(idf) FROM axis_background"
        f" WHERE acquisition_class=? AND axis=? AND bin IN ({placeholders})",
        (acquisition_class, axis, *bins)).fetchone()
    return float(row[0]) if row and row[0] is not None else 0.0
```

- [ ] **Step 4: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_background.py -q`
Expected: PASS（7 件）

- [ ] **Step 5: コミット**

```bash
cd C:/Users/yuu18/massbank-context && git add -A && git commit -m "feat: 背景分布 IDF の算出（df は unique InChIKey）"
```

---

### Task 9: 索引ビルド CLI

ビン幅は全レコードを読まないと決まらないため **2 パス**にする。
パス 1 で TTL を全部読んで `SpectrumRecord` をメモリに持ち、質量誤差を測ってビン幅を決め、
パス 2 でビン化して書き込む。5 万件程度なら常識的なメモリに収まる。

**Files:**
- Create: `massbank_context/index/__main__.py`
- Test: `tests/test_index_cli.py`

**Interfaces:**
- Consumes: `read_ttl_file`（Task 4）、`measure_mass_error` / `choose_bin_width`（Task 6）、`IndexWriter`（Task 7）、`build_background`（Task 8）
- Produces:
  - `build_index(ttl_dir: str | Path, output: str | Path, resolution_filter: str = "high", limit: int | None = None) -> dict`
    — `{"scanned": int, "parsed": int, "accepted": int, "rejected": int, "bin_width": float, "background_rows": int}`
  - `main(argv: list[str] | None = None) -> int`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_index_cli.py`:

```python
import shutil
import sqlite3
from pathlib import Path

from massbank_context.index.__main__ import build_index, main
from massbank_context.index.writer import read_meta

FIXTURES = Path(__file__).parent / "fixtures" / "ttl"


def _ttl_dir(tmp_path: Path) -> Path:
    target = tmp_path / "ttl"
    target.mkdir()
    shutil.copy(FIXTURES / "gen_a" / "MSBNK-TEST-A0001.ttl", target)
    shutil.copy(FIXTURES / "gen_b" / "MSBNK-TEST-B0001.ttl", target)
    return target


def test_build_index_reads_both_generations(tmp_path):
    stats = build_index(_ttl_dir(tmp_path), tmp_path / "index.sqlite")
    assert stats["scanned"] == 2
    assert stats["parsed"] == 2
    assert stats["accepted"] == 2


def test_build_index_writes_meta(tmp_path):
    build_index(_ttl_dir(tmp_path), tmp_path / "index.sqlite")
    conn = sqlite3.connect(str(tmp_path / "index.sqlite"))
    meta = read_meta(conn)
    conn.close()
    assert meta["n_records"] == 2
    assert meta["resolution_filter"] == "high"
    assert meta["bin_width"] > 0
    assert meta["acquisition_classes"] == ["LC-ESI-MS2"]


def test_build_index_populates_background(tmp_path):
    stats = build_index(_ttl_dir(tmp_path), tmp_path / "index.sqlite")
    assert stats["background_rows"] > 0
    conn = sqlite3.connect(str(tmp_path / "index.sqlite"))
    count = conn.execute("SELECT COUNT(*) FROM axis_background").fetchone()[0]
    conn.close()
    assert count == stats["background_rows"]


def test_build_index_rejects_low_resolution(tmp_path):
    """低分解能レコードは既定で弾かれ、rejected に数えられる。"""
    ttl_dir = _ttl_dir(tmp_path)
    source = (FIXTURES / "gen_a" / "MSBNK-TEST-A0001.ttl").read_text(
        encoding="utf-8")
    (ttl_dir / "low.ttl").write_text(
        source.replace("MSBNK-TEST-A0001", "MSBNK-TEST-A0002")
              .replace("LC-ESI-ITFT", "LC-ESI-QQ"), encoding="utf-8")
    stats = build_index(ttl_dir, tmp_path / "index.sqlite")
    assert stats["scanned"] == 3
    assert stats["accepted"] == 2
    assert stats["rejected"] == 1


def test_build_index_is_idempotent(tmp_path):
    """同じ出力先に 2 回走らせても件数が倍にならない。"""
    ttl_dir = _ttl_dir(tmp_path)
    output = tmp_path / "index.sqlite"
    build_index(ttl_dir, output)
    build_index(ttl_dir, output)
    conn = sqlite3.connect(str(output))
    count = conn.execute("SELECT COUNT(*) FROM records").fetchone()[0]
    conn.close()
    assert count == 2


def test_main_returns_zero(tmp_path, capsys):
    code = main(["build", "--ttl-dir", str(_ttl_dir(tmp_path)),
                 "--output", str(tmp_path / "index.sqlite")])
    assert code == 0
    assert "accepted" in capsys.readouterr().out
```

- [ ] **Step 2: テストが落ちることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_index_cli.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'massbank_context.index.__main__'`

- [ ] **Step 3: 実装する**

`massbank_context/index/__main__.py`:

```python
"""索引ビルド CLI（spec §7.1）。

TTL 全走査は分〜時間を要する。MCP ツールにすると LLM が会話の途中で起動して
長時間固まるため、CLI に出す。MCP 側は index_status() で状態を読むだけ。

使い方:
    python -m massbank_context.index build --ttl-dir <dir> --output <sqlite>

ビン幅は全レコードの質量誤差を測らないと決まらないため 2 パス構成にする。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from massbank_context.core.serialization import json_payload
from massbank_context.index.background import build_background
from massbank_context.index.mass_error import choose_bin_width, measure_mass_error
from massbank_context.index.writer import IndexWriter, open_index
from massbank_context.ttl.reader import read_ttl_file


def build_index(ttl_dir: str | Path, output: str | Path,
                resolution_filter: str = "high",
                limit: int | None = None) -> dict:
    """TTL ディレクトリを走査して索引を作る。"""
    ttl_dir = Path(ttl_dir)
    paths = sorted(ttl_dir.rglob("*.ttl"))
    if limit is not None:
        paths = paths[:limit]

    # --- パス 1: 全部読む ---
    records = []
    failed = 0
    for path in paths:
        try:
            records.extend(read_ttl_file(path))
        except Exception:           # 1 件の異常でビルド全体を止めない
            failed += 1

    # --- ビン幅を実測から決める（spec §6.4）---
    profile = measure_mass_error(records)
    bin_width = choose_bin_width(profile)

    # --- パス 2: ビン化して書く ---
    conn = open_index(output, create=True)
    writer = IndexWriter(conn, bin_width=bin_width,
                         resolution_filter=resolution_filter)
    accepted = 0
    for record in records:
        if writer.add_record(record):
            accepted += 1
    writer.write_mass_error_profile(profile)
    writer.finalize({"ttl_source": str(ttl_dir)})
    background_rows = build_background(conn)
    conn.close()

    return {
        "scanned": len(paths),
        "parsed": len(records),
        "failed": failed,
        "accepted": accepted,
        "rejected": len(records) - accepted,
        "bin_width": bin_width,
        "background_rows": background_rows,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m massbank_context.index")
    sub = parser.add_subparsers(dest="command", required=True)

    build = sub.add_parser("build", help="TTL から SQLite 索引を作る")
    build.add_argument("--ttl-dir", required=True)
    build.add_argument("--output", required=True)
    build.add_argument("--resolution-filter", default="high",
                       help="'high' / 'low' / '' (無効化)")
    build.add_argument("--limit", type=int, default=None,
                       help="先頭 N ファイルだけ処理する（動作確認用）")

    args = parser.parse_args(argv)
    stats = build_index(args.ttl_dir, args.output,
                        resolution_filter=args.resolution_filter,
                        limit=args.limit)
    print(json_payload(stats))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_index_cli.py -q`
Expected: PASS（6 件）

- [ ] **Step 5: 全テストを走らせる**

Run: `C:/Python314/python.exe -m pytest tests -q`
Expected: PASS（ここまでの累計 72 件）

- [ ] **Step 6: コミット**

```bash
cd C:/Users/yuu18/massbank-context && git add -A && git commit -m "feat: 索引ビルド CLI（2 パス・ビン幅を実測から決定）"
```

---

### Task 10: クエリスペクトルのリーダ

**Files:**
- Create: `massbank_context/query/__init__.py`
- Create: `massbank_context/query/msp.py`
- Create: `massbank_context/query/dcl.py`
- Create: `massbank_context/query/loader.py`
- Create: `tests/fixtures/query/sample.msp`
- Test: `tests/test_query_loader.py`

`.dcl` リーダは Lipidmix `lipidmix/dcl/reader.py` の `deserialize_dcl()` を移植する。
返り値は `[{"dcl_index", "precursor_mz", "ion_mode", "rt", "msms_spectrum": [[mz, intensity], ...]}, ...]`。
`ion_mode` は int で、`0` が POSITIVE。

**Interfaces:**
- Consumes: `QuerySpectrum` / `Peak`（Task 2）
- Produces:
  - `read_msp(path: str | Path) -> list[QuerySpectrum]` — `.msp` / `.mgf` 両対応
  - `read_dcl(path: str | Path) -> list[QuerySpectrum]`
  - `load_spectra_file(path: str | Path) -> list[QuerySpectrum]` — 拡張子で振り分ける

- [ ] **Step 1: fixture を作る**

`tests/fixtures/query/sample.msp`:

```
NAME: TestCompound1
PRECURSORMZ: 230.0552
PRECURSORTYPE: [M+H]+
IONMODE: Positive
RETENTIONTIME: 0.876
Num Peaks: 3
86.0349 592827
171.0069 1271685
230.0552 5619710

NAME: TestCompound2
PRECURSORMZ: 267.0855
IONMODE: Negative
Num Peaks: 2
110.0349 250000
135.0434 999000
```

- [ ] **Step 2: 失敗するテストを書く**

`tests/test_query_loader.py`:

```python
from pathlib import Path

import pytest

from massbank_context.query.loader import load_spectra_file
from massbank_context.query.msp import read_msp

FIXTURE = Path(__file__).parent / "fixtures" / "query" / "sample.msp"


def test_reads_two_spectra():
    assert len(read_msp(FIXTURE)) == 2


def test_first_spectrum_fields():
    spectrum = read_msp(FIXTURE)[0]
    assert spectrum.precursor_mz == 230.0552
    assert spectrum.polarity == "POSITIVE"
    assert spectrum.retention_time == 0.876
    assert len(spectrum.peaks) == 3


def test_polarity_is_normalised_to_upper():
    """MSP の 'Positive' / 'Negative' を索引側の表記に合わせる。"""
    spectra = read_msp(FIXTURE)
    assert spectra[0].polarity == "POSITIVE"
    assert spectra[1].polarity == "NEGATIVE"


def test_relative_intensity_is_normalised_to_999():
    """装置間で絶対強度の桁が違うため、基準ピーク 999 に正規化する。"""
    spectrum = read_msp(FIXTURE)[0]
    assert max(p.relative_intensity for p in spectrum.peaks) == 999.0


def test_peaks_sorted_by_mz():
    spectrum = read_msp(FIXTURE)[0]
    mzs = [p.mz for p in spectrum.peaks]
    assert mzs == sorted(mzs)


def test_spectrum_ids_are_unique_and_stable():
    spectra = read_msp(FIXTURE)
    ids = [s.spectrum_id for s in spectra]
    assert ids == ["sample:1", "sample:2"]


def test_loader_dispatches_by_extension():
    assert len(load_spectra_file(FIXTURE)) == 2


def test_loader_rejects_unknown_extension(tmp_path):
    path = tmp_path / "x.csv"
    path.write_text("a,b\n", encoding="utf-8")
    with pytest.raises(ValueError, match="対応していない"):
        load_spectra_file(path)
```

- [ ] **Step 3: テストが落ちることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_query_loader.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'massbank_context.query'`

- [ ] **Step 4: 実装する**

`massbank_context/query/__init__.py` は空ファイル。

`massbank_context/query/msp.py`:

```python
"""MSP / MGF 形式の実験スペクトルを読む。

MassBank 側の relative_intensity は基準ピーク 999 のスケールなので、
クエリ側も同じスケールに正規化してから照合する。
"""
from __future__ import annotations

from pathlib import Path

from massbank_context.core.models import Peak, QuerySpectrum

BASE_PEAK_SCALE = 999.0

_FIELD_ALIASES = {
    "PRECURSORMZ": "precursor_mz",
    "PEPMASS": "precursor_mz",
    "IONMODE": "polarity",
    "POLARITY": "polarity",
    "RETENTIONTIME": "retention_time",
    "RTINSECONDS": "retention_time",
}


def _normalise_polarity(raw: str | None) -> str | None:
    if not raw:
        return None
    text = raw.strip().upper()
    if text.startswith("P") or text == "1" or text == "+":
        return "POSITIVE"
    if text.startswith("N") or text == "0" or text == "-":
        return "NEGATIVE"
    return None


def _to_peaks(pairs: list[tuple[float, float]]) -> tuple[Peak, ...]:
    if not pairs:
        return ()
    top = max(intensity for _mz, intensity in pairs) or 1.0
    peaks = [Peak(mz=mz,
                  relative_intensity=intensity / top * BASE_PEAK_SCALE,
                  formula=None, neutral_loss=None)
             for mz, intensity in pairs]
    peaks.sort(key=lambda p: p.mz)
    return tuple(peaks)


def read_msp(path: str | Path) -> list[QuerySpectrum]:
    """MSP / MGF を読む。空行または BEGIN/END IONS でレコードが区切られる。"""
    path = Path(path)
    stem = path.stem
    spectra: list[QuerySpectrum] = []
    fields: dict[str, str] = {}
    pairs: list[tuple[float, float]] = []

    def flush() -> None:
        if not pairs and not fields:
            return
        spectra.append(QuerySpectrum(
            spectrum_id=f"{stem}:{len(spectra) + 1}",
            precursor_mz=float(fields["precursor_mz"].split()[0])
            if "precursor_mz" in fields else None,
            polarity=_normalise_polarity(fields.get("polarity")),
            peaks=_to_peaks(pairs),
            source=str(path),
            retention_time=float(fields["retention_time"])
            if "retention_time" in fields else None,
        ))
        fields.clear()
        pairs.clear()

    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line or line.upper() == "END IONS":
            flush()
            continue
        if line.upper() == "BEGIN IONS":
            continue
        if ":" in line and not line[0].isdigit():
            key, _sep, value = line.partition(":")
            alias = _FIELD_ALIASES.get(key.strip().upper())
            if alias:
                fields[alias] = value.strip()
            continue
        parts = line.replace("\t", " ").split()
        if len(parts) >= 2:
            try:
                pairs.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    flush()
    return spectra
```

`massbank_context/query/dcl.py`:

```python
"""MS-DIAL .dcl（デコンボリューション済み MS/MS）を読む。

Lipidmix lipidmix/dcl/reader.py の deserialize_dcl() を移植したもの。
バイナリレイアウトは MS-DIAL の MSDecResult に対応する。推測で変えないこと。
"""
from __future__ import annotations

import struct
from pathlib import Path

from massbank_context.core.models import Peak, QuerySpectrum
from massbank_context.query.msp import BASE_PEAK_SCALE

DCL_MAGIC = b"DC"


def read_dcl(path: str | Path) -> list[QuerySpectrum]:
    path = Path(path)
    raw = path.read_bytes()
    if raw[:2] != DCL_MAGIC:
        raise ValueError(f"{path} は .dcl 形式ではありません（先頭が 'DC' でない）。")

    count = struct.unpack_from("<i", raw, 7)[0]
    seekpoints = struct.unpack_from(f"<{count}q", raw, 11)

    spectra: list[QuerySpectrum] = []
    for index, seekpoint in enumerate(seekpoints):
        offset = seekpoint
        precursor_mz = struct.unpack_from("<d", raw, offset + 16)[0]
        ion_mode = struct.unpack_from("<i", raw, offset + 24)[0]
        retention_time = struct.unpack_from("<d", raw, offset + 28)[0]
        offset += 60 + 40 + 20                     # Scan / Quant / Scoring
        spectra_n = struct.unpack_from("<3i", raw, offset)[0]
        offset += 12

        pairs: list[tuple[float, float]] = []
        for _ in range(spectra_n):
            mz, intensity = struct.unpack_from("<2d", raw, offset)
            offset += 20                            # double + double + int32
            pairs.append((mz, intensity))

        top = max((intensity for _mz, intensity in pairs), default=0.0) or 1.0
        peaks = tuple(sorted(
            (Peak(mz=mz, relative_intensity=intensity / top * BASE_PEAK_SCALE,
                  formula=None, neutral_loss=None)
             for mz, intensity in pairs),
            key=lambda p: p.mz))

        spectra.append(QuerySpectrum(
            spectrum_id=f"{path.stem}:{index}",
            precursor_mz=precursor_mz,
            polarity="POSITIVE" if ion_mode == 0 else "NEGATIVE",
            peaks=peaks,
            source=str(path),
            retention_time=retention_time,
        ))
    return spectra
```

`massbank_context/query/loader.py`:

```python
"""実験スペクトルの読み込み入口。拡張子で振り分ける。"""
from __future__ import annotations

from pathlib import Path

from massbank_context.core.models import QuerySpectrum
from massbank_context.query.dcl import read_dcl
from massbank_context.query.msp import read_msp

SUPPORTED_SUFFIXES = (".msp", ".mgf", ".dcl")


def load_spectra_file(path: str | Path) -> list[QuerySpectrum]:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in (".msp", ".mgf"):
        return read_msp(path)
    if suffix == ".dcl":
        return read_dcl(path)
    raise ValueError(
        f"対応していない拡張子です: {suffix}。"
        f"対応形式: {', '.join(SUPPORTED_SUFFIXES)}")
```

- [ ] **Step 5: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_query_loader.py -q`
Expected: PASS（8 件）

- [ ] **Step 6: コミット**

```bash
cd C:/Users/yuu18/massbank-context && git add -A && git commit -m "feat: MSP/MGF/.dcl のクエリスペクトルリーダ"
```

---

### Task 11: 照合エンジン

spec §6.8。**IDF 重み付きカバレッジを主スコアとし、cosine は併記に留める。**
`band` は一致したプロダクトイオンの本数だけで決まる。

**Files:**
- Create: `massbank_context/spectral/__init__.py`
- Create: `massbank_context/spectral/matcher.py`
- Test: `tests/test_matcher.py`

**Interfaces:**
- Consumes: `QuerySpectrum`（Task 2）、`probe_bins` / `to_bin`（Task 5）、`lookup_idf`（Task 8）、`read_meta`（Task 7）
- Produces:
  - `MatchCandidate` dataclass — `accession, inchikey, matched_fragments, total_query_peaks, weighted_score, cosine, unexplained_peak_ratio, band, precursor_delta`
  - `match_spectrum(conn, spectrum: QuerySpectrum, acquisition_class: str, precursor_tolerance: float = 0.01, top_k: int = 10) -> list[MatchCandidate]`
  - `derive_band(matched_product_ions: int) -> str`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_matcher.py`:

```python
import shutil
import sqlite3
from pathlib import Path

import pytest

from massbank_context.core.models import Peak, QuerySpectrum
from massbank_context.index.__main__ import build_index
from massbank_context.spectral.matcher import derive_band, match_spectrum

FIXTURES = Path(__file__).parent / "fixtures" / "ttl"


@pytest.fixture()
def conn(tmp_path):
    ttl_dir = tmp_path / "ttl"
    ttl_dir.mkdir()
    shutil.copy(FIXTURES / "gen_a" / "MSBNK-TEST-A0001.ttl", ttl_dir)
    build_index(ttl_dir, tmp_path / "index.sqlite")
    connection = sqlite3.connect(str(tmp_path / "index.sqlite"))
    yield connection
    connection.close()


def _query(peaks, precursor_mz=230.0552) -> QuerySpectrum:
    return QuerySpectrum(
        spectrum_id="q1", precursor_mz=precursor_mz, polarity="POSITIVE",
        peaks=tuple(Peak(mz, intensity, None, None) for mz, intensity in peaks),
        source="test")


@pytest.mark.parametrize("matched,expected", [
    (0, "PRECURSOR_ONLY"),
    (1, "WEAK"),
    (2, "PASS"),
    (7, "PASS"),
])
def test_derive_band(matched, expected):
    assert derive_band(matched) == expected


def test_exact_match_scores_high(conn):
    query = _query([(86.0349, 11.0), (171.0069, 226.0), (230.0552, 999.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert results[0].accession == "MSBNK-TEST-A0001"
    assert results[0].matched_fragments == 3
    assert results[0].band == "PASS"
    assert results[0].unexplained_peak_ratio == 0.0


def test_precursor_only_match_is_banded(conn):
    """precursor だけ合ってプロダクトイオンが合わない候補は PRECURSOR_ONLY。"""
    query = _query([(999.9999, 100.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert results[0].band == "PRECURSOR_ONLY"
    assert results[0].matched_fragments == 0
    assert results[0].unexplained_peak_ratio == 1.0


def test_precursor_mismatch_returns_nothing(conn):
    query = _query([(86.0349, 11.0)], precursor_mz=500.0)
    assert match_spectrum(conn, query, "LC-ESI-MS2") == []


def test_polarity_mismatch_returns_nothing(conn):
    query = QuerySpectrum(spectrum_id="q1", precursor_mz=230.0552,
                          polarity="NEGATIVE",
                          peaks=(Peak(86.0349, 11.0, None, None),),
                          source="test")
    assert match_spectrum(conn, query, "LC-ESI-MS2") == []


def test_boundary_peak_still_matches(conn):
    """ビン境界をまたぐ僅かなずれでも ±1 探索で拾える（Task 5）。"""
    query = _query([(86.0350, 11.0), (171.0070, 226.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert results[0].matched_fragments == 2


def test_unexplained_ratio_counts_query_peaks(conn):
    query = _query([(86.0349, 11.0), (400.0000, 50.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert results[0].matched_fragments == 1
    assert abs(results[0].unexplained_peak_ratio - 0.5) < 1e-9


def test_weighted_score_is_between_zero_and_one(conn):
    query = _query([(86.0349, 11.0), (400.0000, 50.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert 0.0 <= results[0].weighted_score <= 1.0


def test_top_k_limits_results(conn):
    query = _query([(86.0349, 11.0)])
    assert len(match_spectrum(conn, query, "LC-ESI-MS2", top_k=1)) <= 1


def test_results_sorted_by_weighted_score_only(conn):
    """並び順を決めるのは spectral_evidence だけ（spec §5.3 性質 1）。"""
    query = _query([(86.0349, 11.0), (171.0069, 226.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    scores = [r.weighted_score for r in results]
    assert scores == sorted(scores, reverse=True)
```

- [ ] **Step 2: テストが落ちることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_matcher.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'massbank_context.spectral'`

- [ ] **Step 3: 実装する**

`massbank_context/spectral/__init__.py` は空ファイル。

`massbank_context/spectral/matcher.py`:

```python
"""スペクトル照合（spec §6.8）。

主スコアは IDF 重み付きカバレッジ:
    weighted_score = Σ idf(f) for f in matched / Σ idf(f) for f in query_peaks

cosine ではなくカバレッジを主にする理由:
  1. unexplained_peak_ratio と直接対になり、「クエリのピークのうち何割が説明され、
     残りは何か」が一続きの説明になる。cosine の 0.81 は説明できない。
  2. cosine は強度比に敏感で、装置間・CE 間の差を拾う。カバレッジは頑健。

cosine も併記するが、並び順には使わない。
"""
from __future__ import annotations

import math
import sqlite3
from dataclasses import dataclass

from massbank_context.core.binning import PRODUCT_ION, probe_bins, to_bin
from massbank_context.core.models import QuerySpectrum
from massbank_context.index.background import lookup_idf
from massbank_context.index.writer import read_meta

BAND_PASS = "PASS"
BAND_WEAK = "WEAK"
BAND_PRECURSOR_ONLY = "PRECURSOR_ONLY"

DEFAULT_PRECURSOR_TOLERANCE = 0.01


@dataclass(frozen=True)
class MatchCandidate:
    """1 候補の spectral_evidence（spec §5.3）。

    メタデータは一切含まない。並び順を決めるのはこのブロックだけである。
    """
    accession: str
    inchikey: str
    matched_fragments: int
    total_query_peaks: int
    weighted_score: float
    cosine: float
    unexplained_peak_ratio: float
    band: str
    precursor_delta: float | None


def derive_band(matched_product_ions: int) -> str:
    """一致したプロダクトイオンの本数から band を決める（spec §5.3）。

    PRECURSOR_ONLY と PASS を同一視してはならない。精密質量が合っているだけの
    候補を、フラグメントが合った候補と同じ顔で並べないための区別である。
    """
    if matched_product_ions >= 2:
        return BAND_PASS
    if matched_product_ions == 1:
        return BAND_WEAK
    return BAND_PRECURSOR_ONLY


def _cosine(query_intensities: dict[int, float],
            record_intensities: dict[int, float]) -> float:
    shared = set(query_intensities) & set(record_intensities)
    if not shared:
        return 0.0
    dot = sum(query_intensities[b] * record_intensities[b] for b in shared)
    norm_q = math.sqrt(sum(v * v for v in query_intensities.values()))
    norm_r = math.sqrt(sum(v * v for v in record_intensities.values()))
    if norm_q == 0 or norm_r == 0:
        return 0.0
    return dot / (norm_q * norm_r)


def match_spectrum(conn: sqlite3.Connection,
                   spectrum: QuerySpectrum,
                   acquisition_class: str,
                   precursor_tolerance: float = DEFAULT_PRECURSOR_TOLERANCE,
                   top_k: int = 10) -> list[MatchCandidate]:
    """1 本のクエリスペクトルに対する候補を返す（weighted_score の降順）。"""
    meta = read_meta(conn)
    if meta is None:
        return []
    bin_width = meta["bin_width"]

    # --- 手順 2: precursor m/z と polarity で候補レコードを絞る ---
    if spectrum.precursor_mz is None:
        rows = conn.execute(
            "SELECT record_id, accession, inchikey, precursor_mz FROM records"
            " WHERE acquisition_class=? AND (polarity=? OR ? IS NULL)",
            (acquisition_class, spectrum.polarity, spectrum.polarity)).fetchall()
    else:
        rows = conn.execute(
            "SELECT record_id, accession, inchikey, precursor_mz FROM records"
            " WHERE acquisition_class=? AND (polarity=? OR ? IS NULL)"
            "   AND precursor_mz BETWEEN ? AND ?",
            (acquisition_class, spectrum.polarity, spectrum.polarity,
             spectrum.precursor_mz - precursor_tolerance,
             spectrum.precursor_mz + precursor_tolerance)).fetchall()
    if not rows:
        return []

    # --- クエリ側のビンと idf を先に用意する ---
    query_bins = [probe_bins(peak.mz, bin_width) for peak in spectrum.peaks]
    query_idf = [lookup_idf(conn, acquisition_class, PRODUCT_ION, bins)
                 for bins in query_bins]
    query_intensities = {to_bin(peak.mz, bin_width): peak.relative_intensity
                         for peak in spectrum.peaks}
    total_idf = sum(query_idf)
    total_peaks = len(spectrum.peaks)

    candidates: list[MatchCandidate] = []
    for record_id, accession, inchikey, record_precursor in rows:
        record_bins = {
            bin_no: intensity
            for bin_no, intensity in conn.execute(
                "SELECT bin, relative_intensity FROM peaks"
                " WHERE record_id=? AND axis=?", (record_id, PRODUCT_ION))
        }
        matched_idf = 0.0
        matched = 0
        for bins, idf in zip(query_bins, query_idf):
            if any(b in record_bins for b in bins):
                matched += 1
                matched_idf += idf

        candidates.append(MatchCandidate(
            accession=accession,
            inchikey=inchikey,
            matched_fragments=matched,
            total_query_peaks=total_peaks,
            weighted_score=(matched_idf / total_idf) if total_idf > 0 else 0.0,
            cosine=_cosine(query_intensities, record_bins),
            unexplained_peak_ratio=(
                1.0 - matched / total_peaks) if total_peaks else 1.0,
            band=derive_band(matched),
            precursor_delta=(
                spectrum.precursor_mz - record_precursor)
            if spectrum.precursor_mz is not None and record_precursor is not None
            else None,
        ))

    # 並び順を決めるのは weighted_score だけ。同点は accession で決定的にする。
    candidates.sort(key=lambda c: (-c.weighted_score, c.accession))
    return candidates[:top_k]
```

- [ ] **Step 4: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_matcher.py -q`
Expected: PASS（13 件）

- [ ] **Step 5: コミット**

```bash
cd C:/Users/yuu18/massbank-context && git add -A && git commit -m "feat: 照合エンジン（IDF 重み付きカバレッジ・band 判定）"
```

---

### Task 12: MCP サーバとツール 4 本

**Files:**
- Create: `massbank_context/core/session.py`
- Create: `massbank_context/core/config.py`
- Create: `massbank_context/tools/__init__.py`
- Create: `massbank_context/tools/index_tools.py`
- Create: `massbank_context/tools/spectra_tools.py`
- Create: `massbank_context/tools/match_tools.py`
- Create: `massbank_context/server.py`
- Create: `README.md`
- Test: `tests/test_tools.py`

**Interfaces:**
- Consumes: 全タスクの成果物
- Produces:
  - `index_status() -> str` / `load_spectra(path: str) -> str` / `inspect_spectrum(spectrum_id: str) -> str` / `match_spectra(...) -> str`
  - `session` — `massbank_context.core.session` のモジュール変数。**可変状態の正準はここ**

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_tools.py`:

```python
import json
import shutil
from pathlib import Path

import pytest

from massbank_context.core import config, session as session_module
from massbank_context.index.__main__ import build_index
from massbank_context.tools.index_tools import index_status
from massbank_context.tools.match_tools import match_spectra
from massbank_context.tools.spectra_tools import inspect_spectrum, load_spectra

TTL = Path(__file__).parent / "fixtures" / "ttl" / "gen_a" / "MSBNK-TEST-A0001.ttl"
MSP = Path(__file__).parent / "fixtures" / "query" / "sample.msp"


@pytest.fixture(autouse=True)
def clean_session(tmp_path, monkeypatch):
    """索引パスとセッションをテストごとに作り直す。実ファイルに依存させない。"""
    monkeypatch.setattr(config, "INDEX_PATH", tmp_path / "index.sqlite")
    session_module.session.reset()
    yield


def _build(tmp_path):
    ttl_dir = tmp_path / "ttl"
    ttl_dir.mkdir()
    shutil.copy(TTL, ttl_dir)
    build_index(ttl_dir, config.INDEX_PATH)


def test_index_status_without_index():
    payload = json.loads(index_status())
    assert payload["error"]["code"] == "missing_index"
    assert "index build" in payload["error"]["required_action"]


def test_index_status_with_index(tmp_path):
    _build(tmp_path)
    payload = json.loads(index_status())
    assert payload["n_records"] == 1
    assert payload["acquisition_classes"] == ["LC-ESI-MS2"]
    assert payload["bin_width"] > 0


def test_load_spectra_returns_summary_not_peaks():
    payload = json.loads(load_spectra(str(MSP)))
    assert payload["n_spectra"] == 2
    assert payload["polarity_counts"] == {"POSITIVE": 1, "NEGATIVE": 1}
    assert "peaks" not in json.dumps(payload)


def test_inspect_spectrum_requires_load():
    payload = json.loads(inspect_spectrum("sample:1"))
    assert payload["error"]["code"] == "missing_state"
    assert payload["error"]["required_tools"] == ["load_spectra"]


def test_inspect_spectrum_returns_rounded_peaks():
    load_spectra(str(MSP))
    payload = json.loads(inspect_spectrum("sample:1"))
    assert payload["spectrum_id"] == "sample:1"
    assert len(payload["peaks"]) == 3
    for mz, _intensity in payload["peaks"]:
        assert len(str(mz).split(".")[-1]) <= 4


def test_match_spectra_requires_index():
    load_spectra(str(MSP))
    payload = json.loads(match_spectra())
    assert payload["error"]["code"] == "missing_index"


def test_match_spectra_requires_spectra(tmp_path):
    _build(tmp_path)
    payload = json.loads(match_spectra())
    assert payload["error"]["code"] == "missing_state"


def test_match_spectra_rejects_incompatible_class(tmp_path):
    """測定原理の異なる層を跨いだ照合は拒否する（spec §6.1）。"""
    _build(tmp_path)
    load_spectra(str(MSP))
    payload = json.loads(match_spectra(acquisition_class="GC-EI-MS"))
    assert payload["error"]["code"] == "incompatible_index"
    assert payload["error"]["available"] == ["LC-ESI-MS2"]


def test_match_spectra_returns_tsv(tmp_path):
    _build(tmp_path)
    load_spectra(str(MSP))
    payload = json.loads(match_spectra())
    header = payload["candidates"].splitlines()[0].split("\t")
    assert header == ["spectrum_id", "accession", "inchikey",
                      "matched_fragments", "total_query_peaks",
                      "weighted_score", "cosine", "unexplained_peak_ratio",
                      "band", "precursor_delta"]


def test_match_spectra_reports_bias_not_applied(tmp_path):
    """objective 未設定はエラーにせず bias.applied=false で返す（spec §8）。"""
    _build(tmp_path)
    load_spectra(str(MSP))
    payload = json.loads(match_spectra())
    assert payload["bias"]["applied"] is False
    assert payload["bias"]["reason"]
```

- [ ] **Step 2: テストが落ちることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_tools.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'massbank_context.core.config'`

- [ ] **Step 3: 設定とセッションを書く**

`massbank_context/core/config.py`:

```python
"""パス解決と環境変数（すべて任意・上書き用）。

可変状態の正準はこのモジュール。参照側は `config.INDEX_PATH` のように
module 修飾で読むこと。`from ... import INDEX_PATH` はスナップショット束縛に
なり、テストでの差し替えが効かなくなる。
"""
from __future__ import annotations

import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent.parent

INDEX_PATH = Path(os.environ.get(
    "MASSBANK_CONTEXT_INDEX", BASE_DIR / "index" / "massbank.sqlite"))

BUILD_COMMAND = ("CLI: C:/Python314/python.exe -m massbank_context.index build"
                 " --ttl-dir <dir> --output <sqlite>")
```

`massbank_context/core/session.py`:

```python
"""セッション状態。可変状態の正準はこのモジュールの `session`。"""
from __future__ import annotations

from massbank_context.core.models import QuerySpectrum


class Session:
    """読み込んだクエリスペクトルを保持する。

    ピーク列はここに置き、ツールの戻り値には載せない（spec §12）。
    """

    def __init__(self) -> None:
        self.spectra: dict[str, QuerySpectrum] = {}
        self.source: str | None = None

    def reset(self) -> None:
        self.spectra = {}
        self.source = None

    def set_spectra(self, spectra: list[QuerySpectrum], source: str) -> None:
        self.spectra = {s.spectrum_id: s for s in spectra}
        self.source = source


session = Session()
```

- [ ] **Step 4: ツールを書く**

`massbank_context/tools/__init__.py` は空ファイル。

`massbank_context/tools/index_tools.py`:

```python
"""索引の状態を返すツール。"""
from __future__ import annotations

import sqlite3

from massbank_context.core import config
from massbank_context.core.errors import missing_index
from massbank_context.core.serialization import json_payload, round_floats
from massbank_context.index.writer import read_meta

__all__ = ["index_status"]


def index_status() -> str:
    """索引の有無と、その層・件数・ビン幅を返す。"""
    if not config.INDEX_PATH.exists():
        return missing_index(config.BUILD_COMMAND)
    conn = sqlite3.connect(str(config.INDEX_PATH))
    try:
        meta = read_meta(conn)
        if meta is None:
            return missing_index(config.BUILD_COMMAND)
        polarities = sorted(
            row[0] for row in conn.execute(
                "SELECT DISTINCT polarity FROM records WHERE polarity IS NOT NULL"))
    finally:
        conn.close()
    return json_payload(round_floats({**meta, "polarities": polarities}))
```

`massbank_context/tools/spectra_tools.py`:

```python
"""実験スペクトルの読み込みと確認。"""
from __future__ import annotations

from collections import Counter

from massbank_context.core.errors import missing_state
from massbank_context.core.serialization import MZ_DIGITS, json_payload, round_floats
from massbank_context.core.session import session
from massbank_context.query.loader import load_spectra_file

__all__ = ["load_spectra", "inspect_spectrum"]


def load_spectra(path: str) -> str:
    """スペクトルを読み込み、要約だけを返す。ピーク列はセッションに残す。"""
    try:
        spectra = load_spectra_file(path)
    except (OSError, ValueError) as exc:
        return json_payload({"error": {"code": "load_failed",
                                       "message": str(exc)}})
    session.set_spectra(spectra, path)

    precursors = [s.precursor_mz for s in spectra if s.precursor_mz is not None]
    return json_payload(round_floats({
        "source": path,
        "n_spectra": len(spectra),
        "polarity_counts": dict(Counter(
            s.polarity for s in spectra if s.polarity)),
        "precursor_mz_range": [min(precursors), max(precursors)]
        if precursors else None,
        "spectrum_ids": [s.spectrum_id for s in spectra],
    }, MZ_DIGITS))


def inspect_spectrum(spectrum_id: str) -> str:
    """1 本のピーク列を丸めて返す。必要なときだけ呼ぶこと。"""
    if not session.spectra:
        return missing_state(
            "spectra", ["load_spectra"],
            "実験スペクトルが読み込まれていません。先に load_spectra を実行してください。")
    spectrum = session.spectra.get(spectrum_id)
    if spectrum is None:
        return json_payload({"error": {
            "code": "unknown_spectrum_id",
            "message": f"{spectrum_id} は読み込まれていません。",
            "available": sorted(session.spectra)[:20],
        }})
    return json_payload(round_floats({
        "spectrum_id": spectrum.spectrum_id,
        "precursor_mz": spectrum.precursor_mz,
        "polarity": spectrum.polarity,
        "retention_time": spectrum.retention_time,
        "peaks": [[p.mz, p.relative_intensity] for p in spectrum.peaks],
    }, MZ_DIGITS))
```

`massbank_context/tools/match_tools.py`:

```python
"""照合ツール。Phase 1 ではバイアスを適用しない。"""
from __future__ import annotations

import json
import sqlite3

from massbank_context.core import config
from massbank_context.core.errors import (
    incompatible_index,
    missing_index,
    missing_state,
)
from massbank_context.core.serialization import as_tsv, json_payload
from massbank_context.core.session import session
from massbank_context.index.writer import read_meta
from massbank_context.spectral.matcher import (
    DEFAULT_PRECURSOR_TOLERANCE,
    match_spectrum,
)

__all__ = ["match_spectra"]

_COLUMNS = ["spectrum_id", "accession", "inchikey", "matched_fragments",
            "total_query_peaks", "weighted_score", "cosine",
            "unexplained_peak_ratio", "band", "precursor_delta"]

_DEFAULT_ACQUISITION_CLASS = "LC-ESI-MS2"


def match_spectra(spectrum_ids: list[str] | None = None,
                  acquisition_class: str = _DEFAULT_ACQUISITION_CLASS,
                  precursor_tolerance: float = DEFAULT_PRECURSOR_TOLERANCE,
                  top_k: int = 10) -> str:
    """読み込んだスペクトルを索引と照合し、候補を TSV で返す。

    Phase 1 ではバイアス層が無いため bias.applied は常に false。
    objective 未設定はエラーにしない（spec §8）。
    """
    if not config.INDEX_PATH.exists():
        return missing_index(config.BUILD_COMMAND)
    if not session.spectra:
        return missing_state(
            "spectra", ["load_spectra"],
            "実験スペクトルが読み込まれていません。先に load_spectra を実行してください。")

    conn = sqlite3.connect(str(config.INDEX_PATH))
    try:
        meta = read_meta(conn)
        if meta is None:
            return missing_index(config.BUILD_COMMAND)
        available = meta.get("acquisition_classes", [])
        if acquisition_class not in available:
            return incompatible_index(
                acquisition_class, available,
                f"索引には {available} の層しかありません。測定原理の異なる層を"
                f"跨いだ照合は結果が無意味になるため行いません。")

        targets = (spectrum_ids if spectrum_ids is not None
                   else sorted(session.spectra))
        rows = []
        for spectrum_id in targets:
            spectrum = session.spectra.get(spectrum_id)
            if spectrum is None:
                continue
            for candidate in match_spectrum(
                    conn, spectrum, acquisition_class,
                    precursor_tolerance=precursor_tolerance, top_k=top_k):
                rows.append([
                    spectrum_id, candidate.accession, candidate.inchikey,
                    candidate.matched_fragments, candidate.total_query_peaks,
                    round(candidate.weighted_score, 4),
                    round(candidate.cosine, 4),
                    round(candidate.unexplained_peak_ratio, 4),
                    candidate.band,
                    None if candidate.precursor_delta is None
                    else round(candidate.precursor_delta, 4),
                ])
    finally:
        conn.close()

    return json_payload({
        "acquisition_class": acquisition_class,
        "n_spectra": len(targets),
        "n_candidates": len(rows),
        "candidates": as_tsv(_COLUMNS, rows),
        "bias": {
            "applied": False,
            "reason": "Phase 1 ではバイアス層が未実装。並び順は weighted_score のみで決まる。",
        },
        "note": "band が PRECURSOR_ONLY の候補は精密質量しか合っていない。"
                "PASS と同一視しないこと。",
    })
```

- [ ] **Step 5: サーバのファサードと README を書く**

`massbank_context/server.py`:

```python
"""MCP サーバのファサード。

import 副作用で @mcp.tool を登録し、各モジュールの __all__ を再エクスポートする。
新しいツールを足すときは、実体を該当モジュールに書き __all__ に載せること。
"""
from __future__ import annotations

from mcp.server.fastmcp import FastMCP

from massbank_context.tools.index_tools import index_status as _index_status
from massbank_context.tools.match_tools import match_spectra as _match_spectra
from massbank_context.tools.spectra_tools import (
    inspect_spectrum as _inspect_spectrum,
    load_spectra as _load_spectra,
)

mcp = FastMCP("massbank-context")

# structured_output=False は必須。付けないと FastMCP が outputSchema を導出し、
# MCP が同じ内容を content と structuredContent の両方で送る（＝2 倍）。
mcp.tool(structured_output=False)(_index_status)
mcp.tool(structured_output=False)(_load_spectra)
mcp.tool(structured_output=False)(_inspect_spectrum)
mcp.tool(structured_output=False)(_match_spectra)


if __name__ == "__main__":
    mcp.run()
```

`README.md`:

```markdown
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
```

- [ ] **Step 6: テストが通ることを確認する**

Run: `C:/Python314/python.exe -m pytest tests/test_tools.py -q`
Expected: PASS（10 件）

- [ ] **Step 7: 全テストを走らせる**

Run: `C:/Python314/python.exe -m pytest tests -q`
Expected: PASS（累計 103 件）

- [ ] **Step 8: サーバが起動することを確認する**

Run: `C:/Python314/python.exe -c "import massbank_context.server as s; print(len(s.mcp._tool_manager.list_tools()))"`
Expected: `4`

- [ ] **Step 9: コミット**

```bash
cd C:/Users/yuu18/massbank-context && git add -A && git commit -m "feat: MCP サーバとツール 4 本（index_status / load_spectra / inspect_spectrum / match_spectra）"
```

---

## Phase 1 完了後の受け入れ確認

実 TTL での動作確認。fixture ではなく本物で 1 回通す。

- [ ] **MassBank-RDF を取得する**

```bash
git clone --depth 1 https://github.com/systemsomicslab/MassBank-RDF.git C:/Users/yuu18/MassBank-RDF
```

- [ ] **小さく試す（先頭 500 ファイル）**

```bash
C:/Python314/python.exe -m massbank_context.index build --ttl-dir C:/Users/yuu18/MassBank-RDF/draft_ttl --output C:/Users/yuu18/massbank-context/index/sample.sqlite --limit 500
```

期待: `accepted` が 0 でないこと。`bin_width` が 0.001〜0.05 の範囲に入ること
（この範囲を外れたら `measure_mass_error` のクラスタリングを疑う）。

- [ ] **全件でビルドする**

```bash
C:/Python314/python.exe -m massbank_context.index build --ttl-dir C:/Users/yuu18/MassBank-RDF/draft_ttl --output C:/Users/yuu18/massbank-context/index/massbank.sqlite
```

期待: `scanned` が 54,183。`acquisition_classes` に `LC-ESI-MS2` が入ること。
`Fac_Eng_Univ_Tokyo` の 12,378 件は `GC-EI-MS` かつ `resolution_class` が
`unknown` なので既定では弾かれ、`rejected` に数えられる。

- [ ] **既知化合物で逆引きが当たることを確認する**

`draft_ttl` から任意の LC-ESI-MS2 レコードを 1 件選び、そのピーク列を `.msp` に書き出して
`load_spectra` → `match_spectra` を通し、**そのレコード自身が 1 位に来ること**を確認する。
来なければ照合エンジンかビン幅を疑う。

---

## Self-Review

**Spec coverage:**

| spec の節 | 実装するタスク |
|---|---|
| §3.3 世代A/B の食い違い | Task 3（世代判定・世代A）/ Task 4（世代B・振り分け） |
| §3.4 世代A のフィールド | Task 3 |
| §5.1 モジュール境界 | 全タスク（File Structure の通り） |
| §5.3 証拠バンドの契約・band 定義 | Task 11（`MatchCandidate` / `derive_band`）/ Task 12（TSV の列） |
| §6.1 acquisition_class での分割 | Task 1（導出）/ Task 7（索引に格納）/ Task 8（背景を層別）/ Task 12（跨ぎを拒否） |
| §6.2 resolution_class の導出 | Task 1 / Task 7（既定 high のみ受け入れ） |
| §6.3 3 軸（うち product_ion / neutral_loss） | Task 5 |
| §6.4 ビン幅を実測から決める | Task 6 / Task 9（2 パス構成） |
| §6.5 IDF の df は unique InChIKey | Task 8 |
| §6.6 濃縮度の背景 | **Phase 3 の範囲**。Task 7 で `compound_class` テーブルだけ先に埋める |
| §6.7 SQLite スキーマ | Task 7 |
| §6.8 IDF 重み付きカバレッジ・band | Task 11 |
| §7.1 索引ビルドは CLI | Task 9 |
| §7.2 ツール 1・2・3・5 | Task 12 |
| §7.3 `set_objective` | **Phase 3 の範囲**。Phase 1 では扱わない |
| §7.4 想起拡張 | **Phase 3 の範囲** |
| §8 エラー封筒（4 種） | Task 2（封筒）/ Task 12（適用） |
| §10 テスト戦略 | 各タスクのテスト。実 TTL 同梱は Task 3、合成 gen-B は Task 4 |
| §12 payload 規律 | Task 2（`json_payload` / `round_floats` / `as_tsv`）/ Task 12（`structured_output=False`） |

Phase 1 のスコープ外は §6.6 の統計・§7.3・§7.4 のみで、いずれも spec §9 で
Phase 3 に割り当てられている。

**`formula` 軸について:** spec §6.3 は 3 軸を挙げるが、`formula` は一次キーではなく
「存在するときだけ付ける注釈」と定義されている。Task 5 / Task 7 で `peaks.formula` に
格納し、照合には使わない。これは spec の定義通りで、欠落ではない。

**Placeholder scan:** TBD・省略・「Task N と同様」は無い。各ステップは実際に
貼り付けて動くコードを持つ。

**Type consistency:**

- `to_bin(mz, bin_width) -> int` / `probe_bins(mz, bin_width) -> tuple[int,int,int]` —
  Task 5 で定義、Task 11 で使用。引数順一致
- `lookup_idf(conn, acquisition_class, axis, bins)` — Task 8 で定義、Task 11 で使用
- `read_meta(conn) -> dict | None` — Task 7 で定義、Task 11・Task 12 で使用。
  `bin_width` は float、`acquisition_classes` は list[str] に復元される
- `IndexWriter.add_record(record) -> bool` — Task 7 で定義、Task 9 で使用
- `derive_acquisition_class` / `derive_resolution_class` — Task 1 で定義、
  Task 2（`SpectrumRecord` プロパティ）・Task 7 で使用
- `BASE_PEAK_SCALE` — Task 10 の `msp.py` で定義し `dcl.py` が import。
  MassBank 側の `relative_intensity`（基準 999）とスケールが揃う
- `PRODUCT_ION` / `NEUTRAL_LOSS` — Task 5 で定義、Task 8・Task 11 で文字列直書きせず参照

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-08-28-massbank-context-phase1.md`.
