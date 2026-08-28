import json
import shutil
import struct
from pathlib import Path

import pytest

from massbank_context.core import config, session as session_module
from massbank_context.index.__main__ import build_index
from massbank_context.tools import match_tools
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


def test_load_spectra_returns_envelope_for_truncated_dcl(tmp_path):
    """壊れた .dcl はトレースバックではなく load_failed の封筒で返す（spec §8）。"""
    # ヘッダは完全だが本体が無い .dcl（書き込み中断・コピー中断で生じる）。
    # マジック 'DC' + version + has_annotation + レコード数 2 + seekpoint 2 個。
    header = (b"DC" + struct.pack("<i", 1) + bytes(1) + struct.pack("<i", 2)
              + struct.pack("<2q", 27, 200))
    path = tmp_path / "truncated.dcl"
    path.write_bytes(header)
    payload = json.loads(load_spectra(str(path)))
    assert payload["error"]["code"] == "load_failed"
    assert "truncated.dcl" in payload["error"]["message"]


def _write_non_sqlite_index():
    """索引パスに SQLite でないファイルを置く（ビルド中断・コピー中断の再現）。"""
    config.INDEX_PATH.parent.mkdir(parents=True, exist_ok=True)
    config.INDEX_PATH.write_bytes(b"not a sqlite database at all" * 8)


def test_index_status_reports_missing_index_for_corrupt_file():
    """壊れた索引は使えない。CLI でしか直せないので missing_index を返す。"""
    _write_non_sqlite_index()
    payload = json.loads(index_status())
    assert payload["error"]["code"] == "missing_index"


def test_match_spectra_reports_missing_index_for_corrupt_file():
    _write_non_sqlite_index()
    load_spectra(str(MSP))
    payload = json.loads(match_spectra())
    assert payload["error"]["code"] == "missing_index"


# --- load_spectra の正直さ（finding 10 / 14） --------------------------------


def _msp_with_an_empty_spectrum(tmp_path):
    """ピークが 1 本も無いレコードを含む MSP を作る。

    MS-DIAL の .dcl は MS/MS の無い MSDecResult も全件持つので、実データでは
    ピーク 0 本のスペクトルが多数を占める。形式に依らない経路を検証したいので
    ここでは MSP で作る。
    """
    path = tmp_path / "with_empty.msp"
    path.write_text(
        "NAME: HasPeaks\n"
        "PRECURSORMZ: 230.0552\n"
        "IONMODE: Positive\n"
        "Num Peaks: 2\n"
        "86.0349 100\n"
        "171.0069 200\n"
        "\n"
        "NAME: NoPeaks\n"
        "PRECURSORMZ: 300.0000\n"
        "IONMODE: Positive\n"
        "Num Peaks: 0\n"
        "\n", encoding="utf-8")
    return path


def test_load_spectra_skips_empty_spectra_and_reports_the_count(tmp_path):
    """ピーク 0 本のスペクトルは載せない。ただし黙って減らさない（finding 10）。"""
    payload = json.loads(load_spectra(str(_msp_with_an_empty_spectrum(tmp_path))))
    assert payload["n_spectra"] == 1
    assert payload["n_skipped_empty"] == 1
    assert payload["spectrum_ids"] == ["with_empty:1"]
    # セッションにも空スペクトルは入らない（match_spectrum の
    # total_query_peaks == 0 の経路に流さない）。
    assert set(session_module.session.spectra) == {"with_empty:1"}


def test_load_spectra_reports_peak_counts(tmp_path):
    """件数だけでは、丸ごと読み違えたファイルが成功に見える（finding 14）。"""
    payload = json.loads(load_spectra(str(MSP)))
    assert payload["n_peak_total"] == 5          # 3 本 + 2 本
    assert payload["n_peak_median"] == 2.5


def test_load_spectra_makes_a_wholly_misparsed_file_visible(tmp_path):
    """ピーク行を 1 本も取れない MSP 変種は、件数ではなくピーク数に出る。

    `m/z intensity;` を 1 行に複数並べる変種はピーク行パーサが読み飛ばす
    （設計上、失敗を飲み込む）。n_spectra だけ見ていると成功に見える。
    """
    path = tmp_path / "semicolon.msp"
    path.write_text(
        "NAME: Variant\n"
        "PRECURSORMZ: 230.0552\n"
        "Num Peaks: 2\n"
        "86.0349 100; 171.0069 200;\n"
        "\n", encoding="utf-8")
    payload = json.loads(load_spectra(str(path)))
    assert payload["n_peak_total"] == 0
    assert payload["n_spectra"] == 0
    assert payload["n_skipped_empty"] == 1


# --- match_spectra の payload 規律と正直さ（finding 8 / 9 / 13） -------------


def test_match_spectra_reports_unknown_spectrum_ids(tmp_path):
    """読み込まれていない id は「照合したが 0 件」ではなく、未ロードとして返す。

    黙って飛ばしたうえで n_spectra に数えると、打ち間違えた id が
    「そのスペクトルは何にも当たらなかった」と読めてしまう（finding 9）。
    """
    _build(tmp_path)
    load_spectra(str(MSP))
    payload = json.loads(match_spectra(spectrum_ids=["sample:1", "typo:99"]))
    assert payload["unknown_spectrum_ids"] == ["typo:99"]
    assert payload["n_spectra"] == 1              # 実際に照合したのは 1 本
    assert payload["n_spectra_requested"] == 2


def test_match_spectra_omits_unknown_ids_key_when_all_are_known(tmp_path):
    """余計なキーで payload を太らせない（spec §12）。"""
    _build(tmp_path)
    load_spectra(str(MSP))
    payload = json.loads(match_spectra())
    assert "unknown_spectrum_ids" not in payload
    assert "truncated" not in payload
    assert payload["n_spectra"] == payload["n_spectra_requested"] == 2


def test_match_spectra_truncates_loudly_when_over_the_row_cap(tmp_path,
                                                              monkeypatch):
    """行数上限に当たったら必ず明示する。黙って切らない（finding 8）。

    実物の .dcl は MSDecResult を数千本持ち、既定の targets は「読み込んだ
    全スペクトル」なので、既定呼び出しがそのまま数万行になりうる。
    """
    monkeypatch.setattr(match_tools, "_MAX_CANDIDATE_ROWS", 1)
    _build(tmp_path)
    load_spectra(str(MSP))
    payload = json.loads(match_spectra())
    assert payload["truncated"] is True
    assert payload["n_spectra"] == 1              # 実際に報告したのは 1 本
    assert payload["n_spectra_requested"] == 2
    assert "1" in payload["truncated_reason"]


def test_match_spectra_reports_the_tolerances_that_produced_it(tmp_path):
    """数字を再現・検算できるよう、効いたトレランスを返す（finding 13）。

    fragment_tolerance は ±1 ビン探索を Da に直したもので、他のどこにも
    出てこない。index_status の bin_width を実際に使える値にする。
    """
    _build(tmp_path)
    load_spectra(str(MSP))
    payload = json.loads(match_spectra(precursor_tolerance=0.02))
    assert payload["precursor_tolerance"] == 0.02
    meta = json.loads(index_status())
    assert payload["fragment_tolerance"] == pytest.approx(
        meta["bin_width"] * 1.5)
