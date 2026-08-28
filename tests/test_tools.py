import json
import shutil
import struct
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
