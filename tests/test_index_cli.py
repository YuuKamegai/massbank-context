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
