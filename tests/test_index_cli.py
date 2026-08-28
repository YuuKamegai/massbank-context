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


def test_build_index_survives_corrupt_ttl(tmp_path, caplog):
    """壊れた TTL 1 個が走査全体を止めない。failed に計上され、健全な分は索引に入る。

    これはこのタスクが単独で負う正しさ要件そのもの。read_ttl_file は
    graph.parse() を包まない設計で、ファイル単位の隔離は build_index の責務。
    """
    ttl_dir = _ttl_dir(tmp_path)
    (ttl_dir / "broken.ttl").write_text(
        "this is not valid turtle @@@ <<<", encoding="utf-8")

    with caplog.at_level("WARNING"):
        stats = build_index(ttl_dir, tmp_path / "index.sqlite")

    assert stats["scanned"] == 3          # 壊れたファイルも走査対象に数える
    assert stats["failed"] == 1
    assert stats["parsed"] == 2           # 健全な 2 件は読めている
    assert stats["accepted"] == 2         # 索引にも入っている
    assert len(caplog.records) == 1
    assert "broken.ttl" in caplog.records[0].message


def test_build_index_is_idempotent(tmp_path):
    """同じ出力先に 2 回走らせても、依存する 4 表のいずれも二重にならない。"""
    ttl_dir = _ttl_dir(tmp_path)
    output = tmp_path / "index.sqlite"

    build_index(ttl_dir, output)
    conn = sqlite3.connect(str(output))
    first = {
        "records": conn.execute("SELECT COUNT(*) FROM records").fetchone()[0],
        "peaks": conn.execute("SELECT COUNT(*) FROM peaks").fetchone()[0],
        "background": conn.execute(
            "SELECT COUNT(*) FROM axis_background").fetchone()[0],
    }
    first_meta = read_meta(conn)
    conn.close()

    build_index(ttl_dir, output)
    conn = sqlite3.connect(str(output))
    second = {
        "records": conn.execute("SELECT COUNT(*) FROM records").fetchone()[0],
        "peaks": conn.execute("SELECT COUNT(*) FROM peaks").fetchone()[0],
        "background": conn.execute(
            "SELECT COUNT(*) FROM axis_background").fetchone()[0],
    }
    second_meta = read_meta(conn)
    conn.close()

    assert second == first
    assert second_meta["n_records"] == first_meta["n_records"]
    assert second_meta["n_inchikey"] == first_meta["n_inchikey"]
    assert second_meta["acquisition_classes"] == first_meta["acquisition_classes"]


def test_main_returns_zero(tmp_path, capsys):
    code = main(["build", "--ttl-dir", str(_ttl_dir(tmp_path)),
                 "--output", str(tmp_path / "index.sqlite")])
    assert code == 0
    assert "accepted" in capsys.readouterr().out
