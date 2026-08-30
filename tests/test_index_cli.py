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


# --- 世代ごとの内訳とピーク 0 本の検出（finding 7 / 12） ---------------------

#: 世代B と判定されるが、has_peak からは SIO ノードに辿り着けない TTL。
#: 実物の世代B が想定と違う入れ子だった場合に起きる状態を再現する
#: （spec §11 確認事項 1。同梱 fixture は編集せず、ここで組み立てる）。
def _write_peakless_gen_b(ttl_dir: Path) -> None:
    source = (FIXTURES / "gen_b" / "MSBNK-TEST-B0001.ttl").read_text(
        encoding="utf-8")
    text = (source.replace("MSBNK-TEST-B0001", "MSBNK-TEST-B0002")
                  .replace("mbo:has_peak [", "mbo:has_peaks ["))
    (ttl_dir / "MSBNK-TEST-B0002.ttl").write_text(text, encoding="utf-8")


def test_build_index_reports_counts_per_generation(tmp_path):
    """どの世代を何ファイル読み、ピークを何本拾ったかが一目で分かること。"""
    stats = build_index(_ttl_dir(tmp_path), tmp_path / "index.sqlite")
    assert stats["generations"] == {
        "gen_a": {"files": 1, "records": 1, "peaks": 3},
        "gen_b": {"files": 1, "records": 1, "peaks": 2},
    }


def test_build_index_counts_peakless_records_separately(tmp_path):
    """ピーク 0 本のレコードは accepted に混ぜず no_peaks に数える。"""
    ttl_dir = _ttl_dir(tmp_path)
    _write_peakless_gen_b(ttl_dir)
    stats = build_index(ttl_dir, tmp_path / "index.sqlite")
    assert stats["parsed"] == 3
    assert stats["accepted"] == 2
    assert stats["no_peaks"] == 1
    # 世代ごとの内訳を見れば、ピークを拾えていないのが gen_b 側だと分かる。
    assert stats["generations"]["gen_b"] == {
        "files": 2, "records": 2, "peaks": 2}


def test_build_index_records_ttl_generation_in_meta(tmp_path):
    """spec §6.7 の index_meta.ttl_generation。混在ビルドでは全部を記録する。"""
    build_index(_ttl_dir(tmp_path), tmp_path / "index.sqlite")
    conn = sqlite3.connect(str(tmp_path / "index.sqlite"))
    meta = read_meta(conn)
    conn.close()
    assert meta["ttl_generation"] == "gen_a,gen_b"


def test_build_index_logs_bin_width_at_info(tmp_path, caplog):
    """採用ビン幅と装置クラス別プロファイルを運用者が目視できること。"""
    with caplog.at_level("INFO", logger="massbank_context.index.__main__"):
        build_index(_ttl_dir(tmp_path), tmp_path / "index.sqlite")
    messages = [record.getMessage() for record in caplog.records]
    assert any("採用ビン幅" in message for message in messages)
    assert any("TTL 世代 gen_a" in message for message in messages)


def test_bin_width_ignores_records_excluded_from_the_index(tmp_path):
    """ビン幅は索引に入るレコードだけで測る。

    実コーパスでは、resolution_filter で弾かれる低分解能クラス
    （GC-EI-Q / LC-APPI-QQ 等）の散らばりが最大値経由でビン幅を決めていた。
    索引に 1 件も入らないレコードが索引の照合窓を決めるのは誤りである。
    """
    ttl_dir = _ttl_dir(tmp_path)
    source = (FIXTURES / "gen_a" / "MSBNK-TEST-A0001.ttl").read_text(
        encoding="utf-8")
    # 低分解能（弾かれる）レコードを 2 件。同一 InChIKey・同じフラグメントを
    # 大きくずらして置き、質量誤差の標本として最大幅を作る。
    for i, shifted in enumerate(("8.60349e+01", "8.64349e+01"), start=1):
        (ttl_dir / f"low{i}.ttl").write_text(
            source.replace("MSBNK-TEST-A0001", f"MSBNK-TEST-A010{i}")
                  .replace("LC-ESI-ITFT", "LC-ESI-QQ")
                  .replace("mbo:mz 8.60349e+01", f"mbo:mz {shifted}"),
            encoding="utf-8")

    build_index(ttl_dir, tmp_path / "index.sqlite")
    conn = sqlite3.connect(str(tmp_path / "index.sqlite"))
    classes = {row[0] for row in conn.execute(
        "SELECT DISTINCT instrument_class FROM mass_error_profile")}
    indexed = {row[0] for row in conn.execute(
        "SELECT DISTINCT instrument_type FROM records")}
    conn.close()
    assert "LC-ESI-QQ" not in classes
    assert classes <= indexed
