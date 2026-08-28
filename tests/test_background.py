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


def test_lookup_idf_returns_nonzero_value(conn):
    """lookup_idf は非ゼロの idf を返す（単に 0.0 でない、というテスト）。"""
    _insert(conn, "A1", "IK1", [1000])
    _insert(conn, "B1", "IK2", [1000])
    _insert(conn, "C1", "IK3", [5000])
    build_background(conn)
    value = lookup_idf(conn, "LC-ESI-MS2", "product_ion", (999, 1000, 1001))
    assert value > 0


def test_lookup_idf_returns_max_not_sum_of_probe_bins(conn):
    """探索ビンの idf は最大を採る。合計にすると同じフラグメントを重複計上する。

    ±1 展開後の df がビンごとに異なる配置にしないと、MAX と SUM の差が
    定数倍にしかならず、SUM への退行を検出できない。
    """
    _insert(conn, "A1", "IK1", [1000])
    _insert(conn, "A2", "IK2", [1000])
    _insert(conn, "A3", "IK3", [1000])
    _insert(conn, "A4", "IK4", [1002])
    build_background(conn)

    # 展開後: bin 1001 は IK1..IK4 の 4 件(df=4, idf=log(4/4)=0)、
    #         bin 1002 と 1003 は IK4 のみ(df=1, idf=log(4/1))。
    # 探索 (1001, 1002, 1003) の MAX は log(4)、SUM は 2*log(4) で明確に異なる。
    value = lookup_idf(conn, "LC-ESI-MS2", "product_ion", (1001, 1002, 1003))
    assert abs(value - math.log(4)) < 1e-9


def test_ubiquitous_fragment_is_written_with_zero_idf(conn):
    """全化合物に出るフラグメントは、表から消えるのではなく idf=0.0 の行として残る。"""
    _insert(conn, "A1", "IK1", [1000])
    _insert(conn, "A2", "IK2", [1000])
    build_background(conn)
    row = conn.execute(
        "SELECT df_inchikey, idf FROM axis_background WHERE bin=1000").fetchone()
    assert row is not None            # 行が消えていないこと
    assert row[0] == 2                # df == N
    assert row[1] == 0.0              # log(2/2)


def test_lookup_idf_returns_zero_for_unknown_bin(conn):
    _insert(conn, "A1", "IK1", [1000])
    build_background(conn)
    assert lookup_idf(conn, "LC-ESI-MS2", "product_ion", (90000,)) == 0.0
