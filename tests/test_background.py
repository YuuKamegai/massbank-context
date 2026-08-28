import math
import sqlite3

import pytest

from massbank_context.index.background import (
    POOLED_POLARITY,
    build_background,
    lookup_idf,
)
from massbank_context.index.writer import create_schema


@pytest.fixture()
def conn():
    connection = sqlite3.connect(":memory:")
    create_schema(connection)
    yield connection
    connection.close()


def _insert(conn, accession, inchikey, bins, axis="product_ion",
            acquisition_class="LC-ESI-MS2", polarity="POSITIVE"):
    cursor = conn.execute(
        "INSERT INTO records (accession, inchikey, acquisition_class,"
        " resolution_class, polarity) VALUES (?,?,?,'high',?)",
        (accession, inchikey, acquisition_class, polarity))
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
    value = lookup_idf(conn, "LC-ESI-MS2", "POSITIVE", "product_ion", 1000)
    assert value > 0


def _four_inchikey_layout(conn):
    """IK1..IK3 が bin 1000、IK4 が bin 1002。N=4。

    ±1 展開後の各行の意味:
      bin 1000: IK1..IK3           (df=3)
      bin 1001: IK1..IK3 + IK4     (df=4)  ← 全 InChIKey が当たる
      bin 1002: IK4                (df=1)
      bin 1003: IK4                (df=1)  ← 4 件中 1 件しか当たらない
    """
    _insert(conn, "A1", "IK1", [1000])
    _insert(conn, "A2", "IK2", [1000])
    _insert(conn, "A3", "IK3", [1000])
    _insert(conn, "A4", "IK4", [1002])
    build_background(conn)


def test_ubiquitous_query_bin_weighs_nothing(conn):
    """全 InChIKey が持つビンに落ちるクエリピークは何も弁別しない＝ idf 0.0。

    背景表は既に ±1 展開済みなので、bin 1001 の行は「bin 1001 のクエリピークが
    ±1 探索で当たる InChIKey 数」そのものである。ここで再度 ±1 に広げると
    最も稀な隣（bin 1002 の df=1）の重みを継承してしまい、4 件中 4 件に
    当たるピークが最大重みを受け取る。
    """
    _four_inchikey_layout(conn)
    assert lookup_idf(conn, "LC-ESI-MS2", "POSITIVE", "product_ion", 1001) == 0.0


def test_rare_query_bin_weighs_log_n(conn):
    """N 件中 1 件しか持たないビンに落ちるクエリピークは log(N)。"""
    _four_inchikey_layout(conn)
    value = lookup_idf(conn, "LC-ESI-MS2", "POSITIVE", "product_ion", 1003)
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
    assert lookup_idf(conn, "LC-ESI-MS2", "POSITIVE", "product_ion", 90000) == 0.0


# --- 極性による層別（spec §6.8 手順 1） --------------------------------------
#
# 候補プールは (acquisition_class, polarity) で絞られる。重みの出所である
# 背景がそれと違う集団（極性混合）だと、「いま探している集団の中での稀さ」
# という idf の意味そのものが壊れる。


def _mixed_polarity_layout(conn):
    """POSITIVE 2 件・NEGATIVE 1 件。プールごとに df/N が変わる配置。

      POSITIVE プール (N=2): bin 1000 -> IK1 のみ, bin 2000 -> IK2 のみ
      NEGATIVE プール (N=1): bin 1000 -> IK3 のみ
      プール "*"    (N=3): bin 1000 -> IK1,IK3, bin 2000 -> IK2
    """
    _insert(conn, "P1", "IK1", [1000], polarity="POSITIVE")
    _insert(conn, "P2", "IK2", [2000], polarity="POSITIVE")
    _insert(conn, "N1", "IK3", [1000], polarity="NEGATIVE")
    build_background(conn)


def test_background_is_stratified_by_polarity(conn):
    _mixed_polarity_layout(conn)
    polarities = {row[0] for row in conn.execute(
        "SELECT DISTINCT polarity FROM axis_background")}
    assert polarities == {"POSITIVE", "NEGATIVE", POOLED_POLARITY}


def test_polarity_stratified_idf_uses_only_its_own_pool(conn):
    """POSITIVE プールでは bin 1000 は 2 件中 1 件 -> log(2)。"""
    _mixed_polarity_layout(conn)
    assert abs(lookup_idf(conn, "LC-ESI-MS2", "POSITIVE", "product_ion", 1000)
               - math.log(2)) < 1e-9
    # NEGATIVE プールは 1 件中 1 件なので何も弁別しない。
    assert lookup_idf(conn, "LC-ESI-MS2", "NEGATIVE", "product_ion", 1000) == 0.0


def test_pooled_row_covers_all_polarities(conn):
    """極性不明のクエリ（.dcl の ion_mode=2）のプールは両極性にまたがる。

    その背景も両極性をまたいだ集団で数えた行でなければならない。
    """
    _mixed_polarity_layout(conn)
    row = conn.execute(
        "SELECT df_inchikey, idf FROM axis_background"
        " WHERE polarity=? AND axis='product_ion' AND bin=1000",
        (POOLED_POLARITY,)).fetchone()
    assert row[0] == 2                               # IK1(POS) と IK3(NEG)
    assert abs(row[1] - math.log(3 / 2)) < 1e-9      # N=3


def test_lookup_idf_with_polarity_none_uses_pooled_rows(conn):
    _mixed_polarity_layout(conn)
    assert abs(lookup_idf(conn, "LC-ESI-MS2", None, "product_ion", 1000)
               - math.log(3 / 2)) < 1e-9


def test_polarity_less_record_is_counted_only_in_the_pooled_rows(conn):
    """polarity NULL のレコードは極性で絞ったプールからは引けない。

    候補側のクエリ（polarity=? OR ? IS NULL）が NULL 極性のレコードを
    POSITIVE クエリで拾わない以上、背景側も拾ってはならない。
    """
    _insert(conn, "P1", "IK1", [1000], polarity="POSITIVE")
    _insert(conn, "U1", "IK2", [1000], polarity=None)
    build_background(conn)
    positive = conn.execute(
        "SELECT df_inchikey FROM axis_background"
        " WHERE polarity='POSITIVE' AND axis='product_ion' AND bin=1000").fetchone()
    pooled = conn.execute(
        "SELECT df_inchikey FROM axis_background"
        " WHERE polarity=? AND axis='product_ion' AND bin=1000",
        (POOLED_POLARITY,)).fetchone()
    assert positive[0] == 1
    assert pooled[0] == 2
