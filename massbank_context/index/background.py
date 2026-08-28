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

from massbank_context.core.binning import neighbour_bins


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
        for neighbour in neighbour_bins(bin_no):
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
               bin_no: int) -> float:
    """クエリピークが落ちた**中心ビン 1 つ**の idf を返す。該当が無ければ 0.0。

    ここで ±1 に広げてはならない。build_background() が既に各レコードピークを
    ±1 に展開して数えているため、bin k の行は「bin k のクエリピークが ±1 探索で
    当たる unique InChIKey 数」そのものである。すなわち中心ビン 1 行が、
    そのクエリピークの当たりやすさの確率そのものを表している。

    ここで再び (k-1, k, k+1) を引いて MAX(idf)（＝MIN(df)）を採ると、展開済みの
    表の上でもう一度展開することになり、実質 ±2 の窓になる。その結果、
    全 InChIKey に当たるビン（idf=0、何も弁別しない）のクエリピークが、
    たまたま隣にある稀なビンの重みを継承して最大重みを受け取る。
    weighted_score は Σ idf(matched) / Σ idf(query_peaks) なので、これは
    クエリピーク**間**の重み配分を歪め、候補の並び順を変える。

    照合そのものは probe_bins() で ±1 を探り続ける（境界損失を防ぐため）。
    変わるのは重みの引き方だけであり、両者を混同しないこと。
    """
    row = conn.execute(
        "SELECT idf FROM axis_background"
        " WHERE acquisition_class=? AND axis=? AND bin=?",
        (acquisition_class, axis, bin_no)).fetchone()
    return float(row[0]) if row and row[0] is not None else 0.0
