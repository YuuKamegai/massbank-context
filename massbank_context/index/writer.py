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
#: 比較は大小文字非依存で行う（"N/A" と "n/a" のような表記ゆれを両方落とすため）。
_USELESS_CLASS_LABELS = frozenset({"N/A", "NA", "", "ENVIRONMENTAL STANDARD"})

_META_KEYS_JSON = frozenset({"acquisition_classes"})

#: finalize() が DB / writer 自身の状態から導出するキー。caller が meta 引数で
#: これらを渡してはならない。渡された値は索引の中身についての事実ではなく
#: caller の主張にすぎず、無条件に上書きを許すと訂正1が閉じたバグ
#: （メタデータが実際の索引の中身と食い違う）が別の経路で再発する。
_DERIVED_META_KEYS = frozenset({
    "n_records", "n_inchikey", "acquisition_classes", "bin_width",
    "built_at", "resolution_filter",
})


def create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(DDL)
    conn.commit()


def open_index(path: str | Path, create: bool = False) -> sqlite3.Connection:
    """索引を開く。create=True なら無ければ作る。"""
    path = Path(path)
    if not create and not path.exists():
        raise FileNotFoundError(path)
    if create:
        path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(path))
    if create:
        create_schema(conn)
    return conn


def read_meta(conn: sqlite3.Connection) -> dict | None:
    """index_meta を dict にして返す。空なら None。"""
    try:
        rows = conn.execute("SELECT key, value FROM index_meta").fetchall()
    except sqlite3.DatabaseError:
        # OperationalError（表が無い＝未ビルド）の親を捕まえる。索引パスに
        # SQLite でないファイルが在ると "file is not a database" の
        # DatabaseError になり、ビルド中断・コピー中断で実際に生じる。
        # どちらも「索引が使えない」であり、ツールでは直せないので、
        # 呼び出し側は missing_index の封筒に落とす（spec §8）。
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
        #: 層の条件は満たすがピークが 1 本も無かったレコードの件数。
        #: 「受け入れた件数」だけを報告すると、パーサがピークを 1 本も
        #: 拾えていない状態でも索引ビルドは成功したように見える。
        self.n_no_peaks = 0

    def add_record(self, record: SpectrumRecord) -> bool:
        """受け入れたら True、層の条件で弾いたら False を返す。"""
        if not record.inchikey:
            return False
        if record.acquisition_class == "other":
            return False
        if (self._resolution_filter
                and record.resolution_class != self._resolution_filter):
            return False
        if not record.peaks:
            # ピークの無いレコードは照合にも IDF の df にも一切寄与しない。
            # 入れても索引が無意味に膨らむだけで、しかも「受け入れた」件数に
            # 混ざって、パーサがピークを拾えていない失敗を隠してしまう。
            # 層の条件を通ったあとで数えるので、この件数は
            # 「本来入るはずだったのにピークが無い」ものだけを指す。
            self.n_no_peaks += 1
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
            if label and label.upper() not in _USELESS_CLASS_LABELS:
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
        """index_meta を書いてコミットする。

        n_records / n_inchikey / acquisition_classes は DB から数え直す。
        add_record 呼び出しのたびに加算するインスタンスのカウンタでは、
        既存索引への再構築で INSERT OR IGNORE が無視した既存レコードが
        数に入らず、実際には数千件ある索引に n_records=0 と書いてしまう
        （このメタデータを読んだツールがその値を LLM にそのまま報告する）。

        meta に _DERIVED_META_KEYS と衝突するキーがあれば ValueError を投げる。
        DB から導出した事実を caller の文字列で黙って上書きさせない。
        """
        conflicting = _DERIVED_META_KEYS & set(meta)
        if conflicting:
            raise ValueError(
                f"meta に導出済みのキーを渡してはならない: {sorted(conflicting)}。"
                "これらは DB と writer の状態から求める事実であり、上書きすると"
                "索引の中身と食い違うメタデータが記録される。")

        n_records = self._conn.execute(
            "SELECT COUNT(*) FROM records").fetchone()[0]
        n_inchikey = self._conn.execute(
            "SELECT COUNT(DISTINCT inchikey) FROM records").fetchone()[0]
        acquisition_classes = [row[0] for row in self._conn.execute(
            "SELECT DISTINCT acquisition_class FROM records ORDER BY acquisition_class")]

        values = {
            "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "bin_width": repr(self._bin_width),
            "n_records": str(n_records),
            "n_inchikey": str(n_inchikey),
            "acquisition_classes": json.dumps(acquisition_classes),
            "resolution_filter": self._resolution_filter,
        }
        values.update({k: str(v) for k, v in meta.items()})
        for key, value in values.items():
            self._conn.execute(
                "INSERT OR REPLACE INTO index_meta (key, value) VALUES (?,?)",
                (key, value))
        self._conn.commit()
