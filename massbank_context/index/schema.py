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

-- polarity は候補プールの絞り込み条件（spec §6.8 手順 1）なので、背景も同じ
-- 軸で層別する。極性を跨いで数えた df は「いま探している集団の中での稀さ」に
-- ならない。極性不明のクエリ（.dcl の ion_mode=2）はプールが両極性にまたがる
-- ため、その背景として polarity='*'（POOLED_POLARITY）の行を別に持つ。
CREATE TABLE IF NOT EXISTS axis_background (
    acquisition_class TEXT NOT NULL,
    polarity          TEXT NOT NULL,
    axis              TEXT NOT NULL,
    bin               INTEGER NOT NULL,
    df_inchikey       INTEGER NOT NULL,
    idf               REAL NOT NULL,
    PRIMARY KEY (acquisition_class, polarity, axis, bin)
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
