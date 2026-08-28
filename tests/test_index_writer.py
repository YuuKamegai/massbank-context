import sqlite3

import pytest

from massbank_context.core.models import Peak, SpectrumRecord
from massbank_context.index.writer import IndexWriter, create_schema, read_meta


def _record(accession="A1", inchikey="IK1", instrument_type="LC-ESI-QTOF",
            ms_type="MS2", precursor_mz=230.0552) -> SpectrumRecord:
    return SpectrumRecord(
        accession=accession, inchikey=inchikey, precursor_mz=precursor_mz,
        precursor_type="[M+H]+", polarity="POSITIVE",
        instrument_type=instrument_type, ms_type=ms_type, resolution_value=35000.0,
        collision_energy="35", splash="splash10-xxx", license="CC0",
        molecular_formula="C6H8ClN7O", exact_mass=229.0479,
        compound_names=("Amiloride",), ch_compound_class="Natural Product; Lipid",
        peaks=(Peak(86.0349, 11.0, "C2H4N3O+", None),
               Peak(171.0069, 226.0, None, None)),
    )


@pytest.fixture()
def conn():
    connection = sqlite3.connect(":memory:")
    create_schema(connection)
    yield connection
    connection.close()


def test_schema_creates_all_tables(conn):
    names = {row[0] for row in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    assert {"records", "peaks", "axis_background", "compound_class",
            "compounds", "index_meta", "mass_error_profile"} <= names


def test_add_record_stores_layers(conn):
    writer = IndexWriter(conn, bin_width=0.01)
    assert writer.add_record(_record()) is True
    row = conn.execute(
        "SELECT acquisition_class, resolution_class FROM records").fetchone()
    assert row == ("LC-ESI-MS2", "high")


def test_add_record_writes_both_axes(conn):
    writer = IndexWriter(conn, bin_width=0.01)
    writer.add_record(_record())
    axes = {row[0] for row in conn.execute("SELECT DISTINCT axis FROM peaks")}
    assert axes == {"product_ion", "neutral_loss"}


def test_low_resolution_is_rejected_by_default(conn):
    """既定は resolution_class == 'high' のみ（spec §6.2）。"""
    writer = IndexWriter(conn, bin_width=0.01)
    assert writer.add_record(_record(instrument_type="LC-ESI-QQ")) is False
    assert conn.execute("SELECT COUNT(*) FROM records").fetchone()[0] == 0


def test_other_acquisition_class_is_rejected(conn):
    writer = IndexWriter(conn, bin_width=0.01)
    assert writer.add_record(_record(ms_type="MS3")) is False


def test_record_without_inchikey_is_rejected(conn):
    """InChIKey が無いレコードは IDF の分母にも候補にもならない。"""
    writer = IndexWriter(conn, bin_width=0.01)
    assert writer.add_record(_record(inchikey=None)) is False


def test_compound_class_is_split_on_semicolon(conn):
    writer = IndexWriter(conn, bin_width=0.01)
    writer.add_record(_record())
    labels = {row[0] for row in conn.execute(
        "SELECT label FROM compound_class")}
    assert labels == {"Natural Product", "Lipid"}


def test_na_compound_class_is_not_stored(conn):
    writer = IndexWriter(conn, bin_width=0.01)
    record = _record()
    writer.add_record(SpectrumRecord(**{
        **record.__dict__, "ch_compound_class": "N/A; Environmental Standard"}))
    labels = {row[0] for row in conn.execute("SELECT label FROM compound_class")}
    assert "N/A" not in labels


def test_finalize_writes_meta(conn):
    writer = IndexWriter(conn, bin_width=0.004)
    writer.add_record(_record())
    writer.finalize({"ttl_source": "/tmp/ttl", "ttl_generation": "gen_a",
                     "resolution_filter": "high"})
    meta = read_meta(conn)
    assert meta["bin_width"] == 0.004
    assert meta["n_records"] == 1
    assert meta["n_inchikey"] == 1
    assert meta["acquisition_classes"] == ["LC-ESI-MS2"]
    assert meta["built_at"]


def test_read_meta_returns_none_when_absent(conn):
    assert read_meta(conn) is None


def test_finalize_after_rebuild_counts_from_database(conn):
    """再構築時、2 回目の add_record は INSERT OR IGNORE で無視され False を返すが、
    finalize() が書く n_records / n_inchikey / acquisition_classes は
    インスタンスのカウンタではなく DB の実データから導出されなければならない
    （補正1: 既存索引への再構築で n_records=0 と誤記されるバグの防止）。
    """
    writer_a = IndexWriter(conn, bin_width=0.01)
    writer_a.add_record(_record())
    writer_a.finalize({"resolution_filter": "high"})

    # 別インスタンス（例: 再実行のプロセス）で同じレコードを再投入する。
    writer_b = IndexWriter(conn, bin_width=0.01)
    assert writer_b.add_record(_record()) is False  # 既存 accession なので無視される
    writer_b.finalize({"resolution_filter": "high"})

    meta = read_meta(conn)
    assert meta["n_records"] == 1
    assert meta["n_inchikey"] == 1
    assert meta["acquisition_classes"] == ["LC-ESI-MS2"]
