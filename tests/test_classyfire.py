import sqlite3

import pytest

from massbank_context.core.models import Peak, SpectrumRecord
from massbank_context.index.classyfire import CLASSYFIRE_LEVELS, read_classyfire_tsv
from massbank_context.index.writer import IndexWriter, create_schema

# 実物と同じ 6 列・ヘッダ無し。同一 InChIKey がレコード数だけ重複して現れる。
SAMPLE = (
    "https://massbank.eu/MassBank/RecordDisplay.jsp?id=WA002994\t"
    "JFPVXVDWJQMJEE-IZRZKJBUSA-N\tOrganic compounds\t"
    "Organoheterocyclic compounds\tLactams\tBeta lactams\n"
    "https://massbank.eu/MassBank/RecordDisplay.jsp?id=WA002993\t"
    "JFPVXVDWJQMJEE-IZRZKJBUSA-N\tOrganic compounds\t"
    "Organoheterocyclic compounds\tLactams\tBeta lactams\n"
    "https://massbank.eu/MassBank/RecordDisplay.jsp?id=XX000001\t"
    "AAAAAAAAAAAAAA-UHFFFAOYSA-N\tOrganic compounds\tLipids\tN/A\t\n"
)


@pytest.fixture()
def tsv(tmp_path):
    path = tmp_path / "classyfire.tsv"
    path.write_text(SAMPLE, encoding="utf-8")
    return path


def _record(inchikey="JFPVXVDWJQMJEE-IZRZKJBUSA-N") -> SpectrumRecord:
    return SpectrumRecord(
        accession="A1", inchikey=inchikey, precursor_mz=230.0552,
        precursor_type="[M+H]+", polarity="POSITIVE",
        instrument_type="LC-ESI-QTOF", ms_type="MS2", resolution_value=35000.0,
        collision_energy="35", splash=None, license="CC0",
        molecular_formula="C6H8ClN7O", exact_mass=229.0479,
        compound_names=("X",), ch_compound_class=None,
        peaks=(Peak(86.0349, 11.0, None, None),),
    )


@pytest.fixture()
def conn():
    connection = sqlite3.connect(":memory:")
    create_schema(connection)
    yield connection
    connection.close()


def test_levels_become_separate_sources(tsv):
    """階層は source に持たせる。DDL を変えずに 4 レベルを区別できる。"""
    entries = set(read_classyfire_tsv(tsv))
    assert ("JFPVXVDWJQMJEE-IZRZKJBUSA-N", "classyfire:kingdom",
            "Organic compounds") in entries
    assert ("JFPVXVDWJQMJEE-IZRZKJBUSA-N", "classyfire:subclass",
            "Beta lactams") in entries
    assert CLASSYFIRE_LEVELS == ("kingdom", "superclass", "class", "subclass")


def test_duplicate_rows_collapse(tsv):
    """同一 InChIKey がレコード数だけ重複して現れるので必ず畳む。"""
    entries = list(read_classyfire_tsv(tsv))
    assert len(entries) == len(set(entries))
    assert sum(1 for ik, _, _ in entries
               if ik == "JFPVXVDWJQMJEE-IZRZKJBUSA-N") == 4


def test_na_and_empty_labels_are_dropped(tsv):
    """'N/A' と空欄はラベルではない。層別の軸に載せると偽のクラスができる。"""
    labels = {label for ik, _, label in read_classyfire_tsv(tsv)
              if ik == "AAAAAAAAAAAAAA-UHFFFAOYSA-N"}
    assert labels == {"Organic compounds", "Lipids"}


def test_import_only_covers_indexed_inchikeys(conn, tsv):
    """索引に入っていない InChIKey のクラスは持たない。

    TSV は MassBank 全レコードを覆うが、索引は層の条件で絞られている。
    入っていない化合物のクラスを持っても層別背景には使えず、
    compound_class が何を表す表なのかが曖昧になる。
    """
    writer = IndexWriter(conn, bin_width=0.01)
    writer.add_record(_record())
    stats = writer.add_compound_classes(read_classyfire_tsv(tsv))

    assert stats == {"written": 4, "skipped_unknown_inchikey": 2}
    keys = {row[0] for row in conn.execute(
        "SELECT DISTINCT inchikey FROM compound_class")}
    assert keys == {"JFPVXVDWJQMJEE-IZRZKJBUSA-N"}


def test_import_is_idempotent(conn, tsv):
    writer = IndexWriter(conn, bin_width=0.01)
    writer.add_record(_record())
    writer.add_compound_classes(read_classyfire_tsv(tsv))
    writer.add_compound_classes(read_classyfire_tsv(tsv))
    assert conn.execute("SELECT COUNT(*) FROM compound_class").fetchone()[0] == 4
