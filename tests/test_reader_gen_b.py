from pathlib import Path

import rdflib

from massbank_context.ttl.generation import detect_generation
from massbank_context.ttl.reader import read_ttl_file
from massbank_context.ttl.reader_gen_b import parse_gen_b

FIXTURES = Path(__file__).parent / "fixtures" / "ttl"
GEN_B = FIXTURES / "gen_b" / "MSBNK-TEST-B0001.ttl"
GEN_A = FIXTURES / "gen_a" / "MSBNK-TEST-A0001.ttl"


def _graph(path) -> rdflib.Graph:
    g = rdflib.Graph()
    g.parse(path, format="turtle")
    return g


def test_detects_generation_b():
    assert detect_generation(_graph(GEN_B)) == "gen_b"


def test_parses_metadata():
    rec = parse_gen_b(_graph(GEN_B))[0]
    assert rec.accession == "MSBNK-TEST-B0001"
    assert rec.inchikey == "XHKOSGCDDIGCOM-HIAJHCFCSA-N"
    assert rec.precursor_mz == 267.0855
    assert rec.polarity == "NEGATIVE"
    assert rec.acquisition_class == "LC-ESI-MS2"


def test_reads_neutral_loss_directly():
    """世代B は mbo:neutral_loss を持つので算出不要(spec §6.3)。"""
    rec = parse_gen_b(_graph(GEN_B))[0]
    by_mz = {round(p.mz, 4): p.neutral_loss for p in rec.peaks}
    assert by_mz[135.0434] == 132.0421
    assert by_mz[110.0349] == 157.0506


def test_peaks_sorted_and_complete():
    rec = parse_gen_b(_graph(GEN_B))[0]
    mzs = [p.mz for p in rec.peaks]
    assert mzs == sorted(mzs)
    assert len(mzs) == 2


def test_read_ttl_file_dispatches_gen_a():
    records = read_ttl_file(GEN_A)
    assert records[0].accession == "MSBNK-TEST-A0001"


def test_read_ttl_file_dispatches_gen_b():
    records = read_ttl_file(GEN_B)
    assert records[0].accession == "MSBNK-TEST-B0001"


def test_read_ttl_file_returns_empty_for_unknown(tmp_path):
    path = tmp_path / "empty.ttl"
    path.write_text("@prefix ex: <http://example.org/> .\n", encoding="utf-8")
    assert read_ttl_file(path) == []
