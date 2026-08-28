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


def test_neutral_loss_is_read_not_computed():
    """保存された mbo:neutral_loss を読むのであって precursor_mz - mz を計算しない(spec §6.3)。

    fixture の値は precursor_mz - mz とほぼ一致するため、丸め誤差でしか両者を判別できない。
    ここでは計算では絶対に得られない値を保存し、読み取り経路であることを曖昧さなく示す。
    """
    ttl = """
    @prefix mbo: <http://www.massbank.jp/ontology/> .
    @prefix sio: <http://semanticscience.org/resource/> .

    [] a mbo:MassSpectrum ;
        mbo:mass_spectram_meta_data [ a mbo:MassSpectramMetaData ;
                mbo:precursor_mz_value 300.0 ] ;
        mbo:has_peak [ sio:SIO_000255 [
                    sio:SIO_000255 [ a mbo:mz ; sio:SIO_000300 100.0 ],
                                    [ a mbo:neutral_loss ; sio:SIO_000300 999.0 ] ] ] .
    """
    graph = rdflib.Graph()
    graph.parse(data=ttl, format="turtle")
    rec = parse_gen_b(graph)[0]
    assert rec.peaks[0].neutral_loss == 999.0     # 300.0 - 100.0 = 200.0 ではない


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


def test_polarity_is_normalised_to_upper_case(tmp_path):
    """世代B も同様に極性を正規化する（世代Aのテストと同じ理由）。"""
    source = GEN_B.read_text(encoding="utf-8").replace(
        'mbo:ion_mode "NEGATIVE"', 'mbo:ion_mode "negative"')
    path = tmp_path / "lowercase.ttl"
    path.write_text(source, encoding="utf-8")
    assert parse_gen_b(_graph(path))[0].polarity == "NEGATIVE"
