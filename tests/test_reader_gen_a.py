import logging
from pathlib import Path

import rdflib

from massbank_context.ttl.generation import detect_generation
from massbank_context.ttl.reader_gen_a import _number, parse_gen_a

FIXTURE = Path(__file__).parent / "fixtures" / "ttl" / "gen_a" / "MSBNK-TEST-A0001.ttl"


def _graph() -> rdflib.Graph:
    g = rdflib.Graph()
    g.parse(FIXTURE, format="turtle")
    return g


def test_detects_generation_a():
    assert detect_generation(_graph()) == "gen_a"


def test_parses_single_record():
    records = parse_gen_a(_graph())
    assert len(records) == 1


def test_record_metadata():
    rec = parse_gen_a(_graph())[0]
    assert rec.accession == "MSBNK-TEST-A0001"
    assert rec.inchikey == "XSDQTOBWRPYKKA-UHFFFAOYSA-N"
    assert rec.precursor_mz == 230.0552
    assert rec.precursor_type == "[M+H]+"
    assert rec.polarity == "POSITIVE"
    assert rec.instrument_type == "LC-ESI-ITFT"
    assert rec.ms_type == "MS2"
    assert rec.resolution_value == 15000.0
    assert rec.molecular_formula == "C6H8ClN7O"
    assert rec.splash.startswith("splash10-")


def test_peaks_are_sorted_by_mz():
    rec = parse_gen_a(_graph())[0]
    mzs = [p.mz for p in rec.peaks]
    assert mzs == sorted(mzs)
    assert len(mzs) == 3


def test_tentative_formula_is_joined_by_mz():
    """has_peak と has_peak_annotations は別ブロックで、m/z で突き合わせる。"""
    rec = parse_gen_a(_graph())[0]
    by_mz = {round(p.mz, 4): p.formula for p in rec.peaks}
    assert by_mz[86.0349] == "C2H4N3O+"
    assert by_mz[171.0069] == "C5H4ClN4O+"
    assert by_mz[230.0552] is None      # アノテーションが無いピーク


def test_neutral_loss_is_none_before_derivation():
    """世代A では mbo:neutral_loss が無い。算出は Task 5 が行う。"""
    rec = parse_gen_a(_graph())[0]
    assert all(p.neutral_loss is None for p in rec.peaks)


def test_layers_are_derived():
    rec = parse_gen_a(_graph())[0]
    assert rec.acquisition_class == "LC-ESI-MS2"
    assert rec.resolution_class == "high"


def test_remaining_record_fields_are_read():
    """fixture に値があるのにこれまで検証されていなかったフィールド。"""
    rec = parse_gen_a(_graph())[0]
    assert rec.collision_energy == "35% (nominal)"
    assert rec.exact_mass == 229.0479
    assert rec.compound_names == ("Amiloride",)
    assert rec.ch_compound_class == "N/A; Environmental Standard"
    assert rec.license is not None and "creativecommons.org" in rec.license


def test_number_logs_warning_on_invalid_literal(caplog):
    """値が存在するのに float に変換できない場合だけ warning を出す。

    欠損（None）は正常系なので warning は出ない。異常なリテラルだけを
    区別して記録できることを確認する。
    """
    with caplog.at_level(logging.WARNING,
                          logger="massbank_context.ttl.reader_gen_a"):
        assert _number("not-a-number") is None
        assert _number(None) is None

    assert len(caplog.records) == 1
    assert "float に変換できない" in caplog.records[0].message
