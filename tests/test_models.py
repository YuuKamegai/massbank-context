import json

from massbank_context.core.errors import (
    incompatible_index,
    missing_index,
    missing_state,
)
from massbank_context.core.models import Peak, QuerySpectrum, SpectrumRecord
from massbank_context.core.serialization import json_payload, round_floats, as_tsv


def _record(**overrides) -> SpectrumRecord:
    base = dict(
        accession="MSBNK-TEST-0001",
        inchikey="XSDQTOBWRPYKKA-UHFFFAOYSA-N",
        precursor_mz=230.0552,
        precursor_type="[M+H]+",
        polarity="POSITIVE",
        instrument_type="LC-ESI-ITFT",
        ms_type="MS2",
        resolution_value=15000.0,
        collision_energy="35% (nominal)",
        splash="splash10-001i-0290000000-657972e0af96e9da0162",
        license="CC0",
        molecular_formula="C6H8ClN7O",
        exact_mass=229.0479,
        compound_names=("Amiloride",),
        ch_compound_class="N/A; Environmental Standard",
        peaks=(Peak(mz=86.0349, relative_intensity=11.0,
                    formula="C2H4N3O+", neutral_loss=None),),
    )
    base.update(overrides)
    return SpectrumRecord(**base)


def test_record_derives_layers():
    rec = _record()
    assert rec.acquisition_class == "LC-ESI-MS2"
    assert rec.resolution_class == "high"


def test_record_without_precursor_is_gc_ei():
    rec = _record(instrument_type="EI-B", ms_type="MS", precursor_mz=None,
                  precursor_type=None)
    assert rec.acquisition_class == "GC-EI-MS"


def test_json_payload_has_no_whitespace():
    """戻り値の空白は LLM の文脈をそのまま食う（spec §12）。"""
    out = json_payload({"a": 1, "b": [1, 2]})
    assert out == '{"a":1,"b":[1,2]}'
    assert " " not in out


def test_round_floats_preserves_bool():
    """bool は int の派生。数値として丸めると True が 1 に劣化する。"""
    result = round_floats({"x": 1.23456789, "flag": True})
    assert result["x"] == 1.234568
    assert result["flag"] is True
    assert type(result["flag"]) is bool


def test_missing_index_envelope():
    payload = json.loads(missing_index("CLI: python -m massbank_context.index build"))
    assert payload["error"]["code"] == "missing_index"
    assert "required_action" in payload["error"]


def test_missing_state_envelope():
    payload = json.loads(missing_state("spectra", ["load_spectra"], "説明"))
    assert payload["error"]["code"] == "missing_state"
    assert payload["error"]["required_tools"] == ["load_spectra"]


def test_incompatible_index_envelope():
    payload = json.loads(
        incompatible_index("LC-ESI-MS2", ["GC-EI-MS"], "説明"))
    assert payload["error"]["code"] == "incompatible_index"
    assert payload["error"]["query_acquisition_class"] == "LC-ESI-MS2"
    assert payload["error"]["available"] == ["GC-EI-MS"]


def test_query_spectrum_holds_source():
    q = QuerySpectrum(spectrum_id="q1", precursor_mz=230.0552,
                      polarity="POSITIVE",
                      peaks=(Peak(86.0349, 11.0, None, None),),
                      source="sample.msp")
    assert q.source == "sample.msp"
    assert len(q.peaks) == 1


def test_as_tsv_header_emitted_once():
    """ヘッダー行は 1 回だけ、None は空フィールドで返す。"""
    rows = [
        ["mz1", "int1"],
        ["mz2", None],
        ["mz3", "int3"],
    ]
    result = as_tsv(["m/z", "intensity"], rows)
    lines = result.split("\n")
    assert lines[0] == "m/z\tintensity"
    assert lines[1] == "mz1\tint1"
    assert lines[2] == "mz2\t"
    assert lines[3] == "mz3\tint3"


def test_as_tsv_escapes_tab_and_newline():
    """値中のタブ・改行で列や行がずれない。"""
    result = as_tsv(["name", "n"], [["a\tb\nc\\d", 1]])
    lines = result.split("\n")
    assert len(lines) == 2
    assert lines[1].split("\t") == ["a\\tb\\nc\\\\d", "1"]


def test_rows_tsv_unions_columns_and_serializes_nested():
    from massbank_context.core.serialization import rows_tsv
    result = rows_tsv([{"a": 1, "b": {"x": 1}}, {"a": 2, "c": [1, 2]}])
    assert result.split("\n") == ["a\tb\tc", '1\t{"x":1}\t', "2\t\t[1,2]"]


def test_round_floats_keep_small_preserves_tiny_values():
    assert round_floats({"p": 2.1e-13}) == {"p": 0.0}
    assert round_floats({"p": 2.1e-13, "fc": 1.23456789}, keep_small=True) == {"p": 2.1e-13, "fc": 1.234568}
