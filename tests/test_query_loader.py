from pathlib import Path

import pytest

from massbank_context.query.dcl import _polarity_from_ion_mode
from massbank_context.query.loader import load_spectra_file
from massbank_context.query.msp import read_msp

FIXTURE = Path(__file__).parent / "fixtures" / "query" / "sample.msp"


def test_reads_two_spectra():
    assert len(read_msp(FIXTURE)) == 2


def test_first_spectrum_fields():
    spectrum = read_msp(FIXTURE)[0]
    assert spectrum.precursor_mz == 230.0552
    assert spectrum.polarity == "POSITIVE"
    assert spectrum.retention_time == 0.876
    assert len(spectrum.peaks) == 3


def test_polarity_is_normalised_to_upper():
    """MSP の 'Positive' / 'Negative' を索引側の表記に合わせる。"""
    spectra = read_msp(FIXTURE)
    assert spectra[0].polarity == "POSITIVE"
    assert spectra[1].polarity == "NEGATIVE"


def test_relative_intensity_is_normalised_to_999():
    """装置間で絶対強度の桁が違うため、基準ピーク 999 に正規化する。"""
    spectrum = read_msp(FIXTURE)[0]
    assert max(p.relative_intensity for p in spectrum.peaks) == 999.0


def test_peaks_sorted_by_mz():
    spectrum = read_msp(FIXTURE)[0]
    mzs = [p.mz for p in spectrum.peaks]
    assert mzs == sorted(mzs)


def test_spectrum_ids_are_unique_and_stable():
    spectra = read_msp(FIXTURE)
    ids = [s.spectrum_id for s in spectra]
    assert ids == ["sample:1", "sample:2"]


def test_loader_dispatches_by_extension():
    assert len(load_spectra_file(FIXTURE)) == 2


def test_loader_rejects_unknown_extension(tmp_path):
    path = tmp_path / "x.csv"
    path.write_text("a,b\n", encoding="utf-8")
    with pytest.raises(ValueError, match="対応していない"):
        load_spectra_file(path)


@pytest.mark.parametrize(
    ("ion_mode", "expected"),
    [
        (0, "POSITIVE"),
        (1, "NEGATIVE"),
        (2, None),
        (99, None),
    ],
)
def test_polarity_from_ion_mode(ion_mode, expected):
    """MS-DIAL の IonMode 列挙は 0=Positive/1=Negative/2=Both の3値。

    Both と未知の値は極性で絞り込まないよう None にする（docs/output_format/dcl.md）。
    """
    assert _polarity_from_ion_mode(ion_mode) == expected
