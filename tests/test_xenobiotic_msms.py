"""MS/MS の共有の筋と順向きの順位。"""
import pytest

from massbank_context.core.binning import to_bin
from massbank_context.core.models import Peak, QuerySpectrum
from massbank_context.index.writer import open_index_readonly
from massbank_context.xenobiotic.msms import forward_rank, reference_bins, shared_evidence
from xenobiotic_helpers import KETAMINE_KEY as KETAMINE, make_index


@pytest.fixture
def index_conn(tmp_path):
    conn = open_index_readonly(make_index(tmp_path / 'index' / 'massbank.sqlite'))
    yield conn
    conn.close()


def _spectrum(precursor, mzs):
    return QuerySpectrum('q', precursor, 'POSITIVE', tuple(Peak(mz, 100.0) for mz in mzs))


def test_reference_bins_exclude_precursor_and_rank_top(index_conn):
    ref = reference_bins(index_conn, KETAMINE, 'POSITIVE', 0.0075, 10)
    assert ref.n_records == 2
    assert ('product_ion', to_bin(238.0993, 0.0075)) not in ref.bins
    assert ref.bins[('product_ion', to_bin(125.0153, 0.0075))].formula == 'C7H6Cl+'
    assert len(reference_bins(index_conn, KETAMINE, 'POSITIVE', 0.0075, 2).top) == 4


def test_parent_and_metabolite_scores(index_conn):
    ref = reference_bins(index_conn, KETAMINE, 'POSITIVE', 0.0075, 10)
    parent = shared_evidence(ref, [_spectrum(238.0990, [125.0153, 207.0571, 220.0888, 238.0990])], 0.0075)
    hnk = shared_evidence(ref, [_spectrum(240.0786, [125.0153, 89.0386, 222.0680])], 0.0075)
    assert parent.score == pytest.approx(0.8) and parent.n_shared == 6
    assert hnk.score == pytest.approx(0.55) and hnk.n_shared == 3
    assert hnk.formulas == ('C7H5+', 'C7H6Cl+')


def test_one_reference_bin_is_used_once(index_conn):
    ref = reference_bins(index_conn, KETAMINE, 'POSITIVE', 0.0075, 10)
    twice = shared_evidence(ref, [_spectrum(None, [125.0150, 125.0156])], 0.0075)
    assert twice.n_shared == 1


def test_empty_and_precursorless_spectra(index_conn):
    ref = reference_bins(index_conn, KETAMINE, 'POSITIVE', 0.0075, 10)
    assert shared_evidence(ref, [_spectrum(238.0990, [])], 0.0075) is None
    no_precursor = shared_evidence(ref, [_spectrum(None, [125.0153])], 0.0075)
    assert no_precursor.n_shared == 1 and no_precursor.score == pytest.approx(0.4)


def test_zero_idf_denominator_gives_zero_score(tmp_path):
    import sqlite3
    path = make_index(tmp_path / 'zero' / 'massbank.sqlite')
    writable = sqlite3.connect(str(path))
    writable.execute("UPDATE axis_background SET idf = 0")
    writable.commit()
    writable.close()
    conn = open_index_readonly(path)
    ref = reference_bins(conn, KETAMINE, 'POSITIVE', 0.0075, 10)
    result = shared_evidence(ref, [_spectrum(238.0990, [125.0153])], 0.0075)
    assert result.score == 0.0 and result.n_shared == 2
    conn.close()


def test_unknown_inchikey_has_no_reference(index_conn):
    assert reference_bins(index_conn, 'AAAAAAAAAAAAAA-AAAAAAAAAA-A', 'POSITIVE', 0.0075, 10) is None


def test_forward_rank(index_conn):
    assert forward_rank(index_conn, [_spectrum(238.0990, [125.0153, 207.0571])], KETAMINE) == 1
    assert forward_rank(index_conn, [_spectrum(500.3000, [125.0153])], KETAMINE) is None
