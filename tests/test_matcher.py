import shutil
import sqlite3
from pathlib import Path

import pytest

from massbank_context.core.models import Peak, QuerySpectrum
from massbank_context.index.__main__ import build_index
from massbank_context.spectral.matcher import derive_band, match_spectrum

FIXTURES = Path(__file__).parent / "fixtures" / "ttl"


@pytest.fixture()
def conn(tmp_path):
    ttl_dir = tmp_path / "ttl"
    ttl_dir.mkdir()
    shutil.copy(FIXTURES / "gen_a" / "MSBNK-TEST-A0001.ttl", ttl_dir)

    # 訂正: 索引を 1 レコードだけで作ると、単一の InChIKey しかないため
    # 全ビンの idf が log(1/1)=0 になり、weighted_score が常に 0.0 になって
    # 主スコアの計算式そのものが一切検証されない。2 レコード目を InChIKey・
    # precursor_mz を変えて追加し、IDF の背景を意味のあるものにする。
    # precursor_mz は既存クエリ（precursor 230.0552）の候補には決して現れない
    # 値にする。500.0 は test_precursor_mismatch_returns_nothing が「合わない
    # precursor」の代表値として既に使っているため、それと衝突しない 600.0 を選ぶ。
    source = (FIXTURES / "gen_a" / "MSBNK-TEST-A0001.ttl").read_text(
        encoding="utf-8")
    second = (source
              .replace("MSBNK-TEST-A0001", "MSBNK-TEST-A0002")
              .replace("XSDQTOBWRPYKKA-UHFFFAOYSA-N", "YYYYYYYYYYYYYY-UHFFFAOYSA-N")
              .replace("230.0552", "600.0"))
    # ピーク列を「86.0349（レコード1と共有）+ 300.0（レコード2固有）」に絞る。
    second = second.replace(
        """mbo:has_peak [ mbo:intensity 5.92827e+04 ;
            mbo:mz 8.60349e+01 ;
            mbo:peak_sequential_number 1 ;
            mbo:relative_intensity 11 ],
        [ mbo:intensity 1.271685e+06 ;
            mbo:mz 1.710069e+02 ;
            mbo:peak_sequential_number 2 ;
            mbo:relative_intensity 226 ],
        [ mbo:intensity 5.61971e+06 ;
            mbo:mz 2.300552e+02 ;
            mbo:peak_sequential_number 3 ;
            mbo:relative_intensity 999 ] ;
    mbo:has_peak_annotations [ mbo:mz 8.60349e+01 ;
            mbo:tentative_formula "C2H4N3O+" ],
        [ mbo:mz 1.710069e+02 ;
            mbo:tentative_formula "C5H4ClN4O+" ] ;""",
        """mbo:has_peak [ mbo:intensity 5.92827e+04 ;
            mbo:mz 8.60349e+01 ;
            mbo:peak_sequential_number 1 ;
            mbo:relative_intensity 11 ],
        [ mbo:intensity 1.271685e+06 ;
            mbo:mz 3.000000e+02 ;
            mbo:peak_sequential_number 2 ;
            mbo:relative_intensity 226 ] ;
    mbo:has_peak_annotations [ mbo:mz 8.60349e+01 ;
            mbo:tentative_formula "C2H4N3O+" ] ;""")
    second = second.replace("mbo:pk_num_peak 3 ;", "mbo:pk_num_peak 2 ;")
    (ttl_dir / "MSBNK-TEST-A0002.ttl").write_text(second, encoding="utf-8")

    build_index(ttl_dir, tmp_path / "index.sqlite")
    connection = sqlite3.connect(str(tmp_path / "index.sqlite"))
    yield connection
    connection.close()


def _query(peaks, precursor_mz=230.0552) -> QuerySpectrum:
    return QuerySpectrum(
        spectrum_id="q1", precursor_mz=precursor_mz, polarity="POSITIVE",
        peaks=tuple(Peak(mz, intensity, None, None) for mz, intensity in peaks),
        source="test")


@pytest.mark.parametrize("matched,expected", [
    (0, "PRECURSOR_ONLY"),
    (1, "WEAK"),
    (2, "PASS"),
    (7, "PASS"),
])
def test_derive_band(matched, expected):
    assert derive_band(matched) == expected


def test_exact_match_scores_high(conn):
    query = _query([(86.0349, 11.0), (171.0069, 226.0), (230.0552, 999.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert results[0].accession == "MSBNK-TEST-A0001"
    assert results[0].matched_fragments == 3
    assert results[0].band == "PASS"
    assert results[0].unexplained_peak_ratio == 0.0


def test_precursor_only_match_is_banded(conn):
    """precursor だけ合ってプロダクトイオンが合わない候補は PRECURSOR_ONLY。"""
    query = _query([(999.9999, 100.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert results[0].band == "PRECURSOR_ONLY"
    assert results[0].matched_fragments == 0
    assert results[0].unexplained_peak_ratio == 1.0


def test_precursor_mismatch_returns_nothing(conn):
    query = _query([(86.0349, 11.0)], precursor_mz=500.0)
    assert match_spectrum(conn, query, "LC-ESI-MS2") == []


def test_polarity_mismatch_returns_nothing(conn):
    query = QuerySpectrum(spectrum_id="q1", precursor_mz=230.0552,
                          polarity="NEGATIVE",
                          peaks=(Peak(86.0349, 11.0, None, None),),
                          source="test")
    assert match_spectrum(conn, query, "LC-ESI-MS2") == []


def test_boundary_peak_still_matches(conn):
    """ビン境界をまたぐ僅かなずれでも ±1 探索で拾える（Task 5）。"""
    query = _query([(86.0350, 11.0), (171.0070, 226.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert results[0].matched_fragments == 2


def test_unexplained_ratio_counts_query_peaks(conn):
    query = _query([(86.0349, 11.0), (400.0000, 50.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert results[0].matched_fragments == 1
    assert abs(results[0].unexplained_peak_ratio - 0.5) < 1e-9


def test_weighted_score_is_between_zero_and_one(conn):
    query = _query([(86.0349, 11.0), (400.0000, 50.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert 0.0 <= results[0].weighted_score <= 1.0


def test_top_k_limits_results(conn):
    query = _query([(86.0349, 11.0)])
    assert len(match_spectrum(conn, query, "LC-ESI-MS2", top_k=1)) <= 1


def test_results_sorted_by_weighted_score_only(conn):
    """並び順を決めるのは spectral_evidence だけ（spec §5.3 性質 1）。"""
    query = _query([(86.0349, 11.0), (171.0069, 226.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    scores = [r.weighted_score for r in results]
    assert scores == sorted(scores, reverse=True)


def test_weighted_score_is_positive_for_a_full_match(conn):
    """索引が 1 レコードだと全ビンの idf が log(1/1)=0 になり、主スコアが検証されない。
    2 レコードにして初めて重み付きカバレッジが意味を持つ。"""
    query = _query([(86.0349, 11.0), (171.0069, 226.0), (230.0552, 999.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert results[0].weighted_score > 0


def test_weighted_score_reflects_partial_coverage(conn):
    """クエリのピークのうち、索引が知っていて候補にも在るものの重み比が主スコア。

    171.0069 はレコード1に在る(idf=log2)。300.0 は索引は知っているが
    レコード1には無い(idf=log2)。よって log2 / (log2 + log2) = 0.5。
    """
    query = _query([(171.0069, 226.0), (300.0, 100.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert abs(results[0].weighted_score - 0.5) < 1e-9
