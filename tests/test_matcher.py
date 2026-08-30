import math
import shutil
import sqlite3
from pathlib import Path

import pytest

from massbank_context.core.models import Peak, QuerySpectrum
from massbank_context.index.__main__ import build_index
from massbank_context.spectral.matcher import derive_band, match_spectrum

FIXTURES = Path(__file__).parent / "fixtures" / "ttl"

# gen_a の構造を保ったまま、accession・InChIKey・precursor・ピーク列だけを
# 差し替えた合成レコードを作るテンプレート。3 件目以降のレコードを増やすたびに
# 巨大な文字列置換を書き足すのを避けるため、多候補用の conn_multi フィクスチャで使う。
_TTL_TEMPLATE = """@prefix dcterms: <http://purl.org/dc/terms/> .
@prefix inchikey: <http://identifiers.org/inchikey/> .
@prefix massbank: <https://identifiers.org/massbank/> .
@prefix mbo: <http://www.massbank.jp/ontology/> .
@prefix obo: <http://purl.obolibrary.org/obo/> .
@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
@prefix schema: <https://schema.org/> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .

massbank:{accession} a mbo:MassSpectrum ;
    dcterms:identifier "{accession}" ;
    dcterms:license <https://creativecommons.org/publicdomain/zero/1.0/> ;
    mbo:analytical_methods_and_conditions [ a mbo:AnalyticalMethods ;
            mbo:ac_instrument "LTQ Orbitrap XL Thermo Scientific" ;
            mbo:collision_energy "35% (nominal)" ;
            mbo:fragmentation_mode "CID" ;
            mbo:instrument_type "LC-ESI-ITFT" ;
            mbo:ion_mode "{ion_mode}" ;
            mbo:ms_type "MS2" ;
            mbo:resolution "15000" ] ;
    mbo:compound [ a schema:chemicalsubstance ;
            rdfs:label "Test{accession}" ;
            mbo:ch_compound_class "N/A; Environmental Standard" ;
            mbo:ch_exact_mass 2.290479e+02 ;
            schema:inChIKey inchikey:{inchikey} ;
            schema:molecularFormula "C6H8ClN7O" ;
            schema:name "Test{accession}" ] ;
    {peak_block}
    mbo:mass_spectram_meta_data [ a mbo:MassSpectramMetaData ;
            mbo:base_peak 2.300551e+02 ;
            mbo:precursor_mz_value "{precursor_mz}"^^xsd:float ;
            mbo:precursor_type_value "[M+H]+" ] ;
    mbo:pk_num_peak {n_peaks} ;
    mbo:pk_splash "splash10-001i-0290000000-000000000000" .
"""


def _peak_block(peaks: list[tuple[float, float]]) -> str:
    """[(mz, relative_intensity), ...] から mbo:has_peak ブロックを作る。

    mbo:intensity は reader_gen_a が読まない列なので relative_intensity と
    同じ値を入れておけば足りる。has_peak_annotations（tentative_formula）は
    matcher が使わないため省略する（世代Aでもアノテーション無しレコードは
    実在するので、これ自体が読み込み側の未検証経路にはならない）。
    """
    entries = ",\n        ".join(
        f"[ mbo:intensity {intensity:.6e} ;\n"
        f"            mbo:mz {mz:.6e} ;\n"
        f"            mbo:peak_sequential_number {i} ;\n"
        f"            mbo:relative_intensity {intensity:.0f} ]"
        for i, (mz, intensity) in enumerate(peaks, start=1)
    )
    return f"mbo:has_peak {entries} ;"


def _write_ttl_record(ttl_dir: Path, accession: str, inchikey: str,
                      precursor_mz: float,
                      peaks: list[tuple[float, float]],
                      ion_mode: str = "POSITIVE") -> None:
    text = _TTL_TEMPLATE.format(
        accession=accession, inchikey=inchikey, precursor_mz=precursor_mz,
        peak_block=_peak_block(peaks), n_peaks=len(peaks), ion_mode=ion_mode)
    (ttl_dir / f"{accession}.ttl").write_text(text, encoding="utf-8")


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


@pytest.fixture()
def conn_multi(tmp_path):
    """並び順・top_k の検証専用。候補が常に 1 件しかない conn とは別の索引にする。

    conn の全テストは precursor 230.0552 だけを使い、conn の 2 件目
    (precursor 600.0) は絶対に候補に入らないため、conn では並び順の比較子
    （weighted_score 降順、accession でタイブレーク）も top_k の切り詰めも
    1 要素のリストにしか触れず、実際には何も検証していなかった。
    ここでは precursor 230.0552 の窓に 2 件（レコード1・レコード3）を落とし、
    互いにスコアが異なるようにして初めて比較子を動かす。
    """
    ttl_dir = tmp_path / "ttl"
    ttl_dir.mkdir()
    shutil.copy(FIXTURES / "gen_a" / "MSBNK-TEST-A0001.ttl", ttl_dir)
    # レコード2: precursor 600.0（候補窓の外）。IDF 背景に厚みを持たせるためだけに置く。
    _write_ttl_record(ttl_dir, "MSBNK-TEST-A0002", "YYYYYYYYYYYYYY-UHFFFAOYSA-N",
                      600.0, [(86.0349, 11.0), (300.0, 226.0)])
    # レコード3: precursor 230.0552（レコード1と同じ窓）だが、171.0069 までしか
    # 持たず 230.0552 を欠く。レコード1（フルマッチ）より必ず低いスコアになる。
    _write_ttl_record(ttl_dir, "MSBNK-TEST-A0003", "ZZZZZZZZZZZZZZ-UHFFFAOYSA-N",
                      230.0552, [(86.0349, 11.0), (171.0069, 226.0)])

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
    """ビン境界をまたぐ僅かなずれでも ±1 探索で拾える（Task 5）。

    cosine もカバレッジと同じ ±1 窓で対応を取る。ここでは両ピークとも境界を
    またいで一致するので、cosine の内積が 0 のまま（対応漏れ）になっていない
    ことを検証する。record 側の強度ベクトルは (11, 226, 999)、クエリ側は
    一致した 2 本 (11, 226) だけが対応するので、
        cosine = (11*11 + 226*226) / (sqrt(11^2+226^2) * sqrt(11^2+226^2+999^2))
               = 51197 / (sqrt(51197) * sqrt(1049198))
    """
    query = _query([(86.0350, 11.0), (171.0070, 226.0)])
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert results[0].matched_fragments == 2
    expected_cosine = 51197 / ((51197 ** 0.5) * (1049198 ** 0.5))
    assert abs(results[0].cosine - expected_cosine) < 1e-9


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


def test_orders_multiple_candidates_by_weighted_score(conn_multi):
    """候補が 2 件以上あるときに初めて比較子が働く。1 件では順序を検証できない。"""
    query = _query([(86.0349, 11.0), (171.0069, 226.0), (230.0552, 999.0)])
    results = match_spectrum(conn_multi, query, "LC-ESI-MS2")
    assert len(results) >= 2
    scores = [r.weighted_score for r in results]
    assert scores == sorted(scores, reverse=True)
    assert scores[0] > scores[-1]      # 同点でないこと＝比較子が実際に働いたこと


def test_top_k_limits_results_with_multiple_candidates(conn_multi):
    """候補が 2 件以上あるときに top_k の切り詰めを検証する。"""
    query = _query([(86.0349, 11.0), (171.0069, 226.0), (230.0552, 999.0)])
    all_results = match_spectrum(conn_multi, query, "LC-ESI-MS2")
    assert len(all_results) == 2
    limited = match_spectrum(conn_multi, query, "LC-ESI-MS2", top_k=1)
    assert len(limited) == 1
    assert limited[0].accession == all_results[0].accession


def test_query_without_precursor_returns_all_records_in_class(conn):
    """precursor が無いクエリは precursor で絞れない。曖昧なまま拒否せず全件を候補にする。

    `>= 1` では、絞り込みが部分的に復活して通常の窓から外れるレコード
    （ここでは precursor 600.0 の MSBNK-TEST-A0002）が落ちる退行を見逃す
    （2 件 → 1 件に減っても通ってしまう）。conn フィクスチャの中身は
    決まっているので、実測した accession の並びをそのまま固定する。
    両者とも matched_fragments=1 で weighted_score が 0.0 の同点になるため
    （bin 8603 は 2 レコード両方に出現し idf=log(2/2)=0）、並び順は
    accession のタイブレークで決まり、実行のたびに安定する。
    """
    query = QuerySpectrum(spectrum_id="q", precursor_mz=None, polarity="POSITIVE",
                          peaks=(Peak(86.0349, 11.0, None, None),), source="test")
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert [r.accession for r in results] == [
        "MSBNK-TEST-A0001", "MSBNK-TEST-A0002"]
    assert all(r.precursor_delta is None for r in results)


def test_query_without_precursor_or_polarity_does_not_filter_by_polarity(conn):
    """polarity=None は「極性で絞らない」を意味する（.dcl の ion_mode=2 が正当に生む）。

    上のテストと同じ理由で件数・同一性を固定する。conn フィクスチャは
    2 レコードとも polarity=POSITIVE で登録されているため、polarity で
    絞らなくても結果は同じ 2 件になる。ここで検証したいのは
    「precursor で絞られない」という経路そのものであり、この 2 レコードの
    組み合わせでは polarity 無絞り込みの効果は結果件数には現れない。
    """
    query = QuerySpectrum(spectrum_id="q", precursor_mz=None, polarity=None,
                          peaks=(Peak(86.0349, 11.0, None, None),), source="test")
    results = match_spectrum(conn, query, "LC-ESI-MS2")
    assert [r.accession for r in results] == [
        "MSBNK-TEST-A0001", "MSBNK-TEST-A0002"]
    assert all(r.precursor_delta is None for r in results)


# --- 極性混在の索引（finding 3 / 4 / 5） -------------------------------------


@pytest.fixture()
def conn_mixed(tmp_path):
    """POSITIVE と NEGATIVE が混ざった索引。

    他のフィクスチャは全レコードが POSITIVE なので、
      - POSITIVE のクエリが NEGATIVE のレコードを拾わないこと
      - polarity=None のクエリが両極性を拾うこと
      - polarity=None の重みが極性プール "*" の行から来ていること
    のいずれも実際には検証できない（結果が同じになるため）。

    配置（bin_width は実測できず fallback 0.01）:
      A0001 POSITIVE precursor 230.0552  peaks 86.0349 / 171.0069 / 230.0552
      A0004 NEGATIVE precursor 230.0552  peaks 86.0349 / 300.0
      A0002 POSITIVE precursor 600.0     peaks 86.0349 / 300.0   （背景の厚み用）
    """
    ttl_dir = tmp_path / "ttl"
    ttl_dir.mkdir()
    shutil.copy(FIXTURES / "gen_a" / "MSBNK-TEST-A0001.ttl", ttl_dir)
    _write_ttl_record(ttl_dir, "MSBNK-TEST-A0002", "YYYYYYYYYYYYYY-UHFFFAOYSA-N",
                      600.0, [(86.0349, 11.0), (300.0, 226.0)])
    _write_ttl_record(ttl_dir, "MSBNK-TEST-A0004", "NNNNNNNNNNNNNN-UHFFFAOYSA-N",
                      230.0552, [(86.0349, 11.0), (300.0, 226.0)],
                      ion_mode="NEGATIVE")

    build_index(ttl_dir, tmp_path / "index.sqlite")
    connection = sqlite3.connect(str(tmp_path / "index.sqlite"))
    yield connection
    connection.close()


def test_positive_query_does_not_return_negative_records(conn_mixed):
    """極性は候補プールの絞り込み条件。POSITIVE のクエリに NEGATIVE は出ない。"""
    query = _query([(86.0349, 11.0)])
    results = match_spectrum(conn_mixed, query, "LC-ESI-MS2")
    assert [r.accession for r in results] == ["MSBNK-TEST-A0001"]


def test_polarity_none_query_returns_both_polarities(conn_mixed):
    """polarity=None は「極性で絞らない」。両極性のレコードが候補になる。"""
    query = QuerySpectrum(spectrum_id="q", precursor_mz=230.0552, polarity=None,
                          peaks=(Peak(86.0349, 11.0, None, None),),
                          source="test")
    results = match_spectrum(conn_mixed, query, "LC-ESI-MS2")
    assert {r.accession for r in results} == {
        "MSBNK-TEST-A0001", "MSBNK-TEST-A0004"}
    polarities = {
        row[0] for row in conn_mixed.execute(
            "SELECT polarity FROM records WHERE accession IN (?,?)",
            tuple(r.accession for r in results))}
    assert polarities == {"POSITIVE", "NEGATIVE"}


def test_polarity_none_weights_come_from_the_pooled_background(conn_mixed):
    """極性不明クエリの重みは、プール "*"（全極性）の行から来なければならない。

    候補プールが両極性にまたがるのに背景が片極性だと、idf が「いま探している
    集団の中での稀さ」を意味しなくなる。

    クエリのピークは 171.0069（bin 17101）と 300.0（bin 30000）。
      プール "*" (N=3): bin 17101 -> A0001 のみ(df=1, idf=log3)
                        bin 30000 -> A0002 と A0004(df=2, idf=log(3/2))
      POSITIVE プール (N=2): 両ビンとも df=1 -> idf=log2 で等しく、
                        weighted_score はちょうど 0.5 になってしまう。
    A0001 は 171.0069 だけ一致するので、正しい値は log3/(log3+log(3/2))。
    """
    query = QuerySpectrum(spectrum_id="q", precursor_mz=230.0552, polarity=None,
                          peaks=(Peak(171.0069, 226.0, None, None),
                                 Peak(300.0, 100.0, None, None)),
                          source="test")
    results = match_spectrum(conn_mixed, query, "LC-ESI-MS2")
    top = next(r for r in results if r.accession == "MSBNK-TEST-A0001")
    expected = math.log(3) / (math.log(3) + math.log(3 / 2))
    assert abs(top.weighted_score - expected) < 1e-9
    assert abs(top.weighted_score - 0.5) > 0.1      # 片極性の背景ではない


def test_positive_query_weights_come_from_the_positive_pool(conn_mixed):
    """極性を指定したクエリの重みは、その極性のプールだけで数えた行から来る。

    上のテストと同じ 2 ピークのクエリを POSITIVE で投げる。候補は A0001 だけ。
      POSITIVE プール (N=2): bin 17101 -> A0001 のみ, bin 30000 -> A0002 のみ。
      どちらも df=1 で idf=log2 なので、171.0069 だけ一致する A0001 は
      ちょうど 0.5 になる。
    層別しない（＝極性混合 N=3）背景を引くと log3/(log3+log(3/2))≒0.73 になり、
    この値にはならない。
    """
    query = _query([(171.0069, 226.0), (300.0, 100.0)])
    results = match_spectrum(conn_mixed, query, "LC-ESI-MS2")
    assert [r.accession for r in results] == ["MSBNK-TEST-A0001"]
    assert abs(results[0].weighted_score - 0.5) < 1e-9


@pytest.fixture()
def conn_wide_bin(tmp_path):
    """ビンが広く、複数のクエリピークが 1 本の record ピークに重なる索引。

    実データ（bin_width 0.05）で cosine が 1.0411 になった経路を再現する。
    レコードのピークを離して置き、クエリ側だけを 1 ビンに詰める。
    """
    ttl_dir = tmp_path / "ttl"
    ttl_dir.mkdir()
    _write_ttl_record(ttl_dir, "MSBNK-TEST-W0001", "WWWWWWWWWWWWWW-UHFFFAOYSA-N",
                      230.0552, [(100.0, 999.0), (200.0, 500.0)])
    _write_ttl_record(ttl_dir, "MSBNK-TEST-W0002", "VVVVVVVVVVVVVV-UHFFFAOYSA-N",
                      600.0, [(100.0, 999.0), (300.0, 500.0)])
    build_index(ttl_dir, tmp_path / "index.sqlite")
    connection = sqlite3.connect(str(tmp_path / "index.sqlite"))
    yield connection
    connection.close()


def _crowded_query(bin_width: float) -> QuerySpectrum:
    """1 本の record ピーク(100.0) の探索窓に 4 本のクエリピークを詰めた質問。"""
    step = bin_width / 2
    peaks = [(100.0 + i * step, 999.0) for i in range(4)]
    return _query(peaks)


def test_cosine_never_exceeds_one(conn_wide_bin):
    """cosine は定義上 1 を超えない。

    多対一の一致を許すと、record の強度が内積で複数回加算される一方
    ノルムには 1 回しか現れず、Cauchy-Schwarz が破れる。
    """
    from massbank_context.index.writer import read_meta
    bin_width = read_meta(conn_wide_bin)["bin_width"]
    results = match_spectrum(conn_wide_bin, _crowded_query(bin_width),
                             "LC-ESI-MS2")
    assert results
    for candidate in results:
        assert candidate.cosine <= 1.0 + 1e-9, candidate


def test_matched_fragments_cannot_exceed_record_peaks(conn_wide_bin):
    """一致は 1 対 1。record の 1 本を複数のクエリピークで数え上げない。

    matched_fragments が水増しされると weighted_score（並び順を決める主スコア）も
    同時に上振れするため、これは表示上の問題ではない。
    """
    from massbank_context.core.binning import PRODUCT_ION
    from massbank_context.index.writer import read_meta
    bin_width = read_meta(conn_wide_bin)["bin_width"]
    results = match_spectrum(conn_wide_bin, _crowded_query(bin_width),
                             "LC-ESI-MS2")
    assert results
    for candidate in results:
        record_peaks = conn_wide_bin.execute(
            "SELECT COUNT(*) FROM peaks p JOIN records r"
            " ON r.record_id = p.record_id"
            " WHERE r.accession=? AND p.axis=?",
            (candidate.accession, PRODUCT_ION)).fetchone()[0]
        assert candidate.matched_fragments <= record_peaks, candidate


def test_exact_bin_is_preferred_over_neighbour(conn_wide_bin):
    """中心ビンに当たる相手が居るなら、隣のビンを先に消費しない。

    一致を 1 対 1 にすると割り当ての順序が結果を変える。中心を優先しないと、
    ずれの小さいピークが隣のビンを取ってしまい、本来当たるはずの
    ピークが弾かれる。
    """
    from massbank_context.index.writer import read_meta
    bin_width = read_meta(conn_wide_bin)["bin_width"]
    # 100.0 ちょうど（中心ビン）と、その隣のビンに落ちる 1 本。
    query = _query([(100.0, 999.0), (100.0 + bin_width, 500.0)])
    results = match_spectrum(conn_wide_bin, query, "LC-ESI-MS2")
    top = next(c for c in results if c.accession == "MSBNK-TEST-W0001")
    assert top.matched_fragments == 1
