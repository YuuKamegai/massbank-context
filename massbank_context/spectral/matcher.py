"""スペクトル照合（spec §6.8）。

主スコアは IDF 重み付きカバレッジ:
    weighted_score = Σ idf(f) for f in matched / Σ idf(f) for f in query_peaks

cosine ではなくカバレッジを主にする理由:
  1. unexplained_peak_ratio と直接対になり、「クエリのピークのうち何割が説明され、
     残りは何か」が一続きの説明になる。cosine の 0.81 は説明できない。
  2. cosine は強度比に敏感で、装置間・CE 間の差を拾う。カバレッジは頑健。

cosine も併記するが、並び順には使わない。
"""
from __future__ import annotations

import math
import sqlite3
from dataclasses import dataclass

from massbank_context.core.binning import PRODUCT_ION, probe_bins, to_bin
from massbank_context.core.models import QuerySpectrum
from massbank_context.index.background import lookup_idf
from massbank_context.index.writer import read_meta

BAND_PASS = "PASS"
BAND_WEAK = "WEAK"
BAND_PRECURSOR_ONLY = "PRECURSOR_ONLY"

DEFAULT_PRECURSOR_TOLERANCE = 0.01

#: 一致判定の窓をビン幅の何倍として表すか。クエリピークはまず中心ビン k に
#: 落ち、探索は k-1..k+1 を見る。その 3 ビンが覆う m/z 範囲は
#: [(k-1.5)*w, (k+1.5)*w) なので、ビン中心から見た窓はちょうど ±1.5*w になる。
FRAGMENT_TOLERANCE_BINS = 1.5


def fragment_tolerance(bin_width: float) -> float:
    """±1 ビン探索の実効トレランスを Da で返す。

    索引の bin_width だけを見せても「どれだけずれたピークまで一致とみなすのか」は
    読み取れない。結果を再現・検算できるよう、この値をツールの戻り値に載せる。
    """
    return bin_width * FRAGMENT_TOLERANCE_BINS


def _preference_order(bins: tuple[int, int, int]) -> tuple[int, int, int]:
    """probe_bins の (k-1, k, k+1) を (k, k-1, k+1) に並べ替える。

    一致を 1 対 1 に制約すると、どのビンを先に消費するかが結果を変える。
    中心ビン（クエリピークが実際に落ちたビン）を最優先にしないと、
    わずかにずれたピークが隣のビンを先取りし、そのビンに正対する別の
    ピークが一致を失う。探索する 3 ビンの集合そのものは変えない。
    """
    lower, centre, upper = bins
    return (centre, lower, upper)


@dataclass(frozen=True)
class MatchCandidate:
    """1 候補の spectral_evidence（spec §5.3）。

    メタデータは一切含まない。並び順を決めるのはこのブロックだけである。
    """
    accession: str
    inchikey: str
    matched_fragments: int
    total_query_peaks: int
    weighted_score: float
    cosine: float
    unexplained_peak_ratio: float
    band: str
    precursor_delta: float | None


def derive_band(matched_product_ions: int) -> str:
    """一致したプロダクトイオンの本数から band を決める（spec §5.3）。

    PRECURSOR_ONLY と PASS を同一視してはならない。精密質量が合っているだけの
    候補を、フラグメントが合った候補と同じ顔で並べないための区別である。
    """
    if matched_product_ions >= 2:
        return BAND_PASS
    if matched_product_ions == 1:
        return BAND_WEAK
    return BAND_PRECURSOR_ONLY


def _cosine(matched_pairs: list[tuple[float, float]],
            query_values: list[float],
            record_values: list[float]) -> float:
    """modified cosine（照合と同じ ±1 窓で整列させた版）。

    内積は「クエリのどのビンが record のどのビンに当たったか」という、
    カバレッジ照合が既に持っている対応関係のみを使う（境界を挟んだ対応も含む）。
    ノルムは対応の有無にかかわらず両者の全ピーク強度から取る。
    probe_bins（±1）で拾える一致を to_bin の完全一致だけで鍵付けすると、
    境界をまたいだピーク（カバレッジでは matched）が内積から抜け落ち、
    照合窓と cosine の窓がずれる。
    """
    if not matched_pairs:
        return 0.0
    dot = sum(q * r for q, r in matched_pairs)
    norm_q = math.sqrt(sum(v * v for v in query_values))
    norm_r = math.sqrt(sum(v * v for v in record_values))
    if norm_q == 0 or norm_r == 0:
        return 0.0
    return dot / (norm_q * norm_r)


def match_spectrum(conn: sqlite3.Connection,
                   spectrum: QuerySpectrum,
                   acquisition_class: str,
                   precursor_tolerance: float = DEFAULT_PRECURSOR_TOLERANCE,
                   top_k: int = 10) -> list[MatchCandidate]:
    """1 本のクエリスペクトルに対する候補を返す（weighted_score の降順）。"""
    meta = read_meta(conn)
    if meta is None:
        return []
    bin_width = meta["bin_width"]

    # --- 手順 2: precursor m/z と polarity で候補レコードを絞る ---
    if spectrum.precursor_mz is None:
        rows = conn.execute(
            "SELECT record_id, accession, inchikey, precursor_mz FROM records"
            " WHERE acquisition_class=? AND (polarity=? OR ? IS NULL)",
            (acquisition_class, spectrum.polarity, spectrum.polarity)).fetchall()
    else:
        rows = conn.execute(
            "SELECT record_id, accession, inchikey, precursor_mz FROM records"
            " WHERE acquisition_class=? AND (polarity=? OR ? IS NULL)"
            "   AND precursor_mz BETWEEN ? AND ?",
            (acquisition_class, spectrum.polarity, spectrum.polarity,
             spectrum.precursor_mz - precursor_tolerance,
             spectrum.precursor_mz + precursor_tolerance)).fetchall()
    if not rows:
        return []

    # --- クエリ側のビンと idf を先に用意する ---
    # 一致判定は probe_bins（±1）で探る。境界を挟んだ相手を取りこぼさないため。
    # 一方 idf は中心ビン 1 つだけを引く。背景表は既に ±1 展開済みで、
    # 中心ビンの行が「このクエリピークが ±1 探索で当たる InChIKey 数」だから
    # である（詳細は background.lookup_idf の docstring）。窓を 2 つの目的で
    # 使い分けているので、片方に合わせて統一しないこと。
    query_bins = [probe_bins(peak.mz, bin_width) for peak in spectrum.peaks]
    query_idf = [lookup_idf(conn, acquisition_class, spectrum.polarity,
                            PRODUCT_ION, to_bin(peak.mz, bin_width))
                 for peak in spectrum.peaks]
    query_values = [peak.relative_intensity for peak in spectrum.peaks]
    total_idf = sum(query_idf)
    total_peaks = len(spectrum.peaks)

    candidates: list[MatchCandidate] = []
    for record_id, accession, inchikey, record_precursor in rows:
        # 同一ビンに複数の record ピークが落ちたら強度の最大値を代表にする。
        # dict 内包で後勝ちにすると代表値が任意になり、しかも cosine の
        # ノルムから消えたピークの分だけ norm_r が過小になる。
        record_bins: dict[int, float] = {}
        for bin_no, intensity in conn.execute(
                "SELECT bin, relative_intensity FROM peaks"
                " WHERE record_id=? AND axis=?", (record_id, PRODUCT_ION)):
            value = intensity or 0.0
            if value > record_bins.get(bin_no, float("-inf")):
                record_bins[bin_no] = value
        matched_idf = 0.0
        matched = 0
        cosine_pairs: list[tuple[float, float]] = []
        # 一致は 1 対 1 に制約する。record の 1 本を複数のクエリピークで
        # 使い回すと、(a) matched_fragments が化学的に別のフラグメントを
        # 重複計上し、weighted_score（並び順を決める主スコア）が上振れし、
        # (b) 同じ強度が内積で複数回加算される一方 norm_r には 1 回しか
        # 現れないため cosine が 1 を超える。実データ（bin_width 0.05）で
        # 実際に 1.04 が出た。消費済みのビンは以降のピークに渡さない。
        consumed: set[int] = set()
        for peak, bins, idf in zip(spectrum.peaks, query_bins, query_idf):
            # probe_bins（±1）の中で record が実際に持つビンを探す。カバレッジの
            # matched 判定と同じ窓を使うことで、cosine が別の一致基準を持たない
            # ようにする。中心ビンを最優先で見るのは、1 対 1 にすると割り当ての
            # 順序が結果を変えるため。ずれの小さいピークに隣のビンを先取りさせると、
            # そのビンに正対するピークが弾かれる。
            hit_bin = next(
                (b for b in _preference_order(bins)
                 if b in record_bins and b not in consumed), None)
            if hit_bin is not None:
                consumed.add(hit_bin)
                matched += 1
                matched_idf += idf
                cosine_pairs.append((peak.relative_intensity, record_bins[hit_bin]))

        candidates.append(MatchCandidate(
            accession=accession,
            inchikey=inchikey,
            matched_fragments=matched,
            total_query_peaks=total_peaks,
            weighted_score=(matched_idf / total_idf) if total_idf > 0 else 0.0,
            cosine=_cosine(cosine_pairs, query_values, list(record_bins.values())),
            unexplained_peak_ratio=(
                1.0 - matched / total_peaks) if total_peaks else 1.0,
            band=derive_band(matched),
            precursor_delta=(
                spectrum.precursor_mz - record_precursor)
            if spectrum.precursor_mz is not None and record_precursor is not None
            else None,
        ))

    # 並び順を決めるのは weighted_score だけ。同点は accession で決定的にする。
    candidates.sort(key=lambda c: (-c.weighted_score, c.accession))
    return candidates[:top_k]
