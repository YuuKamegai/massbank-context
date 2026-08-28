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
    query_idf = [lookup_idf(conn, acquisition_class, PRODUCT_ION,
                            to_bin(peak.mz, bin_width))
                 for peak in spectrum.peaks]
    query_values = [peak.relative_intensity for peak in spectrum.peaks]
    total_idf = sum(query_idf)
    total_peaks = len(spectrum.peaks)

    candidates: list[MatchCandidate] = []
    for record_id, accession, inchikey, record_precursor in rows:
        record_bins = {
            bin_no: intensity
            for bin_no, intensity in conn.execute(
                "SELECT bin, relative_intensity FROM peaks"
                " WHERE record_id=? AND axis=?", (record_id, PRODUCT_ION))
        }
        matched_idf = 0.0
        matched = 0
        cosine_pairs: list[tuple[float, float]] = []
        for peak, bins, idf in zip(spectrum.peaks, query_bins, query_idf):
            # probe_bins（±1）の中で record が実際に持つビンを探す。カバレッジの
            # matched 判定と同じ窓を使うことで、cosine が別の一致基準を持たない
            # ようにする。複数当たっても二重計上しないよう最初の 1 つだけ使う。
            hit_bin = next((b for b in bins if b in record_bins), None)
            if hit_bin is not None:
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
