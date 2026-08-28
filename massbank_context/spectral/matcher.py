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


def _cosine(query_intensities: dict[int, float],
            record_intensities: dict[int, float]) -> float:
    shared = set(query_intensities) & set(record_intensities)
    if not shared:
        return 0.0
    dot = sum(query_intensities[b] * record_intensities[b] for b in shared)
    norm_q = math.sqrt(sum(v * v for v in query_intensities.values()))
    norm_r = math.sqrt(sum(v * v for v in record_intensities.values()))
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
    query_bins = [probe_bins(peak.mz, bin_width) for peak in spectrum.peaks]
    query_idf = [lookup_idf(conn, acquisition_class, PRODUCT_ION, bins)
                 for bins in query_bins]
    query_intensities = {to_bin(peak.mz, bin_width): peak.relative_intensity
                         for peak in spectrum.peaks}
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
        for bins, idf in zip(query_bins, query_idf):
            if any(b in record_bins for b in bins):
                matched += 1
                matched_idf += idf

        candidates.append(MatchCandidate(
            accession=accession,
            inchikey=inchikey,
            matched_fragments=matched,
            total_query_peaks=total_peaks,
            weighted_score=(matched_idf / total_idf) if total_idf > 0 else 0.0,
            cosine=_cosine(query_intensities, record_bins),
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
