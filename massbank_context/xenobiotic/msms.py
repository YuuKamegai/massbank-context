"""MS/MS の筋: 索引の読み出し、参照スペクトルとのビンの共有、順向きの順位。

共有の判定は spectral/matcher.py と同じ規則（±1 ビンを探り、中心ビンを優先して
参照のビンを一度だけ使う）を、import せずにここへ書く。
"""
from __future__ import annotations

import sqlite3
from collections.abc import Sequence
from dataclasses import dataclass

from massbank_context.core import config
from massbank_context.core.binning import NEUTRAL_LOSS, PRODUCT_ION, neighbour_bins, to_bin
from massbank_context.core.models import QuerySpectrum
from massbank_context.index.background import lookup_idf
from massbank_context.index.writer import open_index_readonly, read_meta
from massbank_context.spectral.matcher import (
    BAND_PASS, BAND_WEAK, DEFAULT_PRECURSOR_TOLERANCE, match_spectrum_result)

__all__ = ["ACQUISITION_CLASS", "MSMS_ACQUISITION_CLASS", "PRECURSOR_MARGIN", "FORWARD_TOP_K",
           "RefBin", "ReferenceBins", "SharedResult", "ForwardHit", "open_index", "reference_bins",
           "query_axis_bins", "shared_evidence", "forward_rank"]

#: 外来の分子の走査が照合する層。測定原理の異なる層は跨がない（spec §6.1）。
ACQUISITION_CLASS = "LC-ESI-MS2"
MSMS_ACQUISITION_CLASS = ACQUISITION_CLASS
#: 前駆体の残り（前駆体 − この値以上の product ion）は共有に数えない。
PRECURSOR_MARGIN = 0.5
#: 順向きの順位を見る候補の件数（match_spectrum_result の上限）。
FORWARD_TOP_K = 20


def open_index() -> tuple[sqlite3.Connection, dict]:
    """索引を読み取り専用で開き、(接続, index_meta) を返す。

    索引が無い・読めない場合は ``missing_index``、LC-ESI-MS2 の層が無い場合は
    ``incompatible_index`` を ValueError で返す。別の層で黙って代用しない。
    """
    path = config.INDEX_PATH
    if not path.exists():
        raise ValueError(f"missing_index: 索引がありません: {path}。{config.BUILD_COMMAND}")
    conn = open_index_readonly(path)
    meta = read_meta(conn)
    if meta is None:
        conn.close()
        raise ValueError(f"missing_index: 索引が空または壊れています: {path}。{config.BUILD_COMMAND}")
    if ACQUISITION_CLASS not in meta.get("acquisition_classes", []):
        conn.close()
        raise ValueError(f"incompatible_index: 索引に {ACQUISITION_CLASS} の層がありません"
                         f"（acquisition_classes={meta.get('acquisition_classes')}）")
    return conn, meta


@dataclass(frozen=True)
class RefBin:
    """参照の 1 ビン。idf は中心ビン 1 行、強度は全レコード中の最大。"""
    idf: float
    intensity: float
    formula: str | None


@dataclass(frozen=True)
class ReferenceBins:
    """登録した分子の参照スペクトル（同一 InChIKey の全レコードを併せたもの）。"""
    inchikey: str
    n_records: int
    bins: dict[tuple[str, int], RefBin]
    top: frozenset[tuple[str, int]]


@dataclass(frozen=True)
class SharedResult:
    """得点（0〜1）、共有したビン数（和集合全体）、共有したビンの組成式。"""
    score: float
    n_shared: int
    formulas: tuple[str, ...]


def reference_bins(conn: sqlite3.Connection, inchikey: str, index_polarity: str,
                   bin_width: float, top_n: int) -> ReferenceBins | None:
    """InChIKey の参照ビンを作る。LC-ESI-MS2 のレコードが無ければ None。

    product ion のうち m/z が前駆体 − PRECURSOR_MARGIN 以上のものは除く。
    前駆体が NULL のレコードは除外せず残す。idf は中心ビン 1 行だけを引く
    （±1 に広げない。index/background.py の lookup_idf の docstring を参照）。
    """
    records = conn.execute(
        "SELECT record_id, precursor_mz FROM records"
        " WHERE inchikey=? AND acquisition_class=? AND polarity=?",
        (inchikey, ACQUISITION_CLASS, index_polarity)).fetchall()
    if not records:
        return None
    intensity: dict[tuple[str, int], float] = {}
    formula: dict[tuple[str, int], str | None] = {}
    for record_id, precursor in sorted(records):
        rows = conn.execute(
            "SELECT axis, bin, relative_intensity, formula FROM peaks WHERE record_id=?"
            " ORDER BY axis, bin", (record_id,)).fetchall()
        for axis, bin_no, value, name in rows:
            if (axis == PRODUCT_ION and precursor is not None
                    and bin_no * bin_width >= precursor - PRECURSOR_MARGIN):
                continue
            key = (axis, bin_no)
            intensity[key] = max(intensity.get(key, 0.0), value or 0.0)
            if not formula.get(key) and name:
                formula[key] = name
            else:
                formula.setdefault(key, None)
    bins = {key: RefBin(lookup_idf(conn, ACQUISITION_CLASS, index_polarity, key[0], key[1]),
                        intensity[key], formula[key]) for key in intensity}
    top: set[tuple[str, int]] = set()
    for axis in (PRODUCT_ION, NEUTRAL_LOSS):
        ranked = sorted((k for k in bins if k[0] == axis), key=lambda k: (-bins[k].intensity, k[1]))
        top.update(ranked[:top_n])
    return ReferenceBins(inchikey, len(records), bins, frozenset(top))


def query_axis_bins(spectrum: QuerySpectrum, bin_width: float) -> list[tuple[str, int]]:
    """クエリスペクトルの (軸, ビン) の列。

    product ion は前駆体 − PRECURSOR_MARGIN 未満だけ。neutral loss は前駆体 − mz が
    PRECURSOR_MARGIN 以上のものだけで、前駆体が無ければ作らない。
    """
    precursor = spectrum.precursor_mz
    result: list[tuple[str, int]] = []
    for peak in spectrum.peaks:
        if precursor is None or peak.mz < precursor - PRECURSOR_MARGIN:
            result.append((PRODUCT_ION, to_bin(peak.mz, bin_width)))
        if precursor is not None and precursor - peak.mz >= PRECURSOR_MARGIN:
            result.append((NEUTRAL_LOSS, to_bin(precursor - peak.mz, bin_width)))
    return result


def _shared_for(ref: ReferenceBins, spectrum: QuerySpectrum, bin_width: float) -> SharedResult:
    consumed: set[tuple[str, int]] = set()
    for axis, bin_no in query_axis_bins(spectrum, bin_width):
        lower, centre, upper = neighbour_bins(bin_no)
        # 中心ビンを優先し、参照のビンは一度しか使わない（matcher と同じ規則）。
        hit = next(((axis, b) for b in (centre, lower, upper)
                    if (axis, b) in ref.bins and (axis, b) not in consumed), None)
        if hit is not None:
            consumed.add(hit)
    denominator = sum(ref.bins[k].idf for k in ref.top)
    numerator = sum(ref.bins[k].idf for k in consumed if k in ref.top)
    formulas = sorted({ref.bins[k].formula for k in consumed if ref.bins[k].formula})
    return SharedResult(numerator / denominator if denominator > 0 else 0.0,
                        len(consumed), tuple(formulas))


def shared_evidence(ref: ReferenceBins, spectra: Sequence[QuerySpectrum],
                    bin_width: float) -> SharedResult | None:
    """参照とのビンの共有。全スペクトルがピーク 0 本なら None、複数なら得点最大のもの。"""
    results = [_shared_for(ref, s, bin_width) for s in spectra if s.peaks]
    if not results:
        return None
    return max(results, key=lambda r: r.score)


@dataclass(frozen=True)
class ForwardHit:
    """順向きの照合の順位（1 始まり）と、その順位を決めた候補の band。"""
    rank: int
    band: str


def _forward_order(hit: ForwardHit) -> tuple[int, int]:
    strength = {BAND_PASS: 0, BAND_WEAK: 1}.get(hit.band, 2)
    return hit.rank, strength


def forward_rank(conn: sqlite3.Connection, spectra: Sequence[QuerySpectrum],
                 inchikey: str) -> ForwardHit | None:
    """特徴のスペクトルを照合した候補の中での InChIKey の順位と band。

    重複を除いた InChIKey の順に並べ、band はその InChIKey の先頭の候補のものを取る。
    全候補が得点 0 なら順位は accession の文字列順で決まるので、順位だけでは同定の根拠に
    ならない。複数のスペクトルがあれば順位が最小のもの、同順位なら band の強いものを返す。
    上位 FORWARD_TOP_K に無ければ None。
    """
    best: ForwardHit | None = None
    for spectrum in spectra:
        if not spectrum.peaks:
            continue
        result = match_spectrum_result(conn, spectrum, ACQUISITION_CLASS,
                                       DEFAULT_PRECURSOR_TOLERANCE, FORWARD_TOP_K)
        first = {}
        for candidate in result.candidates:
            first.setdefault(candidate.inchikey, candidate)
        if inchikey in first:
            hit = ForwardHit(list(first).index(inchikey) + 1, first[inchikey].band)
            best = hit if best is None else min(best, hit, key=_forward_order)
    return best
