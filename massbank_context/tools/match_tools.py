"""照合ツール。Phase 1 ではバイアスを適用しない。"""
from __future__ import annotations

import sqlite3

from massbank_context.core import config
from massbank_context.core.errors import (
    incompatible_index,
    missing_index,
    missing_state,
)
from massbank_context.core.serialization import as_tsv, json_payload
from massbank_context.core.session import session
from massbank_context.index.writer import read_meta
from massbank_context.spectral.matcher import (
    DEFAULT_PRECURSOR_TOLERANCE,
    match_spectrum,
)

__all__ = ["match_spectra"]

_COLUMNS = ["spectrum_id", "accession", "inchikey", "matched_fragments",
            "total_query_peaks", "weighted_score", "cosine",
            "unexplained_peak_ratio", "band", "precursor_delta"]

_DEFAULT_ACQUISITION_CLASS = "LC-ESI-MS2"


def match_spectra(spectrum_ids: list[str] | None = None,
                  acquisition_class: str = _DEFAULT_ACQUISITION_CLASS,
                  precursor_tolerance: float = DEFAULT_PRECURSOR_TOLERANCE,
                  top_k: int = 10) -> str:
    """読み込んだスペクトルを索引と照合し、候補を TSV で返す。

    Phase 1 ではバイアス層が無いため bias.applied は常に false。
    objective 未設定はエラーにしない（spec §8）。
    """
    if not config.INDEX_PATH.exists():
        return missing_index(config.BUILD_COMMAND)
    if not session.spectra:
        return missing_state(
            "spectra", ["load_spectra"],
            "実験スペクトルが読み込まれていません。先に load_spectra を実行してください。")

    conn = sqlite3.connect(str(config.INDEX_PATH))
    try:
        meta = read_meta(conn)
        if meta is None:
            return missing_index(config.BUILD_COMMAND)
        available = meta.get("acquisition_classes", [])
        if acquisition_class not in available:
            return incompatible_index(
                acquisition_class, available,
                f"索引には {available} の層しかありません。測定原理の異なる層を"
                f"跨いだ照合は結果が無意味になるため行いません。")

        targets = (spectrum_ids if spectrum_ids is not None
                   else sorted(session.spectra))
        rows = []
        for spectrum_id in targets:
            spectrum = session.spectra.get(spectrum_id)
            if spectrum is None:
                continue
            for candidate in match_spectrum(
                    conn, spectrum, acquisition_class,
                    precursor_tolerance=precursor_tolerance, top_k=top_k):
                rows.append([
                    spectrum_id, candidate.accession, candidate.inchikey,
                    candidate.matched_fragments, candidate.total_query_peaks,
                    round(candidate.weighted_score, 4),
                    round(candidate.cosine, 4),
                    round(candidate.unexplained_peak_ratio, 4),
                    candidate.band,
                    None if candidate.precursor_delta is None
                    else round(candidate.precursor_delta, 4),
                ])
    finally:
        conn.close()

    return json_payload({
        "acquisition_class": acquisition_class,
        "n_spectra": len(targets),
        "n_candidates": len(rows),
        "candidates": as_tsv(_COLUMNS, rows),
        "bias": {
            "applied": False,
            "reason": "Phase 1 ではバイアス層が未実装。並び順は weighted_score のみで決まる。",
        },
        "note": "band が PRECURSOR_ONLY の候補は精密質量しか合っていない。"
                "PASS と同一視しないこと。",
    })
