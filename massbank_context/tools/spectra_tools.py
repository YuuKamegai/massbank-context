"""実験スペクトルの読み込みと確認。"""
from __future__ import annotations

from collections import Counter

from massbank_context.core.errors import missing_state
from massbank_context.core.serialization import MZ_DIGITS, json_payload, round_floats
from massbank_context.core.session import session
from massbank_context.query.loader import load_spectra_file

__all__ = ["load_spectra", "inspect_spectrum"]


def load_spectra(path: str) -> str:
    """スペクトルを読み込み、要約だけを返す。ピーク列はセッションに残す。"""
    try:
        spectra = load_spectra_file(path)
    except (OSError, ValueError) as exc:
        return json_payload({"error": {"code": "load_failed",
                                       "message": str(exc)}})
    session.set_spectra(spectra, path)

    precursors = [s.precursor_mz for s in spectra if s.precursor_mz is not None]
    return json_payload(round_floats({
        "source": path,
        "n_spectra": len(spectra),
        "polarity_counts": dict(Counter(
            s.polarity for s in spectra if s.polarity)),
        "precursor_mz_range": [min(precursors), max(precursors)]
        if precursors else None,
        "spectrum_ids": [s.spectrum_id for s in spectra],
    }, MZ_DIGITS))


def inspect_spectrum(spectrum_id: str) -> str:
    """1 本のピーク列を丸めて返す。必要なときだけ呼ぶこと。"""
    if not session.spectra:
        return missing_state(
            "spectra", ["load_spectra"],
            "実験スペクトルが読み込まれていません。先に load_spectra を実行してください。")
    spectrum = session.spectra.get(spectrum_id)
    if spectrum is None:
        return json_payload({"error": {
            "code": "unknown_spectrum_id",
            "message": f"{spectrum_id} は読み込まれていません。",
            "available": sorted(session.spectra)[:20],
        }})
    return json_payload(round_floats({
        "spectrum_id": spectrum.spectrum_id,
        "precursor_mz": spectrum.precursor_mz,
        "polarity": spectrum.polarity,
        "retention_time": spectrum.retention_time,
        "peaks": [[p.mz, p.relative_intensity] for p in spectrum.peaks],
    }, MZ_DIGITS))
