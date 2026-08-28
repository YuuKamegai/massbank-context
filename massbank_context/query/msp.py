"""MSP / MGF 形式の実験スペクトルを読む。

MassBank 側の relative_intensity は基準ピーク 999 のスケールなので、
クエリ側も同じスケールに正規化してから照合する。
"""
from __future__ import annotations

from pathlib import Path

from massbank_context.core.models import Peak, QuerySpectrum

BASE_PEAK_SCALE = 999.0

_FIELD_ALIASES = {
    "PRECURSORMZ": "precursor_mz",
    "PEPMASS": "precursor_mz",
    "IONMODE": "polarity",
    "POLARITY": "polarity",
    "RETENTIONTIME": "retention_time",
    "RTINSECONDS": "retention_time",
}


def _normalise_polarity(raw: str | None) -> str | None:
    if not raw:
        return None
    text = raw.strip().upper()
    if text.startswith("P") or text == "1" or text == "+":
        return "POSITIVE"
    if text.startswith("N") or text == "0" or text == "-":
        return "NEGATIVE"
    return None


def _to_peaks(pairs: list[tuple[float, float]]) -> tuple[Peak, ...]:
    if not pairs:
        return ()
    top = max(intensity for _mz, intensity in pairs) or 1.0
    peaks = [Peak(mz=mz,
                  relative_intensity=intensity / top * BASE_PEAK_SCALE,
                  formula=None, neutral_loss=None)
             for mz, intensity in pairs]
    peaks.sort(key=lambda p: p.mz)
    return tuple(peaks)


def read_msp(path: str | Path) -> list[QuerySpectrum]:
    """MSP / MGF を読む。空行または BEGIN/END IONS でレコードが区切られる。"""
    path = Path(path)
    stem = path.stem
    spectra: list[QuerySpectrum] = []
    fields: dict[str, str] = {}
    pairs: list[tuple[float, float]] = []

    def flush() -> None:
        if not pairs and not fields:
            return
        spectra.append(QuerySpectrum(
            spectrum_id=f"{stem}:{len(spectra) + 1}",
            precursor_mz=float(fields["precursor_mz"].split()[0])
            if "precursor_mz" in fields else None,
            polarity=_normalise_polarity(fields.get("polarity")),
            peaks=_to_peaks(pairs),
            source=str(path),
            retention_time=float(fields["retention_time"])
            if "retention_time" in fields else None,
        ))
        fields.clear()
        pairs.clear()

    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line or line.upper() == "END IONS":
            flush()
            continue
        if line.upper() == "BEGIN IONS":
            continue
        if ":" in line and not line[0].isdigit():
            key, _sep, value = line.partition(":")
            alias = _FIELD_ALIASES.get(key.strip().upper())
            if alias:
                fields[alias] = value.strip()
            continue
        parts = line.replace("\t", " ").split()
        if len(parts) >= 2:
            try:
                pairs.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    flush()
    return spectra
