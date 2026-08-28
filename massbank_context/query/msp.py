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
    "RETENTIONTIME": "retention_time",          # 分
    "RTINSECONDS": "retention_time_seconds",    # 秒。分へ変換する
}


def _split_header(line: str) -> tuple[str, str] | None:
    """ヘッダ行を (キー, 値) に割る。MSP は `KEY: value`、MGF は `KEY=value`。

    最初に現れた区切り文字を採るので、`TITLE=a:b` も `NAME: a=b` も正しく割れる。
    値そのものに区切り文字が含まれていても、キー側は壊れない。
    """
    candidates = [line.find(sep) for sep in (":", "=")]
    positions = [pos for pos in candidates if pos > 0]
    if not positions:
        return None
    cut = min(positions)
    return line[:cut].strip(), line[cut + 1:].strip()


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
        # RT の単位は分に統一する。MSP の RETENTIONTIME は分、
        # MGF の RTINSECONDS は秒。混ざったまま同じフィールドに入れると
        # 後から単位のずれに気付けない。
        if "retention_time" in fields:
            retention_time = float(fields["retention_time"])
        elif "retention_time_seconds" in fields:
            retention_time = float(fields["retention_time_seconds"]) / 60.0
        else:
            retention_time = None
        spectra.append(QuerySpectrum(
            spectrum_id=f"{stem}:{len(spectra) + 1}",
            precursor_mz=float(fields["precursor_mz"].split()[0])
            if "precursor_mz" in fields else None,
            polarity=_normalise_polarity(fields.get("polarity")),
            peaks=_to_peaks(pairs),
            source=str(path),
            retention_time=retention_time,
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
        header = _split_header(line) if not line[0].isdigit() else None
        if header:
            key, value = header
            alias = _FIELD_ALIASES.get(key.strip().upper())
            if alias:
                fields[alias] = value.strip()
            continue
        parts = line.replace("\t", " ").split()
        if len(parts) >= 2:
            try:
                pairs.append((float(parts[0]), float(parts[1])))
            except ValueError:
                # ヘッダでもピークでもない紛れ込んだ行（コメント等）は
                # 意図的に読み飛ばす。
                continue
    flush()
    return spectra
