"""MS-DIAL .dcl（デコンボリューション済み MS/MS）を読む。

Lipidmix lipidmix/dcl/reader.py の deserialize_dcl() を移植したもの。
バイナリレイアウトは MS-DIAL の MSDecResult に対応する。推測で変えないこと。
"""
from __future__ import annotations

import logging
import struct
from pathlib import Path

from massbank_context.core.models import Peak, QuerySpectrum
from massbank_context.query.msp import BASE_PEAK_SCALE

logger = logging.getLogger(__name__)

DCL_MAGIC = b"DC"


def _polarity_from_ion_mode(ion_mode: int) -> str | None:
    """MS-DIAL の IonMode 列挙を極性に写す。

    0=Positive / 1=Negative / 2=Both（Lipidmix docs/output_format/dcl.md）。
    Both と未知の値は None にする。極性は照合の絞り込み条件であり、
    None は「極性で絞らない」を意味する。決め打ちで片方に倒すと候補が不当に消える。
    """
    if ion_mode == 0:
        return "POSITIVE"
    if ion_mode == 1:
        return "NEGATIVE"
    if ion_mode != 2:
        logger.warning("未知の ion_mode を極性不明として扱います: %r", ion_mode)
    return None


def read_dcl(path: str | Path) -> list[QuerySpectrum]:
    path = Path(path)
    raw = path.read_bytes()
    if raw[:2] != DCL_MAGIC:
        raise ValueError(f"{path} は .dcl 形式ではありません（先頭が 'DC' でない）。")
    try:
        return _parse_records(raw, path)
    except (struct.error, IndexError) as exc:
        # struct.error は Exception 直下で ValueError の派生ではない。そのまま
        # 投げるとツール層の except (OSError, ValueError) をすり抜け、機械可読な
        # 封筒（spec §8）ではなくトレースバックになる。書き込み中断やコピー中断で
        # 実際に生じる状態なので、ここで契約に載る型へ変換する。
        raise ValueError(
            f"{path} の .dcl バイナリを読み切れませんでした"
            f"（途中で切れているか、レイアウトが想定と異なります）: {exc}"
        ) from exc


def _parse_records(raw: bytes, path: Path) -> list[QuerySpectrum]:
    count = struct.unpack_from("<i", raw, 7)[0]
    seekpoints = struct.unpack_from(f"<{count}q", raw, 11)

    spectra: list[QuerySpectrum] = []
    for index, seekpoint in enumerate(seekpoints):
        offset = seekpoint
        precursor_mz = struct.unpack_from("<d", raw, offset + 16)[0]
        ion_mode = struct.unpack_from("<i", raw, offset + 24)[0]
        retention_time = struct.unpack_from("<d", raw, offset + 28)[0]
        offset += 60 + 40 + 20                     # Scan / Quant / Scoring
        spectra_n = struct.unpack_from("<3i", raw, offset)[0]
        offset += 12

        pairs: list[tuple[float, float]] = []
        for _ in range(spectra_n):
            mz, intensity = struct.unpack_from("<2d", raw, offset)
            offset += 20                            # double + double + int32
            pairs.append((mz, intensity))

        top = max((intensity for _mz, intensity in pairs), default=0.0) or 1.0
        peaks = tuple(sorted(
            (Peak(mz=mz, relative_intensity=intensity / top * BASE_PEAK_SCALE,
                  formula=None, neutral_loss=None)
             for mz, intensity in pairs),
            key=lambda p: p.mz))

        spectra.append(QuerySpectrum(
            spectrum_id=f"{path.stem}:{index}",
            precursor_mz=precursor_mz,
            polarity=_polarity_from_ion_mode(ion_mode),
            peaks=peaks,
            source=str(path),
            retention_time=retention_time,
        ))
    return spectra
