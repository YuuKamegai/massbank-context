"""セッション状態。可変状態の正準はこのモジュールの `session`。"""
from __future__ import annotations

from massbank_context.core.models import QuerySpectrum


class Session:
    """読み込んだクエリスペクトルを保持する。

    ピーク列はここに置き、ツールの戻り値には載せない（spec §12）。
    """

    def __init__(self) -> None:
        self.spectra: dict[str, QuerySpectrum] = {}
        self.source: str | None = None

    def reset(self) -> None:
        self.spectra = {}
        self.source = None

    def set_spectra(self, spectra: list[QuerySpectrum], source: str) -> None:
        self.spectra = {s.spectrum_id: s for s in spectra}
        self.source = source


session = Session()
