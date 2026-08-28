"""索引とクエリのデータモデル。依存グラフの leaf（classify のみ）。"""
from __future__ import annotations

from dataclasses import dataclass, field

from massbank_context.core.classify import (
    derive_acquisition_class,
    derive_resolution_class,
)


@dataclass(frozen=True)
class Peak:
    """1 本のピーク。

    neutral_loss は世代B TTL では mbo:neutral_loss から読み、世代A では
    precursor_mz - mz で後から算出する（spec §6.3）。算出できない場合は None。
    """
    mz: float
    relative_intensity: float
    formula: str | None = None
    neutral_loss: float | None = None


@dataclass(frozen=True)
class SpectrumRecord:
    """MassBank レコード 1 件（spec §6.7 の records + peaks）。"""
    accession: str
    inchikey: str | None
    precursor_mz: float | None
    precursor_type: str | None
    polarity: str | None
    instrument_type: str | None
    ms_type: str | None
    resolution_value: float | None
    collision_energy: str | None
    splash: str | None
    license: str | None
    molecular_formula: str | None
    exact_mass: float | None
    compound_names: tuple[str, ...] = ()
    ch_compound_class: str | None = None
    peaks: tuple[Peak, ...] = ()

    @property
    def acquisition_class(self) -> str:
        return derive_acquisition_class(self.instrument_type, self.ms_type)

    @property
    def resolution_class(self) -> str:
        return derive_resolution_class(self.instrument_type,
                                       self.resolution_value)


@dataclass(frozen=True)
class QuerySpectrum:
    """ユーザーが持ち込んだ実験スペクトル 1 本。"""
    spectrum_id: str
    precursor_mz: float | None
    polarity: str | None
    peaks: tuple[Peak, ...] = ()
    source: str = ""
    retention_time: float | None = None
