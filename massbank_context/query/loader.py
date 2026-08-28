"""実験スペクトルの読み込み入口。拡張子で振り分ける。"""
from __future__ import annotations

from pathlib import Path

from massbank_context.core.models import QuerySpectrum
from massbank_context.query.dcl import read_dcl
from massbank_context.query.msp import read_msp

SUPPORTED_SUFFIXES = (".msp", ".mgf", ".dcl")


def load_spectra_file(path: str | Path) -> list[QuerySpectrum]:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in (".msp", ".mgf"):
        return read_msp(path)
    if suffix == ".dcl":
        return read_dcl(path)
    raise ValueError(
        f"対応していない拡張子です: {suffix}。"
        f"対応形式: {', '.join(SUPPORTED_SUFFIXES)}")
