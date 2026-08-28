"""instrument_type / ms_type から索引の層を導出する。

spec §6.1 / §6.2。この 2 つの層は索引の分割単位であり、背景分布と濃縮度統計の
計算範囲を決める。測定原理の異なるスペクトルを同じ背景で数えることは無意味なので、
判定不能なものは "other" / "unknown" に落として明示的に扱う。

依存グラフの leaf（stdlib のみ）。
"""
from __future__ import annotations

#: 中性損失軸（precursor - fragment）を使える acquisition_class。
#: GC-EI-MS は precursor という概念自体を持たないため除外する。
NEUTRAL_LOSS_CAPABLE = frozenset({"LC-ESI-MS2", "MALDI-MS2"})

_ACQUISITION_TABLE = {
    ("LC-ESI", "MS2"): "LC-ESI-MS2",
    ("GC-EI", "MS"): "GC-EI-MS",
    ("MALDI", "MS2"): "MALDI-MS2",
}

_HIGH_ANALYZER_SUFFIXES = ("FT", "TOF")
_LOW_ANALYZERS = frozenset({"QQ", "QIT", "IT", "Q"})


def _tokens(instrument_type: str | None) -> list[str]:
    """'LC-ESI-ITFT' -> ['LC', 'ESI', 'ITFT']。

    部分文字列一致だと 'ESI' から 'EI' を誤って拾う経路が生まれるため、
    必ずトークンに割ってから判定する。
    """
    if not instrument_type:
        return []
    return [t for t in instrument_type.upper().replace("_", "-").split("-") if t]


def _ionization_family(tokens: list[str]) -> str | None:
    token_set = set(tokens)
    if "MALDI" in token_set:
        return "MALDI"
    if "ESI" in token_set:
        return "LC-ESI"
    if "EI" in token_set:
        return "GC-EI"
    return None


def derive_acquisition_class(instrument_type: str | None,
                             ms_type: str | None) -> str:
    """索引の分割単位を返す（spec §6.1）。

    判定できないものは "other" を返す。呼び出し側は "other" を索引対象外として
    扱ってよいが、黙って別の層に混ぜてはならない。
    """
    family = _ionization_family(_tokens(instrument_type))
    if family is None or not ms_type:
        return "other"
    return _ACQUISITION_TABLE.get((family, ms_type.strip().upper()), "other")


def derive_resolution_class(instrument_type: str | None,
                            resolution_value: float | None = None) -> str:
    """質量分解能の層を返す（spec §6.2）。

    mbo:resolution は実測で約半数が欠損していたため、判定は instrument_type の
    アナライザ部を主根拠とする。resolution_value は記録・補強のためだけに受け取り、
    判定を覆さない。
    """
    tokens = _tokens(instrument_type)
    if not tokens:
        return "unknown"
    analyzer = tokens[-1]
    if analyzer in _LOW_ANALYZERS:
        return "low"
    if any(analyzer.endswith(suffix) for suffix in _HIGH_ANALYZER_SUFFIXES):
        return "high"
    return "unknown"
