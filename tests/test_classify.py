import pytest

from massbank_context.core.classify import (
    NEUTRAL_LOSS_CAPABLE,
    derive_acquisition_class,
    derive_resolution_class,
)


@pytest.mark.parametrize("instrument_type,ms_type,expected", [
    # spec §3.5 の実測値
    ("EI-B", "MS", "GC-EI-MS"),
    ("LC-ESI-ITFT", "MS2", "LC-ESI-MS2"),
    ("LC-ESI-QQ", "MS2", "LC-ESI-MS2"),
    ("LC-ESI-QTOF", "MS2", "LC-ESI-MS2"),
    ("LC-ESI-QIT", "MS2", "LC-ESI-MS2"),
    ("ESI-QTOF", "MS2", "LC-ESI-MS2"),      # LC 接頭辞が無い版
    ("LC-ESI-QFT", "MS2", "LC-ESI-MS2"),
    ("MALDI-TOFTOF", "MS2", "MALDI-MS2"),
    # 判定不能
    ("LC-ESI-QTOF", "MS3", "other"),
    ("EI-B", "MS2", "other"),               # EI の MS2 は想定外
    (None, "MS2", "other"),
    ("LC-ESI-QTOF", None, "other"),
])
def test_derive_acquisition_class(instrument_type, ms_type, expected):
    assert derive_acquisition_class(instrument_type, ms_type) == expected


def test_esi_is_not_mistaken_for_ei():
    """'ESI' の部分文字列一致で 'EI' を拾ってはならない。"""
    assert derive_acquisition_class("LC-ESI-ITFT", "MS2") == "LC-ESI-MS2"


@pytest.mark.parametrize("instrument_type,expected", [
    ("LC-ESI-ITFT", "high"),     # Orbitrap 系
    ("LC-ESI-QFT", "high"),
    ("LC-ESI-QTOF", "high"),
    ("MALDI-TOFTOF", "high"),
    ("ESI-QTOF", "high"),
    ("LC-ESI-QQ", "low"),        # 三連四重極
    ("LC-ESI-QIT", "low"),       # イオントラップ
    ("EI-B", "unknown"),         # 磁場型。MassBank の EI は公称質量が多い
    (None, "unknown"),
])
def test_derive_resolution_class(instrument_type, expected):
    assert derive_resolution_class(instrument_type) == expected


def test_resolution_value_does_not_override_class():
    """mbo:resolution は補強にのみ使い、instrument_type の判定を覆さない。

    spec §6.2: resolution は約半数で欠損するため主根拠にできない。
    """
    assert derive_resolution_class("LC-ESI-QTOF", 7000.0) == "high"
    assert derive_resolution_class("LC-ESI-QQ", 50000.0) == "low"


def test_neutral_loss_capable_excludes_gc_ei():
    """GC-EI-MS は precursor を持たないため中性損失軸を使えない（spec §6.1）。"""
    assert "LC-ESI-MS2" in NEUTRAL_LOSS_CAPABLE
    assert "GC-EI-MS" not in NEUTRAL_LOSS_CAPABLE
