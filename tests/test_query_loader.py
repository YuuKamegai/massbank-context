import logging
import struct
from pathlib import Path

import pytest

from massbank_context.query.dcl import _polarity_from_ion_mode, read_dcl
from massbank_context.query.loader import load_spectra_file
from massbank_context.query.msp import read_msp

FIXTURE = Path(__file__).parent / "fixtures" / "query" / "sample.msp"
MGF_FIXTURE = Path(__file__).parent / "fixtures" / "query" / "sample.mgf"


def test_reads_two_spectra():
    assert len(read_msp(FIXTURE)) == 2


def test_first_spectrum_fields():
    spectrum = read_msp(FIXTURE)[0]
    assert spectrum.precursor_mz == 230.0552
    assert spectrum.polarity == "POSITIVE"
    assert spectrum.retention_time == 0.876
    assert len(spectrum.peaks) == 3


def test_polarity_is_normalised_to_upper():
    """MSP の 'Positive' / 'Negative' を索引側の表記に合わせる。"""
    spectra = read_msp(FIXTURE)
    assert spectra[0].polarity == "POSITIVE"
    assert spectra[1].polarity == "NEGATIVE"


def test_relative_intensity_is_normalised_to_999():
    """装置間で絶対強度の桁が違うため、基準ピーク 999 に正規化する。"""
    spectrum = read_msp(FIXTURE)[0]
    assert max(p.relative_intensity for p in spectrum.peaks) == 999.0


def test_peaks_sorted_by_mz():
    spectrum = read_msp(FIXTURE)[0]
    mzs = [p.mz for p in spectrum.peaks]
    assert mzs == sorted(mzs)


def test_spectrum_ids_are_unique_and_stable():
    spectra = read_msp(FIXTURE)
    ids = [s.spectrum_id for s in spectra]
    assert ids == ["sample:1", "sample:2"]


def test_loader_dispatches_by_extension():
    assert len(load_spectra_file(FIXTURE)) == 2


def test_loader_rejects_unknown_extension(tmp_path):
    path = tmp_path / "x.csv"
    path.write_text("a,b\n", encoding="utf-8")
    with pytest.raises(ValueError, match="対応していない"):
        load_spectra_file(path)


@pytest.mark.parametrize(
    ("ion_mode", "expected"),
    [
        (0, "POSITIVE"),
        (1, "NEGATIVE"),
        (2, None),
        (99, None),
    ],
)
def test_polarity_from_ion_mode(ion_mode, expected):
    """MS-DIAL の IonMode 列挙は 0=Positive/1=Negative/2=Both の3値。

    Both と未知の値は極性で絞り込まないよう None にする（docs/output_format/dcl.md）。
    """
    assert _polarity_from_ion_mode(ion_mode) == expected


def test_polarity_from_ion_mode_warns_only_on_unknown_value(caplog):
    """2 (Both) は仕様どおりの値なので warning を出さない。未知の値だけ出す。"""
    with caplog.at_level(logging.WARNING, logger="massbank_context.query.dcl"):
        _polarity_from_ion_mode(2)
    assert not caplog.records

    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="massbank_context.query.dcl"):
        _polarity_from_ion_mode(99)
    assert len(caplog.records) == 1
    assert "99" in caplog.records[0].getMessage()


def test_reads_mgf_equals_delimited_headers():
    """MGF は KEY=value 構文。コロンしか見ないとヘッダが全部落ちる。"""
    spectra = read_msp(MGF_FIXTURE)
    assert len(spectra) == 1
    assert spectra[0].precursor_mz == 267.0855      # PEPMASS の第 1 トークン
    assert len(spectra[0].peaks) == 2
    assert spectra[0].polarity is None               # IONMODE が無い


def test_mgf_rtinseconds_is_converted_to_minutes():
    """本プロジェクトの RT は分に統一する。MGF の RTINSECONDS=52.5 秒 -> 0.875 分。"""
    spectrum = read_msp(MGF_FIXTURE)[0]
    assert spectrum.retention_time == pytest.approx(0.875)


def test_msp_retentiontime_stays_in_minutes():
    """MSP の RETENTIONTIME はもともと分なので変換しない。"""
    spectrum = read_msp(FIXTURE)[0]
    assert spectrum.retention_time == pytest.approx(0.876)


def test_loader_dispatches_mgf_by_extension():
    assert len(load_spectra_file(MGF_FIXTURE)) == 1


# --- .dcl: 合成バイナリでの検証 -------------------------------------------
#
# 実物の .dcl ファイルはこの環境に無いため、read_dcl() のバイナリ解析経路
# （seekpoint テーブル・scan/quant/scoring/counts ブロックのオフセット・
# 20 バイト刻みのピーク読み取り）を確かめるすべが無かった。ここではテストの
# 中でバイト列を struct.pack で組み立てて tmp_path に書き、実際に読み返す。
# レイアウトは massbank_context/query/dcl.py の read_dcl() docstring
# （= このタスクブリーフのバイナリレイアウト）と同じ。

_SCAN_BLOCK_SIZE = 60
_QUANT_BLOCK_SIZE = 40
_SCORING_BLOCK_SIZE = 20
_COUNTS_BLOCK_SIZE = 12
_PEAK_STRIDE = 20


def _build_dcl_bytes(records: list[dict]) -> bytes:
    """テスト用に .dcl バイナリを合成する。"""
    header_len = 11 + 8 * len(records)
    record_sizes = [
        _SCAN_BLOCK_SIZE + _QUANT_BLOCK_SIZE + _SCORING_BLOCK_SIZE
        + _COUNTS_BLOCK_SIZE + _PEAK_STRIDE * len(record["peaks"])
        for record in records
    ]
    seekpoints = []
    offset = header_len
    for size in record_sizes:
        seekpoints.append(offset)
        offset += size
    total_len = offset

    buf = bytearray(total_len)
    buf[0:2] = b"DC"
    struct.pack_into("<i", buf, 2, 1)              # version
    buf[6] = 0                                      # has_annotation
    struct.pack_into("<i", buf, 7, len(records))
    struct.pack_into(f"<{len(records)}q", buf, 11, *seekpoints)

    for seekpoint, record in zip(seekpoints, records):
        o = seekpoint
        # scan block（先頭 8 バイトは未使用のまま 0 埋め）
        struct.pack_into("<i", buf, o + 8, record["scan_id"])
        struct.pack_into("<i", buf, o + 12, record["raw_spec_id"])
        struct.pack_into("<d", buf, o + 16, record["precursor_mz"])
        struct.pack_into("<i", buf, o + 24, record["ion_mode"])
        struct.pack_into("<d", buf, o + 28, record["retention_time"])
        o += _SCAN_BLOCK_SIZE
        o += _QUANT_BLOCK_SIZE                      # quant block: 0 埋めのまま
        o += _SCORING_BLOCK_SIZE                    # scoring block: 0 埋めのまま
        struct.pack_into("<i", buf, o, len(record["peaks"]))  # counts[0] = ピーク本数
        o += _COUNTS_BLOCK_SIZE
        for mz, intensity in record["peaks"]:
            struct.pack_into("<d", buf, o, mz)
            struct.pack_into("<d", buf, o + 8, intensity)
            o += _PEAK_STRIDE                        # 末尾 int32（未使用）ぶんスキップ
    return bytes(buf)


_DCL_RECORDS = [
    {
        "scan_id": 1, "raw_spec_id": 100,
        "precursor_mz": 230.0552, "ion_mode": 0, "retention_time": 0.876,
        "peaks": [(86.0349, 500000.0), (171.0069, 1000000.0)],
    },
    {
        "scan_id": 2, "raw_spec_id": 101,
        "precursor_mz": 267.0855, "ion_mode": 1, "retention_time": 1.234,
        "peaks": [(110.0349, 250000.0), (135.0434, 999000.0)],
    },
]


@pytest.fixture
def dcl_path(tmp_path):
    path = tmp_path / "synthetic.dcl"
    path.write_bytes(_build_dcl_bytes(_DCL_RECORDS))
    return path


def test_dcl_reads_all_records(dcl_path):
    assert len(read_dcl(dcl_path)) == 2


def test_dcl_precursor_mz_and_retention_time(dcl_path):
    spectra = read_dcl(dcl_path)
    assert spectra[0].precursor_mz == pytest.approx(230.0552)
    assert spectra[0].retention_time == pytest.approx(0.876)
    assert spectra[1].precursor_mz == pytest.approx(267.0855)
    assert spectra[1].retention_time == pytest.approx(1.234)


def test_dcl_polarity_from_ion_mode(dcl_path):
    spectra = read_dcl(dcl_path)
    assert spectra[0].polarity == "POSITIVE"   # ion_mode 0
    assert spectra[1].polarity == "NEGATIVE"   # ion_mode 1


def test_dcl_peak_count_and_order(dcl_path):
    for spectrum, record in zip(read_dcl(dcl_path), _DCL_RECORDS):
        assert len(spectrum.peaks) == len(record["peaks"])
        mzs = [p.mz for p in spectrum.peaks]
        assert mzs == sorted(mzs)


def test_dcl_base_peak_is_normalised_to_999(dcl_path):
    for spectrum in read_dcl(dcl_path):
        assert max(p.relative_intensity for p in spectrum.peaks) == 999.0


def test_truncated_dcl_raises_value_error_naming_the_file(tmp_path):
    """途中で切れた .dcl は struct.error ではなく ValueError にする（spec §8）。

    struct.error は Exception 直下で ValueError の派生ではないため、
    ツール層の except (OSError, ValueError) をすり抜けてトレースバックになる。
    ビルド中断・コピー中断で実際に生じる状態なので、機械可読な封筒に載せられる
    型へ変換する。ファイル名を必ず含める（どのファイルが壊れているか分かるため）。
    """
    path = tmp_path / "truncated.dcl"
    path.write_bytes(_build_dcl_bytes(_DCL_RECORDS)[:30])
    with pytest.raises(ValueError, match="truncated.dcl"):
        read_dcl(path)
