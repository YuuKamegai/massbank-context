"""質量誤差の実測とビン幅の決定（spec §6.4）。

決め打ちのビン幅には根拠がない。MassBank は装置が混在しており、必要なトレランスは
レコード集団によって違う。同一 InChIKey の複数レコード（装置・CE 違い）で同じ
フラグメントの m/z がどう散らばるかを測れば、必要な幅が実測で出る。
"""
from __future__ import annotations

from collections import defaultdict

from massbank_context.core.models import SpectrumRecord

MASS_ERROR_PERCENTILES = (50.0, 90.0, 95.0, 99.0)

#: これ以上離れたピークは別フラグメントとみなす。装置間差より十分大きく、
#: 隣接フラグメント（多くは 1 Da 以上離れる）より十分小さい値。
DEFAULT_CLUSTER_TOLERANCE = 0.05

DEFAULT_BIN_WIDTH = 0.01


def _percentile(values: list[float], percentile: float) -> float:
    """線形補間なしの単純な順位統計。値が少ないときも壊れない。"""
    if not values:
        return 0.0
    ordered = sorted(values)
    # round() は x.5 を偶数側に丸める（銀行家丸め）。n=2 のとき 50 パーセンタイルは
    # round(0.5)=0 で常に小さいほうを選ぶ。補間しない設計なので許容する。
    index = min(len(ordered) - 1,
                max(0, int(round(percentile / 100.0 * (len(ordered) - 1)))))
    return ordered[index]


def _instrument_class(record: SpectrumRecord) -> str:
    return record.instrument_type or "unknown"


def measure_mass_error(
    records: list[SpectrumRecord],
    cluster_tolerance: float = DEFAULT_CLUSTER_TOLERANCE,
) -> dict[str, list[tuple[float, float]]]:
    """装置種別ごとの m/z 散らばりを返す。

    同一 (instrument_type, inchikey) のレコード群を集め、各ピークを
    cluster_tolerance 内で貪欲にクラスタリングして、**2 つ以上の異なる
    レコード由来の m/z を含む**クラスタの最大幅（max - min）を散らばりの
    標本とする。

    「2 点以上」ではなく「2 レコード以上」であることが要点。測っているのは
    同じフラグメントを別々に測ったときのばらつきであり、1 本のスペクトルの
    中で近接している 2 本は単に別のフラグメントである。cluster_tolerance の
    既定 0.05 は測ろうとしている誤差の 10 倍ほどあるため、その距離を装置誤差
    として数えると、choose_bin_width が採る最大値経由で索引全体のビンが太る。
    """
    # m/z にレコードの識別子を付けて運ぶ。records 内の位置を使うのは、
    # accession が欠けている・重複しているレコードでも必ず一意になるため。
    grouped: dict[tuple[str, str], list[tuple[float, int]]] = defaultdict(list)
    for record_index, record in enumerate(records):
        if not record.inchikey:
            continue
        grouped[(_instrument_class(record), record.inchikey)].extend(
            (peak.mz, record_index) for peak in record.peaks)

    def _flush(instrument: str, cluster: list[tuple[float, int]],
               spreads: dict[str, list[float]]) -> None:
        if len({record_index for _mz, record_index in cluster}) < 2:
            return
        spreads[instrument].append(cluster[-1][0] - cluster[0][0])

    spreads: dict[str, list[float]] = defaultdict(list)
    for (instrument, _inchikey), entries in grouped.items():
        if len({record_index for _mz, record_index in entries}) < 2:
            continue
        cluster: list[tuple[float, int]] = []
        for mz, record_index in sorted(entries):
            if cluster and mz - cluster[0][0] > cluster_tolerance:
                _flush(instrument, cluster, spreads)
                cluster = []
            cluster.append((mz, record_index))
        if cluster:
            _flush(instrument, cluster, spreads)

    profile: dict[str, list[tuple[float, float]]] = {}
    for instrument, values in spreads.items():
        if not values:
            continue
        profile[instrument] = [
            (percentile, _percentile(values, percentile))
            for percentile in MASS_ERROR_PERCENTILES
        ]
    return profile


def choose_bin_width(profile: dict[str, list[tuple[float, float]]],
                     percentile: float = 95.0,
                     fallback: float = DEFAULT_BIN_WIDTH) -> float:
    """索引全体で使う 1 つのビン幅を決める。

    索引は単一のビン幅で作るため、最も散らばる装置に合わせる。
    測定できなかった場合は fallback を返し、index_meta に記録する。
    """
    candidates = [
        value
        for entries in profile.values()
        for pct, value in entries
        if pct == percentile and value > 0
    ]
    # value > 0 は下流の前提条件を守るためのもので、単なる衛生ではない。
    # core/binning.py の to_bin は bin_width <= 0 で ValueError を投げるため、
    # ここで 0.0 を返すと索引ビルド全体が落ちる。
    # 散らばりがちょうど 0 になるのは実在する: MassBank には同一化合物の
    # 同一ピーク列を持つレコード（同じスペクトルの重複登録）が頻出する。
    # その場合は「測定不能」と同じく fallback を返す — 幅 0 のビンは使えないため。
    if not candidates:
        return fallback
    return max(candidates)
