"""質量誤差の実測とビン幅の決定（spec §6.4）。

決め打ちのビン幅には根拠がない。MassBank は装置が混在しており、必要な
トレランスはレコード集団によって違う。同一 InChIKey の複数レコード
（装置・CE 違い）で同じフラグメントの m/z がどう散らばるかを測れば、
必要な幅が実測で出る。

測るのは「同じフラグメントを別レコードで測ったときの差」だけである。ここを
外すと、窓の中にたまたま入った**別**フラグメントまで誤差として数えてしまい、
出力は窓の幅そのものに追随する。実際に踏んだ症状がそれで、標本を集める窓を
25 倍動かすと採用ビン幅もその 85〜90% で追随していた（＝装置誤差を測って
いなかった）。対策は 2 段構え:

1. 対応付けを**相互最近傍**で決める。相手を持たないピークが、既に対に
   なっているピークへ寄生するのを防ぐ
2. それでも残る偶然の近接は、窓の外側（BACKGROUND_FROM 以遠）で密度を
   測り、**一様な背景として差し引く**

MassBank の高分解能スペクトルは 1 レコード平均 26.5 本のピークを持つので、
50 mDa の窓に無関係なフラグメントが入るのは日常茶飯事である。背景を引かずに
上位パーセンタイルを採ると、その混入がそのまま「装置誤差」に化ける。
"""
from __future__ import annotations

import bisect
from collections import defaultdict
from dataclasses import dataclass

from massbank_context.core.models import SpectrumRecord

MASS_ERROR_PERCENTILES = (50.0, 90.0, 95.0, 99.0)

#: 標本を集める窓。これ以上離れた相手は対応候補にしない。
MATCH_WINDOW = 0.05

#: この距離以遠の対応は、装置誤差ではなく別フラグメントの偶然の近接だと
#: みなし、背景密度の推定に使う。高分解能装置の実測誤差（Orbitrap で 1 mDa
#: 前後、QTOF で 10 mDa 前後）より十分外側に置く。
BACKGROUND_FROM = 0.02

#: 背景を差し引くためのヒストグラムの刻み。返す値はこの刻みに量子化され、
#: 該当ビンの上端になる（下端を返すと幅を過小に見積もり、本物の一致を
#: 取りこぼす側に倒れる）。
HISTOGRAM_STEP = 0.0005

DEFAULT_BIN_WIDTH = 0.01

#: これ未満の標本数しか出せなかった装置クラスは、ビン幅の決定に使わない。
#: 少数標本の p95 は装置誤差ではなく外れ値そのもので、それを採ると索引全体の
#: 照合窓が数十件のレコードに引きずられる（実際に踏んだ: 45 レコードの
#: MALDI-TOFTOF が 84,469 レコードの索引の窓を決めていた）。
#: 除外したクラスも mass_error_profile には残るので、「測ったが使わなかった」
#: ことは後から確認できる。
MIN_MASS_ERROR_SAMPLES = 200

#: ビン幅が覆うべき索引レコードの割合。クラス横断の最大値ではなく、
#: レコード数で重み付けした分位点を採る。
DEFAULT_RECORD_COVERAGE = 0.95


@dataclass(frozen=True)
class InstrumentErrorProfile:
    """1 つの装置クラスについて測った質量誤差。

    n_records / n_samples を値と一緒に運ぶ。これが無いと choose_bin_width は
    「その数字がどれだけの母集団から出たか」を知らずに採否を決めることになり、
    ごく少数のレコードしか無いクラスが索引全体のビン幅を決めてしまう。
    """

    instrument_class: str
    n_records: int
    n_samples: int
    percentiles: tuple[tuple[float, float], ...]

    def value_at(self, percentile: float) -> float | None:
        for pct, value in self.percentiles:
            if pct == percentile:
                return value
        return None


def _instrument_class(record: SpectrumRecord) -> str:
    return record.instrument_type or "unknown"


def _nearest(sorted_mz: list[float], value: float) -> tuple[float, float] | None:
    """sorted_mz の中で value に最も近い要素を (距離, 値) で返す。"""
    if not sorted_mz:
        return None
    position = bisect.bisect_left(sorted_mz, value)
    best: tuple[float, float] | None = None
    for index in (position - 1, position):
        if 0 <= index < len(sorted_mz):
            distance = abs(sorted_mz[index] - value)
            if best is None or distance < best[0]:
                best = (distance, sorted_mz[index])
    return best


def _mutual_deltas(left: list[float], right: list[float],
                   window: float) -> list[float]:
    """2 レコード間で相互最近傍になった対の |Δm/z| を返す。

    片方向の最近傍だと、相手を持たないピークが最寄りの他人を掴んでしまう。
    その距離は装置誤差ではなくフラグメント間の距離であり、窓の幅いっぱいまで
    伸びうる。両方向で最近傍が一致することを要求すれば、この寄生は落ちる。
    """
    deltas: list[float] = []
    for mz in left:
        forward = _nearest(right, mz)
        if forward is None or forward[0] > window:
            continue
        backward = _nearest(left, forward[1])
        if backward is None or backward[1] != mz:
            continue
        deltas.append(forward[0])
    return deltas


def _percentiles(deltas: list[float], window: float,
                 background_from: float) -> tuple[tuple[float, float], ...]:
    """背景を差し引いた分布の分位点を返す。

    背景（無関係なフラグメントの偶然の近接）は background_from 以遠の密度で
    測り、一様と仮定して信号領域から差し引く。差し引いた結果が空になった
    場合は「測定不能」として 0.0 を返す — choose_bin_width が 0 を採らずに
    fallback へ回すため、緩すぎる幅が黙って採用されることはない。
    """
    zero = tuple((pct, 0.0) for pct in MASS_ERROR_PERCENTILES)
    if not deltas or max(deltas) == 0.0:
        # 完全に同一のピーク列を持つ重複レコードだけの群。量子化して
        # 0 でない値を返すと、測っていない精度を主張することになる。
        return zero

    n_bins = max(1, int(round(window / HISTOGRAM_STEP)))
    histogram = [0] * n_bins
    for delta in deltas:
        histogram[min(n_bins - 1, int(delta / HISTOGRAM_STEP))] += 1

    first_background = min(n_bins, int(round(background_from / HISTOGRAM_STEP)))
    background_bins = histogram[first_background:]
    density = (sum(background_bins) / len(background_bins)
               if background_bins else 0.0)
    signal = [max(0.0, count - density)
              for count in histogram[:first_background]]
    total = sum(signal)
    if total <= 0:
        return zero

    values = []
    for percentile in MASS_ERROR_PERCENTILES:
        target = percentile / 100.0 * total
        cumulative = 0.0
        value = first_background * HISTOGRAM_STEP
        for index, count in enumerate(signal):
            cumulative += count
            if cumulative >= target:
                value = (index + 1) * HISTOGRAM_STEP
                break
        # 刻みの整数倍を浮動小数で作ると 0.0045000000000000005 のような値に
        # なる。この値はビン幅として index_meta に載り、ツールの戻り値にも
        # 出るので、ここで丸めておく（spec §12）。
        values.append((percentile, round(value, 6)))
    return tuple(values)


def measure_mass_error(
    records: list[SpectrumRecord],
    window: float = MATCH_WINDOW,
    background_from: float = BACKGROUND_FROM,
) -> dict[str, InstrumentErrorProfile]:
    """装置種別ごとの m/z 散らばりを測る。

    同一 (instrument_type, inchikey) のレコードを accession 順に並べ、隣り合う
    レコードの間で相互最近傍の対を作り、その |Δm/z| を標本とする。全ペアでは
    なく隣接ペアだけを見るのは、レコードを大量に持つ化合物（同一化合物の CE
    違いが数百件ある）で比較回数が二乗に膨らむのを避けるため。実測では全ペアと
    隣接ペアで p95 は 0.0090 と 0.0085 でほとんど変わらない。
    """
    grouped: dict[tuple[str, str], list[SpectrumRecord]] = defaultdict(list)
    for record in records:
        if not record.inchikey:
            continue
        grouped[(_instrument_class(record), record.inchikey)].append(record)

    deltas: dict[str, list[float]] = defaultdict(list)
    record_counts: dict[str, int] = defaultdict(int)
    for (instrument, _inchikey), group in grouped.items():
        if len(group) < 2:
            continue
        record_counts[instrument] += len(group)
        ordered = [sorted(peak.mz for peak in record.peaks)
                   for record in sorted(group, key=lambda r: r.accession)]
        for left, right in zip(ordered, ordered[1:]):
            deltas[instrument].extend(_mutual_deltas(left, right, window))

    profiles: dict[str, InstrumentErrorProfile] = {}
    for instrument, values in deltas.items():
        if not values:
            continue
        profiles[instrument] = InstrumentErrorProfile(
            instrument_class=instrument,
            n_records=record_counts[instrument],
            n_samples=len(values),
            percentiles=_percentiles(values, window, background_from),
        )
    return profiles


def choose_bin_width(profiles: dict[str, InstrumentErrorProfile],
                     percentile: float = 95.0,
                     coverage: float = DEFAULT_RECORD_COVERAGE,
                     min_samples: int = MIN_MASS_ERROR_SAMPLES,
                     fallback: float = DEFAULT_BIN_WIDTH) -> float:
    """索引全体で使う 1 つのビン幅を決める。

    索引は単一のビン幅で作るので、どの装置クラスに合わせるかを決める必要が
    ある。以前はクラス横断の最大値を採っていたが、それだと最も散らばる 1
    クラスが — レコード数が索引の 0.05% であっても — 全体の照合窓を決めて
    しまい、残り 99.95% の照合が不必要に緩くなる。

    そこで、レコード数で重み付けした分位点を採る: 値の小さい順にレコード数を
    積み、索引レコードの coverage 割合を覆えた時点の値を返す。標本数が
    min_samples に満たないクラスは統計として成立していないので除外する。

    測定できなかった場合は fallback を返し、index_meta に記録する。
    """
    weighted = [
        (value, profile.n_records)
        for profile in profiles.values()
        if profile.n_samples >= min_samples
        # value > 0 は下流の前提条件を守るためのもので、単なる衛生ではない。
        # core/binning.py の to_bin は bin_width <= 0 で ValueError を投げるため、
        # ここで 0.0 を返すと索引ビルド全体が落ちる。散らばりがちょうど 0 に
        # なるのは実在する（同一スペクトルの重複登録）。
        for value in (profile.value_at(percentile),)
        if value is not None and value > 0
    ]
    if not weighted:
        return fallback
    total = sum(n_records for _value, n_records in weighted)
    target = coverage * total
    cumulative = 0
    for value, n_records in sorted(weighted):
        cumulative += n_records
        if cumulative >= target:
            return value
    return max(value for value, _n in weighted)
