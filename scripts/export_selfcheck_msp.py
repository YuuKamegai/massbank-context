"""索引の元になった MassBank レコードを .msp に書き出す（既知解テスト用）。

照合器の健全性は、正解の分かっている問いでしか測れない。索引に入っている
レコードそのものをクエリとして投げれば、**自分の accession が 1 位に来なければ
おかしい**という強い判定ができる。実データの評価に進む前にこれを通す。

出力は 2 つ:

- `<出力先>.msp`         — クエリ。spectrum_id は読み込み側で `<ファイル名>:<連番>` になる
- `<出力先>.answer.tsv`  — 正解表（spectrum_id / accession / inchikey / instrument_type）

使い方:

    python scripts/export_selfcheck_msp.py --ttl-dir <TTLのディレクトリ> \
        --output C:/Users/yuu18/datasets/_selfcheck/massbank_selfcheck

装置クラスごとに --per-class 件ずつ拾う。CE 違いの同一化合物が索引に多数
あるため、1 位が「自分」ではなく「同一 InChIKey の別レコード」になることは
起こりうる。answer.tsv に inchikey も出しているのはそれを区別するため。
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from massbank_context.index.writer import is_indexable
from massbank_context.ttl.reader import read_ttl_file_with_generation

#: 1 装置クラスあたりに拾う既定件数
DEFAULT_PER_CLASS = 5

#: ピークが少なすぎるレコードは既知解テストの題材にならない
#: （偶然の一致でも 1 位になれてしまう）
MIN_PEAKS = 5


def collect(ttl_dir: Path, per_class: int, resolution_filter: str,
            acquisition_class: str, stride: int) -> list:
    """装置クラスごとに per_class 件だけ拾う。全走査はしない。"""
    paths = sorted(ttl_dir.rglob("*.ttl"))[::stride]
    picked: dict[str, list] = defaultdict(list)
    for path in paths:
        if all(len(v) >= per_class for v in picked.values()) and picked:
            # 既に集まったクラスばかりでも、未出現のクラスがある可能性は
            # 残るので走査は続ける。打ち切りは呼び出し側の stride で行う。
            pass
        try:
            _generation, records = read_ttl_file_with_generation(path)
        except Exception:
            continue
        for record in records:
            if not is_indexable(record, resolution_filter):
                continue
            if record.acquisition_class != acquisition_class:
                continue
            if len(record.peaks) < MIN_PEAKS:
                continue
            instrument = record.instrument_type or "unknown"
            if len(picked[instrument]) >= per_class:
                continue
            picked[instrument].append(record)
    return [record for records in picked.values() for record in records]


def write_msp(records: list, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    msp = output.with_suffix(".msp")
    answer = output.with_suffix(".answer.tsv")
    stem = msp.stem

    blocks = []
    rows = ["spectrum_id\taccession\tinchikey\tinstrument_type\tn_peaks"]
    for number, record in enumerate(records, start=1):
        lines = [f"NAME: {record.accession}"]
        if record.precursor_mz is not None:
            lines.append(f"PRECURSORMZ: {record.precursor_mz}")
        if record.polarity:
            lines.append(f"IONMODE: {record.polarity}")
        lines.append(f"Num Peaks: {len(record.peaks)}")
        lines.extend(f"{peak.mz}\t{peak.relative_intensity}"
                     for peak in record.peaks)
        blocks.append("\n".join(lines))
        rows.append("\t".join([
            f"{stem}:{number}", record.accession, record.inchikey,
            record.instrument_type or "", str(len(record.peaks)),
        ]))

    msp.write_text("\n\n".join(blocks) + "\n", encoding="utf-8")
    answer.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"{msp} に {len(records)} 件")
    print(f"{answer} に正解表")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ttl-dir", required=True)
    parser.add_argument("--output", required=True,
                        help="拡張子を除いた出力先。.msp と .answer.tsv を作る")
    parser.add_argument("--per-class", type=int, default=DEFAULT_PER_CLASS)
    parser.add_argument("--resolution-filter", default="high")
    parser.add_argument("--acquisition-class", default="LC-ESI-MS2")
    parser.add_argument("--stride", type=int, default=200,
                        help="TTL を何ファイルおきに見るか（全走査は分かかる）")
    args = parser.parse_args()

    records = collect(Path(args.ttl_dir), args.per_class,
                      args.resolution_filter, args.acquisition_class,
                      args.stride)
    if not records:
        raise SystemExit("条件に合うレコードが見つかりませんでした。"
                         "--stride を小さくするか層の条件を見直してください。")
    write_msp(records, Path(args.output))


if __name__ == "__main__":
    main()
