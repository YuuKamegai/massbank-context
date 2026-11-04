"""LLM へ返す文字列の直列化。依存グラフの leaf（stdlib のみ）。

戻り値は LLM の文脈をそのまま占める。json.dumps(..., indent=2) は 1 ネストにつき
改行と空白を全キーに掛けるため、Lipidmix の実測では戻り値の 15〜57% が空白だった。
人間向けの整形はクライアント側の責務なので、サーバは最小形で返す。
"""
from __future__ import annotations

import json

DEFAULT_FLOAT_DIGITS = 6
MZ_DIGITS = 4


def json_payload(obj) -> str:
    """ツールの戻り値を最小トークンの JSON 文字列にする。"""
    return json.dumps(obj, ensure_ascii=False, separators=(",", ":"))


SMALL_FLOAT_DIGITS = 8


def round_floats(obj, digits: int = DEFAULT_FLOAT_DIGITS, *, keep_small: bool = False):
    """入れ子の dict/list を辿って float を丸める（非破壊・型は保つ）。

    bool は int の派生だが数値として丸めてはいけないので素通しする。
    keep_small=True なら 10^-digits 未満の非零値（小さな p 値など）を 0 に潰さず、
    有効数字 SMALL_FLOAT_DIGITS 桁で残す。既存ツールの書式を変えないため既定は False。
    """
    if isinstance(obj, bool):
        return obj
    if isinstance(obj, float):
        if keep_small and 0 < abs(obj) < 10 ** -digits:
            return float(format(obj, f".{SMALL_FLOAT_DIGITS}g"))
        return round(obj, digits)
    if isinstance(obj, list):
        return [round_floats(item, digits, keep_small=keep_small) for item in obj]
    if isinstance(obj, tuple):
        return tuple(round_floats(item, digits, keep_small=keep_small) for item in obj)
    if isinstance(obj, dict):
        return {key: round_floats(value, digits, keep_small=keep_small)
                for key, value in obj.items()}
    return obj


_TSV_ESCAPES = str.maketrans({"\\": "\\\\", "\t": "\\t", "\n": "\\n", "\r": "\\r"})


def _tsv_cell(value) -> str:
    if value is None:
        return ""
    # 値中のタブ・改行で列や行がずれないよう、バックスラッシュ表記に置き換える。
    return str(value).translate(_TSV_ESCAPES)


def as_tsv(columns: list[str], rows: list[list]) -> str:
    """行が並ぶ一覧を TSV で返す（列名は 1 回だけ）。

    JSON の配列オブジェクトは行ごとにキー名を繰り返すため、行数に比例して無駄が増える。
    セル中のバックスラッシュ・タブ・改行は \\\\ \\t \\n \\r の表記にエスケープする。
    """
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(_tsv_cell(v) for v in row))
    return "\n".join(lines)


def rows_tsv(rows: list[dict]) -> str:
    """dict の行を TSV にする。列は初出順の和集合、入れ子の値は JSON 文字列にする。"""
    columns = list(dict.fromkeys(key for row in rows for key in row))
    return as_tsv(columns, [[_scalar(row.get(key)) for key in columns] for row in rows])


def _scalar(value):
    if isinstance(value, (dict, list, tuple)):
        return json_payload(value)
    return value
