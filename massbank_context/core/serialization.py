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


def round_floats(obj, digits: int = DEFAULT_FLOAT_DIGITS):
    """入れ子の dict/list を辿って float を丸める（非破壊・型は保つ）。

    bool は int の派生だが数値として丸めてはいけないので素通しする。
    """
    if isinstance(obj, bool):
        return obj
    if isinstance(obj, float):
        return round(obj, digits)
    if isinstance(obj, list):
        return [round_floats(item, digits) for item in obj]
    if isinstance(obj, tuple):
        return tuple(round_floats(item, digits) for item in obj)
    if isinstance(obj, dict):
        return {key: round_floats(value, digits) for key, value in obj.items()}
    return obj


def as_tsv(columns: list[str], rows: list[list]) -> str:
    """行が並ぶ一覧を TSV で返す（列名は 1 回だけ）。

    JSON の配列オブジェクトは行ごとにキー名を繰り返すため、行数に比例して無駄が増える。
    """
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join("" if v is None else str(v) for v in row))
    return "\n".join(lines)
