"""世代A・世代B の両リーダが共有するリテラル取得ヘルパ（spec §3.3〜§3.4）。

`_one` / `_text` は単なる走査の便宜だが、`_number` は上流データ品質に関する
方針そのもの（欠損は無言・不正は warning 1 回・例外は投げない）を埋め込んでいる。
索引ビルドが 54,183 ファイルを走査するため、この方針が世代ごとにばらばらだと
片方だけ静かに退行しても気付けない。ここに正本を 1 つ置き、両リーダから import する。
"""
from __future__ import annotations

import logging

import rdflib

logger = logging.getLogger(__name__)


def _one(graph: rdflib.Graph, subject, predicate):
    for value in graph.objects(subject, predicate):
        return value
    return None


def _text(value) -> str | None:
    return None if value is None else str(value)


def _number(value) -> float | None:
    """リテラルを float にする。欠損は None、変換不能も None を返す。

    欠損（value が None）は正常。mbo:resolution は実測で約半数のレコードに無い。
    一方「値はあるが float にできない」は上流データの異常であり、54,183 ファイルの
    走査で黙って落ちると気付けないため、そちらだけ warning を出す。
    例外は投げない — 不正リテラル 1 個でレコード全体を失うほうが損失が大きい。
    """
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        logger.warning("float に変換できないリテラルを無視しました: %r", value)
        return None
