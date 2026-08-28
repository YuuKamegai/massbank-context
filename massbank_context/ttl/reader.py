"""TTL ファイルを読む入口。世代を判定して適切なリーダへ振り分ける。"""
from __future__ import annotations

from pathlib import Path

import rdflib

from massbank_context.core.models import SpectrumRecord
from massbank_context.ttl.generation import GEN_A, GEN_B, detect_generation
from massbank_context.ttl.reader_gen_a import parse_gen_a
from massbank_context.ttl.reader_gen_b import parse_gen_b


def read_ttl_file(path: str | Path) -> list[SpectrumRecord]:
    """1 ファイルを読む。世代が判定できない場合は空リストを返す。

    索引ビルドは数万ファイルを走査するため、1 件の異常でビルド全体を止めない。
    呼び出し側（CLI）が件数を集計して報告する。
    """
    graph = rdflib.Graph()
    graph.parse(str(path), format="turtle")
    generation = detect_generation(graph)
    if generation == GEN_A:
        return parse_gen_a(graph)
    if generation == GEN_B:
        return parse_gen_b(graph)
    return []
