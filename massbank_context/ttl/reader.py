"""TTL ファイルを読む入口。世代を判定して適切なリーダへ振り分ける。"""
from __future__ import annotations

from pathlib import Path

import rdflib

from massbank_context.core.models import SpectrumRecord
from massbank_context.ttl.generation import GEN_A, GEN_B, detect_generation
from massbank_context.ttl.reader_gen_a import parse_gen_a
from massbank_context.ttl.reader_gen_b import parse_gen_b


def read_ttl_file_with_generation(
        path: str | Path) -> tuple[str, list[SpectrumRecord]]:
    """1 ファイルを (世代, レコード列) で読む。判定できない場合は空リスト。

    索引ビルドは数万ファイルを走査するため、1 件の異常でビルド全体を止めない。
    呼び出し側（CLI）が件数を集計して報告する。

    世代を返すのは、世代B の TTL 実物が入手できておらず（spec §11 確認事項 1）、
    もし実物の SIO ノードの入れ子が想定と違えば parse_gen_b が
    「メタデータは揃っているがピーク 0 本」のレコードを返すため。世代ごとの
    件数を出せば `gen_b: 12000 files, 0 peaks` として一目で分かる。
    """
    graph = rdflib.Graph()
    graph.parse(str(path), format="turtle")
    generation = detect_generation(graph)
    if generation == GEN_A:
        return generation, parse_gen_a(graph)
    if generation == GEN_B:
        return generation, parse_gen_b(graph)
    return generation, []


def read_ttl_file(path: str | Path) -> list[SpectrumRecord]:
    """1 ファイルのレコード列だけを読む。世代が要るなら上の関数を使う。"""
    return read_ttl_file_with_generation(path)[1]
