"""TTL のスキーマ世代を判定する（spec §3.3）。

MassBank-RDF には 2 世代のスキーマが存在する。世代A（draft_ttl/ の実体）は
ピーク値をリテラルで直付けし、世代B（SPARQL/ の例が前提とするもの）は
SIO で具象化したノードに sio:SIO_000300 で値を付ける。判定は述語の存在で行う。
"""
from __future__ import annotations

import rdflib

MBO = rdflib.Namespace("http://www.massbank.jp/ontology/")
SIO = rdflib.Namespace("http://semanticscience.org/resource/")
SCHEMA = rdflib.Namespace("https://schema.org/")

#: 世代B で値を保持する述語。
SIO_HAS_VALUE = SIO["SIO_000300"]

GEN_A = "gen_a"
GEN_B = "gen_b"
UNKNOWN = "unknown"


def detect_generation(graph: rdflib.Graph) -> str:
    """グラフの世代を返す。

    世代B は mz を rdf:type mbo:mz のノードとして持つ。世代A は
    mbo:mz をリテラル述語として使う。両者は排他的なので type 文の有無で分かれる。
    """
    if (None, rdflib.RDF.type, MBO["mz"]) in graph:
        return GEN_B
    for _ in graph.subject_objects(MBO["has_peak"]):
        return GEN_A
    return UNKNOWN
