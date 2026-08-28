"""世代B の TTL を SpectrumRecord に変換する（spec §3.3）。

世代B は値を SIO で具象化する:
  <accession> mbo:has_peak ?peak
  ?peak  sio:SIO_000255 ?annotation
  ?annotation sio:SIO_000255 ?content
  ?content a mbo:mz | mbo:relative_intensity | mbo:intensity | mbo:neutral_loss
  ?content sio:SIO_000300 ?value

世代A と違い、中性損失が既に入っているので算出しない。
"""
from __future__ import annotations

import rdflib

from massbank_context.core.models import Peak, SpectrumRecord
from massbank_context.ttl._literals import _number, _one, _text
from massbank_context.ttl.generation import MBO, SCHEMA, SIO, SIO_HAS_VALUE

SIO_HAS_ATTRIBUTE = SIO["SIO_000255"]


def _typed_value(graph: rdflib.Graph, annotation, type_uri):
    """?annotation sio:SIO_000255 [ a <type_uri> ; sio:SIO_000300 ?v ] の ?v。"""
    for content in graph.objects(annotation, SIO_HAS_ATTRIBUTE):
        if (content, rdflib.RDF.type, type_uri) in graph:
            return _one(graph, content, SIO_HAS_VALUE)
    return None


def _inchikey(graph: rdflib.Graph, compound) -> str | None:
    """世代B では inChIKey も具象化ノードになりうる。両形を受ける。"""
    node = _one(graph, compound, SCHEMA["inChIKey"])
    if node is None:
        return None
    if isinstance(node, rdflib.Literal):
        return str(node)
    value = _one(graph, node, SIO_HAS_VALUE)
    if value is not None:
        return str(value)
    text = str(node)
    return text.rsplit("/", 1)[-1] if "/" in text else text


def parse_gen_b(graph: rdflib.Graph) -> list[SpectrumRecord]:
    records: list[SpectrumRecord] = []
    for subject in graph.subjects(rdflib.RDF.type, MBO["MassSpectrum"]):
        records.append(_parse_one(graph, subject))
    records.sort(key=lambda r: r.accession)
    return records


def _parse_one(graph: rdflib.Graph, subject) -> SpectrumRecord:
    conditions = _one(graph, subject, MBO["analytical_methods_and_conditions"])
    compound = _one(graph, subject, MBO["compound"])
    meta = _one(graph, subject, MBO["mass_spectram_meta_data"])

    peaks: list[Peak] = []
    for peak_node in graph.objects(subject, MBO["has_peak"]):
        for annotation in graph.objects(peak_node, SIO_HAS_ATTRIBUTE):
            mz = _number(_typed_value(graph, annotation, MBO["mz"]))
            if mz is None:
                continue
            peaks.append(Peak(
                mz=mz,
                relative_intensity=_number(_typed_value(
                    graph, annotation, MBO["relative_intensity"])) or 0.0,
                formula=_text(_typed_value(
                    graph, annotation, MBO["tentative_formula"])),
                neutral_loss=_number(_typed_value(
                    graph, annotation, MBO["neutral_loss"])),
            ))
    peaks.sort(key=lambda p: p.mz)

    return SpectrumRecord(
        accession=_text(_one(graph, subject, rdflib.URIRef(
            "http://purl.org/dc/terms/identifier"))) or str(subject).rsplit("/", 1)[-1],
        inchikey=_inchikey(graph, compound) if compound is not None else None,
        precursor_mz=_number(_one(graph, meta, MBO["precursor_mz_value"]))
        if meta is not None else None,
        precursor_type=_text(_one(graph, meta, MBO["precursor_type_value"]))
        if meta is not None else None,
        polarity=_text(_one(graph, conditions, MBO["ion_mode"]))
        if conditions is not None else None,
        instrument_type=_text(_one(graph, conditions, MBO["instrument_type"]))
        if conditions is not None else None,
        ms_type=_text(_one(graph, conditions, MBO["ms_type"]))
        if conditions is not None else None,
        resolution_value=_number(_one(graph, conditions, MBO["resolution"]))
        if conditions is not None else None,
        collision_energy=_text(_one(graph, conditions, MBO["collision_energy"]))
        if conditions is not None else None,
        splash=_text(_one(graph, subject, MBO["pk_splash"])),
        license=_text(_one(graph, subject, rdflib.URIRef(
            "http://purl.org/dc/terms/license"))),
        molecular_formula=_text(_one(graph, compound, SCHEMA["molecularFormula"]))
        if compound is not None else None,
        exact_mass=_number(_one(graph, compound, MBO["ch_exact_mass"]))
        if compound is not None else None,
        compound_names=tuple(sorted(
            str(v) for v in graph.objects(compound, SCHEMA["name"])
        )) if compound is not None else (),
        ch_compound_class=_text(_one(graph, compound, MBO["ch_compound_class"]))
        if compound is not None else None,
        peaks=tuple(peaks),
    )
