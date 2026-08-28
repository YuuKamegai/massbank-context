"""世代A の TTL を SpectrumRecord に変換する（spec §3.4）。

世代A の構造:
  <accession> mbo:has_peak [ mbo:mz ; mbo:intensity ; mbo:relative_intensity ]
              mbo:has_peak_annotations [ mbo:mz ; mbo:tentative_formula ]
              mbo:analytical_methods_and_conditions [ ... ]
              mbo:compound [ schema:inChIKey ; ... ]
              mbo:mass_spectram_meta_data [ mbo:precursor_mz_value ; ... ]

has_peak と has_peak_annotations は別ブロックなので m/z で突き合わせる。
アノテーションは全レコードにあるわけではない（サブミッター単位で全か無か）。
"""
from __future__ import annotations

import rdflib

from massbank_context.core.models import Peak, SpectrumRecord
from massbank_context.ttl.generation import MBO, SCHEMA

#: m/z を突き合わせる際の丸め桁。TTL 内の同一ピークは同一リテラルなので
#: 完全一致でよいが、指数表記の揺れを吸収するために丸めて鍵にする。
_MZ_KEY_DIGITS = 6


def _one(graph: rdflib.Graph, subject, predicate):
    for value in graph.objects(subject, predicate):
        return value
    return None


def _text(value) -> str | None:
    return None if value is None else str(value)


def _number(value) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _inchikey(value) -> str | None:
    """inchikey:XSDQ... という URI から末尾のキーだけを取り出す。"""
    if value is None:
        return None
    text = str(value)
    return text.rsplit("/", 1)[-1] if "/" in text else text


def parse_gen_a(graph: rdflib.Graph) -> list[SpectrumRecord]:
    """グラフに含まれる全レコードを返す(通常は 1 ファイル 1 レコード)。"""
    records: list[SpectrumRecord] = []
    for subject in graph.subjects(rdflib.RDF.type, MBO["MassSpectrum"]):
        records.append(_parse_one(graph, subject))
    records.sort(key=lambda r: r.accession)
    return records


def _parse_one(graph: rdflib.Graph, subject) -> SpectrumRecord:
    conditions = _one(graph, subject, MBO["analytical_methods_and_conditions"])
    compound = _one(graph, subject, MBO["compound"])
    meta = _one(graph, subject, MBO["mass_spectram_meta_data"])

    names = tuple(sorted(
        str(v) for v in graph.objects(compound, SCHEMA["name"])
    )) if compound is not None else ()

    formula_by_mz: dict[float, str] = {}
    for annotation in graph.objects(subject, MBO["has_peak_annotations"]):
        mz = _number(_one(graph, annotation, MBO["mz"]))
        formula = _text(_one(graph, annotation, MBO["tentative_formula"]))
        if mz is not None and formula:
            formula_by_mz[round(mz, _MZ_KEY_DIGITS)] = formula

    peaks: list[Peak] = []
    for node in graph.objects(subject, MBO["has_peak"]):
        mz = _number(_one(graph, node, MBO["mz"]))
        if mz is None:
            continue
        peaks.append(Peak(
            mz=mz,
            relative_intensity=_number(
                _one(graph, node, MBO["relative_intensity"])) or 0.0,
            formula=formula_by_mz.get(round(mz, _MZ_KEY_DIGITS)),
            neutral_loss=None,
        ))
    peaks.sort(key=lambda p: p.mz)

    license_value = _one(graph, subject, rdflib.URIRef(
        "http://purl.org/dc/terms/license"))

    return SpectrumRecord(
        accession=_text(_one(graph, subject, rdflib.URIRef(
            "http://purl.org/dc/terms/identifier"))) or str(subject).rsplit("/", 1)[-1],
        inchikey=_inchikey(_one(graph, compound, SCHEMA["inChIKey"]))
        if compound is not None else None,
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
        license=_text(license_value),
        molecular_formula=_text(_one(graph, compound, SCHEMA["molecularFormula"]))
        if compound is not None else None,
        exact_mass=_number(_one(graph, compound, MBO["ch_exact_mass"]))
        if compound is not None else None,
        compound_names=names,
        ch_compound_class=_text(_one(graph, compound, MBO["ch_compound_class"]))
        if compound is not None else None,
        peaks=tuple(peaks),
    )
