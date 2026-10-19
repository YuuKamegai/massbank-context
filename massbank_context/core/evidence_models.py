"""根拠受け渡しの型。上位の処理に依存しない。"""
from dataclasses import dataclass, field
from pathlib import Path
from massbank_context.core.models import QuerySpectrum


@dataclass(frozen=True)
class Feature:
    feature_id: str
    spot_id: str
    mz: float | None
    rt: float | None
    log2fc: float | None
    p_value: float | None
    q_value: float | None
    mean_a: float | None
    mean_b: float | None
    significant: bool


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    feature_id: str
    source: str
    source_id: str
    spectrum_id: str | None = None
    name: str = ''
    inchikey: str | None = None
    structure_ids: dict = field(default_factory=dict)
    rank: int | None = None
    scores: dict = field(default_factory=dict)
    status: str = 'active'
    same_structure_as: tuple[str, ...] = ()
    provenance: dict = field(default_factory=dict)
    curation: dict = field(default_factory=dict)


@dataclass(frozen=True)
class SpectrumLink:
    feature_id: str
    spectrum_file: str
    spectrum_id: str
    source_analysis_id: str
    mapping_evidence: str


@dataclass(frozen=True)
class EvidenceBundle:
    analysis_id: str
    manifest: dict
    features: tuple[Feature, ...]
    annotations: tuple[Candidate, ...]
    spectrum_links: tuple[SpectrumLink, ...]
    spectra: dict[str, QuerySpectrum]
    population_complete: bool
    manifest_path: Path
