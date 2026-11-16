"""走査: 質量・MS/MS の共有・出方の 3 つの筋から、特徴と entry の組に水準を付ける。

全特徴を TSV に書き、戻り値は要約だけにする。行はスナップショットに残し、
``get_screen_rows`` で引く。名前は照合キーにしない（InChIKey と組成式だけを使う）。
"""
from __future__ import annotations

import base64
import binascii
import csv
import hashlib
import io
import json
import math
import os
from collections import Counter
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Sequence

from massbank_context.core import config
from massbank_context.core.evidence_ids import canonical_json
from massbank_context.core.evidence_models import EvidenceBundle, Feature
from massbank_context.core.models import QuerySpectrum
from massbank_context.core.serialization import as_tsv, round_floats
from massbank_context.evidence.store import EvidenceStore
from massbank_context.experiment.bundle import (
    POLARITIES, parse_intensities, parse_samples, role_path)
from massbank_context.xenobiotic.mass import (
    MassRelation, adduct_recognized, find_mass_relations, validate_ppm_bands, validate_ppm_offset)
from massbank_context.spectral.matcher import BAND_PASS
from massbank_context.xenobiotic.msms import (
    ForwardHit, SharedResult, forward_rank, open_index, reference_bins, shared_evidence)
from massbank_context.xenobiotic.pattern import (
    SampleSides, detection_counts, group_bias, parent_linked, parent_rho, sample_sides)
from massbank_context.xenobiotic.registry import XenobioticError, XenoEntry, load_registration
from massbank_context.xenobiotic.tables import (
    ADDUCT_TABLE_VERSION, CALIBRATION, SHIFT_TABLE_VERSION, THRESHOLDS, THRESHOLDS_VERSION)

__all__ = ["SCREEN_COLUMNS", "LEVELS", "PairEvidence", "assign_level", "resolve_polarity",
           "resolve_test_side", "run_screen", "get_screen_rows"]

SCREEN_COLUMNS = ['feature_id', 'mz', 'rt', 'log2fc', 'q_value', 'significant', 'level', 'entry_id',
                  'role', 'relation', 'theoretical_mz', 'ppm', 'msms_shared_score', 'msms_shared_n',
                  'msms_shared_formulas', 'forward_rank', 'detected_test', 'detected_reference',
                  'rho_parent', 'notes', 'coeluting']
LEVELS = ('A', 'B', 'C', 'D')
#: 共溶出とみなす RT の差の既定（分）。印だけで水準には使わない。
DEFAULT_COELUTION_RT_TOL = 0.1
_RT_EPS = 1e-9  # 浮動小数の境目（1.3 − 1.2 など）を窓の内に入れる

SNAPSHOT_KIND = 'xenobiotic_screen'
MAX_ROWS_LIMIT = 200
LOW_SCORE_SOURCE = 'msdial_sml'
LOW_SCORE_PREFIX = 'low score:'
_FEATURE_ID = SCREEN_COLUMNS.index('feature_id')
_LEVEL = SCREEN_COLUMNS.index('level')
_DETECTED_TEST = SCREEN_COLUMNS.index('detected_test')


@dataclass(frozen=True)
class PairEvidence:
    """特徴と entry の組の筋。``bias`` は判定できなければ None。"""
    mass_parent: bool
    mass_any: bool
    msms_hit: bool
    forward_rank: int | None
    forward_band: str | None
    bias: bool | None
    parent_link: bool
    mass_plain: bool = True  # 損失の無い質量の一致が 1 つでもあるか
    parent_codetected: bool | None = None  # 親の特徴と同じ試料で出ているか（判定できなければ None）


def assign_level(pair: PairEvidence, received_by: str) -> str | None:
    """組の水準（spec §5）。both は群の偏りを使わず、親との連動だけを出方とする。

    順位 1 でも band が PASS でなければ A にしない。PRECURSOR_ONLY は精密質量が合っただけで、
    全候補が得点 0 なら順位は accession の文字列順で決まる。

    損失つきの質量の一致だけの組は、群の偏りでは C にしない。r2 でおとりの C を増やしたのは
    この組だった。MS/MS の共有があれば先に B になるので、ここで効くのは親との連動だけ。

    質量差のある一致だけの組は、MS/MS の共有があっても、親の特徴と同じ試料で出ていなければ
    B にしない（spec §5 の改訂、r4 の偽陽性 4 件）。共起が判定できないときも B にしない。
    """
    if pair.mass_parent and pair.forward_rank == 1 and pair.forward_band == BAND_PASS:
        return 'A'
    if pair.mass_any and pair.msms_hit and (pair.mass_parent or pair.parent_codetected is True):
        return 'B'
    strand = pair.mass_any or pair.msms_hit
    pattern = pair.parent_link if received_by == 'both' else pair.bias is True or pair.parent_link
    loss_only = pair.mass_any and not pair.mass_plain
    if strand and pattern and (not loss_only or pair.parent_link):
        return 'C'
    return 'D' if strand else None


def resolve_polarity(manifest: dict, polarity: str | None) -> str:
    """manifest と引数から極性（小文字）を決める。食い違いは黙って選ばず止める。"""
    if polarity is not None and polarity not in POLARITIES:
        raise XenobioticError(f"polarity は {POLARITIES} のいずれか: {polarity!r}")
    declared = manifest['comparison'].get('polarity')
    if declared is None:
        if polarity is None:
            raise XenobioticError('このバンドルには極性が無いので polarity を指定してください',
                                  code='polarity_required')
        return polarity
    if polarity is not None and polarity != declared:
        raise XenobioticError(f"polarity がバンドルと違います: {polarity} / {declared}",
                              code='polarity_mismatch')
    return declared


def resolve_test_side(manifest: dict, test_group: str) -> str:
    """test_group が comparison のどちらの群か（'a' / 'b'）。名前から推測しない。"""
    comparison = manifest['comparison']
    for side in ('a', 'b'):
        if test_group == comparison[f'group_{side}']:
            return side
    raise XenobioticError(
        f"test_group は {comparison['group_a']} か {comparison['group_b']}: {test_group!r}",
        code='unknown_test_group')


def _validated_bands(ppm, max_shifts, ppm_bands) -> list[dict]:
    if isinstance(ppm, bool) or not isinstance(ppm, (int, float)) or not math.isfinite(ppm) or ppm <= 0:
        raise XenobioticError(f"ppm は正の有限数: {ppm!r}")
    if type(max_shifts) is not int or not 0 <= max_shifts <= 2:
        raise XenobioticError(f"max_shifts は 0〜2 の整数: {max_shifts!r}")
    try:
        return validate_ppm_bands(ppm_bands)
    except (ValueError, TypeError, AttributeError) as exc:
        raise XenobioticError(f"ppm_bands が不正: {exc}") from exc


def _validated_offset(ppm_offset) -> dict | None:
    try:
        return validate_ppm_offset(ppm_offset)
    except ValueError as exc:
        raise XenobioticError(f"ppm_offset が不正: {exc}") from exc


def _validated_rt_tol(value) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise XenobioticError(f"coelution_rt_tol は正の有限数（分）: {value!r}")
    return float(value)


def _coelution(features: Sequence[Feature], relations: dict,
               rt_tol: float) -> dict[tuple, list[tuple[str, str]]]:
    """同じ分子の別のイオンの形に一致し、同じ RT に出る別の特徴（共溶出の印。spec §6）。

    同じ分子とは、同じ entry・同じ質量差の組。別の形とは、付加体か損失が違うこと
    （同じ形に割れたスポットは裏付けにならない）。RT の差が ``rt_tol`` 分以内なら共溶出。
    鍵は (feature_id, entry_id, shifts, ion)、値は (相手の特徴, 相手のイオン) を特徴の順に並べたもの。
    水準 A〜D には使わない。印だけで、較正（spec §10）の結果は変えない。
    """
    rt = {f.feature_id: f.rt for f in features}
    order = {f.feature_id: i for i, f in enumerate(features)}
    groups: dict[tuple, list[tuple[str, str]]] = {}
    for (fid, eid), found in relations.items():
        if rt.get(fid) is None:
            continue
        for relation in found:
            groups.setdefault((eid, relation.shifts), []).append((fid, relation.ion))
    result: dict[tuple, list[tuple[str, str]]] = {}
    for (eid, shifts), members in groups.items():
        for fid, ion in members:
            partners = {(other, other_ion) for other, other_ion in members
                        if other != fid and other_ion != ion
                        and abs(rt[fid] - rt[other]) <= rt_tol + _RT_EPS}
            if partners:
                result[(fid, eid, shifts, ion)] = sorted(partners, key=lambda p: (order[p[0]], p[1]))
    return result


def _coelution_cell(coelution: dict, feature_id: str, entry_id: str,
                    relation: MassRelation) -> str | None:
    partners = coelution.get((feature_id, entry_id, relation.shifts, relation.ion))
    return ';'.join(f'{fid}:{ion}' for fid, ion in partners) if partners else None


def _low_score_adducts(bundle: EvidenceBundle) -> frozenset[str]:
    """付加体の記載が低得点の注釈の候補でしかない特徴。"""
    return frozenset(c.feature_id for c in bundle.annotations
                     if c.source == LOW_SCORE_SOURCE and c.name.startswith(LOW_SCORE_PREFIX))


def _feature_spectra(bundle: EvidenceBundle) -> dict[str, list[QuerySpectrum]]:
    """特徴に結ばれたスペクトル。bundle.spectra に無いリンクは無いものとして扱う。"""
    result: dict[str, list[QuerySpectrum]] = {}
    for link in bundle.spectrum_links:
        spectrum = bundle.spectra.get(f'{link.spectrum_file}::{link.spectrum_id}')
        if spectrum is not None:
            result.setdefault(link.feature_id, []).append(spectrum)
    return result


@dataclass(frozen=True)
class _MsmsStrand:
    unavailable: dict[str, str]  # entry_id -> 理由
    shared: dict[tuple[str, str], SharedResult]
    hits: frozenset[tuple[str, str]]
    ranks: dict[tuple[str, str], ForwardHit | None]
    all_reason: str | None = None  # 筋全体が使えない理由


def _links_unavailable(bundle: EvidenceBundle) -> str | None:
    """links が空で対応づけも確かめられていないとき、MS/MS の筋全体が使えない理由。

    ``verified_alignment`` で空なら、対応づけた結果 MS/MS が無かったという事実なので ``no_msms``。
    それ以外の空は MS/MS の有無が分からないので、全特徴を ``no_msms`` と言わない。
    """
    status = bundle.manifest.get('provenance', {}).get('spectrum_links_status')
    if bundle.spectrum_links or status == 'verified_alignment':
        return None
    return 'no_spectrum_links'


def _msms_strand(conn, meta: dict, entries: Sequence[XenoEntry], polarity: str,
                 spectra: dict[str, list[QuerySpectrum]], relations: dict,
                 all_reason: str | None) -> _MsmsStrand:
    bin_width = meta['bin_width']
    unavailable: dict[str, str] = {}
    refs = {}
    for entry in entries:
        ref = None
        if entry.inchikey is not None:
            ref = reference_bins(conn, entry.inchikey, polarity.upper(), bin_width,
                                 THRESHOLDS['msms_top_n'])
        if ref is None:
            unavailable[entry.entry_id] = 'no_reference_spectra' if entry.inchikey else 'no_inchikey'
        else:
            refs[entry.entry_id] = ref
    shared: dict[tuple[str, str], SharedResult] = {}
    for feature_id, group in spectra.items():
        for entry_id, ref in refs.items():
            result = shared_evidence(ref, group, bin_width)
            if result is not None:
                shared[(feature_id, entry_id)] = result
    hits = frozenset(key for key, result in shared.items()
                     if result.score >= THRESHOLDS['msms_shared_min'])
    # 全組で回すと照合が特徴 × entry 回になるので、筋のある組だけで順位を見る。
    ranks = {key: forward_rank(conn, spectra[key[0]], refs[key[1]].inchikey)
             for key in hits | set(relations) if key[0] in spectra and key[1] in refs}
    return _MsmsStrand(unavailable, shared, hits, ranks, all_reason)


@dataclass(frozen=True)
class _PatternStrand:
    """出方の筋。``reason`` があれば筋全体が使えない。"""
    reason: str | None
    sides: SampleSides | None = None
    matrix: dict | None = None
    counts: dict | None = None


def _pattern_strand(store: EvidenceStore, analysis_id: str, bundle: EvidenceBundle,
                    test_side: str) -> _PatternStrand:
    manifest = bundle.manifest
    intensities = role_path(manifest, 'intensities')
    if intensities is None:
        return _PatternStrand('no_intensities')

    def read(path: str) -> str:
        return store.read_source(analysis_id, path).decode('utf-8-sig')

    samples = parse_samples(read(role_path(manifest, 'samples')))
    ids = {f.feature_id for f in bundle.features}
    matrix = parse_intensities(read(intensities), ids, samples)
    sides = sample_sides(samples, test_side)
    counts = {fid: detection_counts(matrix[fid], sides) for fid in ids}
    # 検出が全部不明なら「未検出」と数えず、値だけの ρ からも水準を上げない。
    unknown = all(pair == (None, None) for pair in counts.values())
    return _PatternStrand('detection_unknown' if unknown else None, sides, matrix, counts)


@dataclass(frozen=True)
class _Pair:
    entry: XenoEntry
    relations: tuple[MassRelation, ...]
    evidence: PairEvidence
    level: str
    shared: SharedResult | None
    rho: float | None
    parent_feature: str | None
    codetection: tuple[int, int] | None = None  # (k, n)。質量差のある一致だけで MS/MS の共有がある組だけ


def _parent_link(feature_id: str, entry: XenoEntry, relations: Sequence[MassRelation],
                 a_features: dict[str, list[str]],
                 pattern: _PatternStrand) -> tuple[bool, float, str] | None:
    """親の水準 A の特徴との連動 (連動か, ρ, 親の特徴)。比べる相手が無ければ None。"""
    if pattern.reason is not None:
        return None
    if entry.parent_entry_id is not None:
        parent_id = entry.parent_entry_id
    elif any(r.shifts for r in relations):
        parent_id = entry.entry_id
    else:
        return None
    scored = []
    for other in a_features.get(parent_id, ()):
        if other == feature_id:
            continue  # 同じ特徴との ρ は常に 1 で、偽の連動になる
        rho, n = parent_rho(pattern.matrix[feature_id], pattern.matrix[other], pattern.sides)
        if rho is not None:
            scored.append((parent_linked(rho, n, THRESHOLDS), rho, other))
    # n の足りない ρ=1 が、連動の成り立つ相手を隠さないよう、連動を先に見る。
    return max(scored, key=lambda item: item[:2], default=None)


def _parent_codetection(feature_id: str, entry: XenoEntry, entries_by_id: dict[str, XenoEntry],
                        plain_features: dict[str, list[str]],
                        pattern: _PatternStrand) -> tuple[int, int] | None:
    """(親の特徴と同じ試料で出た数, この特徴が出た試料の数)。出方が使えなければ None。

    親の特徴は、entry かその祖先（parent_entry_id をたどる）と質量差も損失も無く一致する、
    この特徴以外の特徴（水準は問わない）。試料は test と reference の両方を数え、
    この特徴の detected が True の試料だけを分母にする。
    """
    if pattern.reason is not None:
        return None
    lineage, seen = [], set()
    current: XenoEntry | None = entry
    while current is not None and current.entry_id not in seen:
        seen.add(current.entry_id)
        lineage.append(current.entry_id)
        parent_id = current.parent_entry_id
        current = entries_by_id.get(parent_id) if parent_id is not None else None
    parents = {fid for eid in lineage for fid in plain_features.get(eid, ()) if fid != feature_id}
    row = pattern.matrix[feature_id]
    k = n = 0
    for sid in pattern.sides.test + pattern.sides.reference:
        if sid not in row or row[sid][1] is not True:
            continue
        n += 1
        if any(sid in pattern.matrix[p] and pattern.matrix[p][sid][1] is True for p in parents):
            k += 1
    return k, n


def _first_evidence(found: Sequence[MassRelation], is_hit: bool,
                    forward: ForwardHit | None) -> PairEvidence:
    return PairEvidence(any(not r.shifts for r in found), bool(found), is_hit,
                        forward.rank if forward else None, forward.band if forward else None,
                        None, False, any(not r.has_loss for r in found))


def _judge_pairs(features: Sequence[Feature], entries: Sequence[XenoEntry], relations: dict,
                 msms: _MsmsStrand, pattern: _PatternStrand) -> dict[str, list[_Pair]]:
    """特徴ごとの組を水準の強い順 → entry_id の順で返す。

    一巡目で全組の A を決める（A は質量と順位と band だけで決まる）。二巡目で親の A の特徴と
    比べて連動を見て、残りの水準を決める。
    """
    first: dict[tuple[str, str], PairEvidence] = {}
    a_features: dict[str, list[str]] = {}
    plain_features: dict[str, list[str]] = {}
    entries_by_id = {entry.entry_id: entry for entry in entries}
    for feature in features:
        for entry in entries:
            key = (feature.feature_id, entry.entry_id)
            found = relations.get(key, [])
            if any(not r.shifts and not r.has_loss for r in found):
                plain_features.setdefault(entry.entry_id, []).append(feature.feature_id)
            if not found and key not in msms.hits:
                continue
            first[key] = _first_evidence(found, key in msms.hits, msms.ranks.get(key))
            if assign_level(first[key], entry.received_by) == 'A':
                a_features.setdefault(entry.entry_id, []).append(feature.feature_id)

    result: dict[str, list[_Pair]] = {}
    for feature in features:
        pairs = []
        for entry in entries:
            key = (feature.feature_id, entry.entry_id)
            if key not in first:
                continue
            found = tuple(relations.get(key, ()))
            link = _parent_link(feature.feature_id, entry, found, a_features, pattern)
            bias = None
            if pattern.reason is None:
                bias = group_bias(pattern.counts[feature.feature_id], pattern.sides,
                                  entry.received_by, THRESHOLDS)
            codetection = codetected = None
            base = first[key]
            if base.mass_any and not base.mass_parent and base.msms_hit:
                codetection = _parent_codetection(feature.feature_id, entry, entries_by_id,
                                                  plain_features, pattern)
                if codetection is not None:
                    k, n = codetection
                    codetected = n > 0 and k / n > THRESHOLDS['codetect_min_fraction']
            evidence = replace(base, bias=bias, parent_link=bool(link and link[0]),
                               parent_codetected=codetected)
            pairs.append(_Pair(entry, found, evidence, assign_level(evidence, entry.received_by),
                               msms.shared.get(key), link[1] if link else None,
                               link[2] if link else None, codetection))
        pairs.sort(key=lambda p: (LEVELS.index(p.level), p.entry.entry_id))
        result[feature.feature_id] = pairs
    return result


def _msms_state(feature_id: str, spectra: dict[str, list[QuerySpectrum]]) -> str | None:
    group = spectra.get(feature_id)
    if not group:
        return 'no_msms'
    return None if any(s.peaks for s in group) else 'empty_spectrum'


def _notes(feature: Feature, msms_state: str | None, entry_reason: str | None,
           pattern_reason: str | None, parent_feature: str | None,
           forward_band: str | None = None, n_dropped: int = 0,
           low_score: bool = False, codetection: tuple[int, int] | None = None) -> str:
    notes = []
    if feature.mz is None:
        notes.append('no_mz')
    if msms_state:
        notes.append(msms_state)
    if entry_reason:
        notes.append(f'msms_unavailable:{entry_reason}')
    if pattern_reason:
        notes.append(f'pattern_unavailable:{pattern_reason}')
    if parent_feature:
        notes.append(f'parent_feature:{parent_feature}')
    if forward_band:
        notes.append(f'forward_band:{forward_band}')
    if codetection is not None:
        notes.append(f'parent_codetected:{codetection[0]}/{codetection[1]}')
    if n_dropped:
        notes.append(f'adduct_conflict_dropped:{n_dropped}')
    if low_score and feature.mz is not None and feature.adduct:
        notes.append('adduct_low_score')
    if feature.mz is not None and feature.adduct and not adduct_recognized(feature.adduct):
        notes.append('adduct_unrecognized')
    return ';'.join(notes)


def _build_rows(features: Sequence[Feature], pairs: dict[str, list[_Pair]],
                spectra: dict[str, list[QuerySpectrum]], msms: _MsmsStrand,
                pattern: _PatternStrand, dropped: Counter,
                low_score: frozenset[str], coelution: dict) -> list[list]:
    rows: list[list] = []
    dropped_by_feature: Counter = Counter()
    for (fid, _entry_id), n in dropped.items():
        dropped_by_feature[fid] += n
    for feature in features:
        fid = feature.feature_id
        head = [fid, round_floats(feature.mz), round_floats(feature.rt),
                round_floats(feature.log2fc), round_floats(feature.q_value, keep_small=True),
                'true' if feature.significant else 'false']
        state = None if msms.all_reason else _msms_state(fid, spectra)
        counts = pattern.counts[fid] if pattern.counts else (None, None)
        reason = pattern.reason or ('detection_unknown' if None in counts else None)
        if not pairs[fid]:
            # 検出の例数は特徴の値で entry に依らない。水準の無い特徴の行にも入れる。
            blank = _DETECTED_TEST - len(head)
            rows.append(head + [None] * blank + [counts[0], counts[1], None]
                        + [_notes(feature, state, msms.all_reason, reason, None,
                                  n_dropped=dropped_by_feature[fid],
                                  low_score=fid in low_score), None])
            continue
        for pair in pairs[fid]:
            entry, shared = pair.entry, pair.shared
            tail = [None if shared is None else round_floats(shared.score, 4),
                    None if shared is None else shared.n_shared,
                    None if shared is None else ';'.join(shared.formulas),
                    pair.evidence.forward_rank, counts[0], counts[1], round_floats(pair.rho),
                    _notes(feature, state, msms.all_reason or msms.unavailable.get(entry.entry_id),
                           reason,
                           pair.parent_feature, pair.evidence.forward_band,
                           dropped[(fid, entry.entry_id)], fid in low_score,
                           pair.codetection)]
            start = [pair.level, entry.entry_id, entry.role]
            for relation in pair.relations:
                rows.append(head + start + [relation.relation, round_floats(relation.theoretical_mz),
                                            round_floats(relation.ppm, 4)] + tail
                            + [_coelution_cell(coelution, fid, entry.entry_id, relation)])
            if not pair.relations:
                rows.append(head + start + ['msms_only', None, None] + tail + [None])
    return rows


def _coelution_summary(features: Sequence[Feature], pairs: dict[str, list[_Pair]],
                       coelution: dict, rt_tol: float) -> dict:
    """共溶出の相手がいる特徴の数と、その組の水準（特徴ごとに最も強いもの）。"""
    levels = dict.fromkeys(LEVELS, 0)
    for feature in features:
        for pair in pairs[feature.feature_id]:   # 水準の強い順に並んでいる
            if any((feature.feature_id, pair.entry.entry_id, r.shifts, r.ion) in coelution
                   for r in pair.relations):
                levels[pair.level] += 1
                break
    return {'rt_tol': rt_tol, 'n_features': sum(levels.values()), 'levels': levels}


def _summarize(bundle: EvidenceBundle, entries: Sequence[XenoEntry],
               pairs: dict[str, list[_Pair]], spectra: dict[str, list[QuerySpectrum]],
               msms: _MsmsStrand, pattern: _PatternStrand,
               mass_counts: Counter, n_rows: int) -> dict:
    features, manifest = bundle.features, bundle.manifest
    levels = dict.fromkeys((*LEVELS, 'none'), 0)
    significant = dict(levels)
    per_entry = {e.entry_id: {'A': 0, 'B': 0} for e in entries}
    for feature in features:
        found = pairs[feature.feature_id]
        best = found[0].level if found else 'none'
        levels[best] += 1
        significant[best] += feature.significant
        for pair in found:
            if pair.level in ('A', 'B'):
                per_entry[pair.entry.entry_id][pair.level] += 1
    states = Counter() if msms.all_reason else Counter(
        _msms_state(f.feature_id, spectra) for f in features)
    return {
        'calibration': dict(CALIBRATION),
        'intensities_stage': manifest.get('intensities_stage'),
        'n_features': len(features),
        'levels': levels,
        'significant_levels': significant,
        'per_entry': per_entry,
        'unavailable': {
            'msms': msms.unavailable,
            # 使えないなら A・B は原理的に 0。「薬剤由来なし」に読ませない。
            'msms_all': msms.all_reason,
            'spectrum_links_status': manifest.get('provenance', {}).get('spectrum_links_status'),
            'n_features_without_msms': None if msms.all_reason
            else states['no_msms'] + states['empty_spectrum'],
            'pattern': pattern.reason,
            'detection_unavailable_reason': manifest['population'].get('detection_unavailable_reason'),
        },
        'counts': {key: mass_counts[key]
                   for key in ('adduct_conflict_dropped', 'adduct_low_score_ignored',
                               'adduct_unrecognized', 'no_mz')}
        | {key: states[key] for key in ('no_msms', 'empty_spectrum')},
        'n_rows': n_rows,
    }


def _tsv_bytes(rows: list[list]) -> bytes:
    buffer = io.StringIO(newline='')
    writer = csv.writer(buffer, delimiter='\t', lineterminator='\n')
    writer.writerow(SCREEN_COLUMNS)
    writer.writerows(rows)
    return buffer.getvalue().encode('utf-8')


def _publish(store: EvidenceStore, analysis_id: str, payload: dict, target: Path, data: bytes,
             base_revision: str) -> str:
    """TSV とスナップショットを揃えて残す。どの段で失敗しても前の状態に戻す。

    ファイルの差し替えを全部済ませてから確定する。前の TSV を別のアプリが開いていて
    動かせないときは、スナップショットを進める前にここで止まる。
    """
    if target.is_dir():
        raise XenobioticError(f"output_path がディレクトリです: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + '.tmp')
    backup = target.with_name(target.name + '.bak')
    revision = None
    moved_aside = placed = False
    try:
        temporary.write_bytes(data)
        if target.exists():
            os.replace(target, backup)
            moved_aside = True
        os.replace(temporary, target)
        placed = True
        revision = store.put_snapshot(analysis_id, SNAPSHOT_KIND, payload,
                                      expected_revision=base_revision)
    finally:
        temporary.unlink(missing_ok=True)
        if revision is not None:
            backup.unlink(missing_ok=True)
        elif moved_aside:
            os.replace(backup, target)
        elif placed:
            target.unlink(missing_ok=True)
    return revision


def run_screen(store: EvidenceStore, analysis_id: str, *, test_group: str, ppm: float,
               ppm_bands: list[dict] | None = None, polarity: str | None = None,
               max_shifts: int = 2, output_path: str | None = None,
               ppm_offset: dict | None = None,
               coelution_rt_tol: float = DEFAULT_COELUTION_RT_TOL) -> dict:
    """走査して全特徴の TSV とスナップショットを書き、要約を返す（行は載せない）。

    ``ppm_offset`` は質量の窓の中心のずれ（ms-data-parser と同じ形）。渡したときだけ
    ``params`` に残す（渡さなければ記録は前と同じ）。``coelution_rt_tol`` は共溶出の印の
    RT の幅（分）で、水準は変えない。
    """
    bundle, base_revision, _ = store.read_view(analysis_id)
    polarity = resolve_polarity(bundle.manifest, polarity)
    test_side = resolve_test_side(bundle.manifest, test_group)
    bands = _validated_bands(ppm, max_shifts, ppm_bands)
    offset = _validated_offset(ppm_offset)
    rt_tol = _validated_rt_tol(coelution_rt_tol)
    entries, registration_id = load_registration(store, analysis_id)

    features = bundle.features
    spectra = _feature_spectra(bundle)
    low_score = _low_score_adducts(bundle)
    relations, mass_counts, dropped = find_mass_relations(features, entries, polarity, ppm, bands,
                                                          max_shifts, low_score, offset)
    conn, meta = open_index()
    index = {'bin_width': meta['bin_width'], 'built_at': meta.get('built_at')}
    try:
        msms = _msms_strand(conn, meta, entries, polarity, spectra, relations,
                            _links_unavailable(bundle))
    finally:
        conn.close()
    pattern = _pattern_strand(store, analysis_id, bundle, test_side)
    pairs = _judge_pairs(features, entries, relations, msms, pattern)
    coelution = _coelution(features, relations, rt_tol)
    rows = _build_rows(features, pairs, spectra, msms, pattern, dropped, low_score, coelution)
    summary = _summarize(bundle, entries, pairs, spectra, msms, pattern, mass_counts, len(rows))
    summary['coelution'] = _coelution_summary(features, pairs, coelution, rt_tol)

    if output_path is None:
        # analysis_id の ':' は Windows のファイル名に使えない
        name = f"{analysis_id.split(':', 1)[1][:16]}_{polarity}.tsv"
        target = config.ARTIFACTS_PATH / 'xenobiotic' / name
    else:
        target = Path(output_path)
    data = _tsv_bytes(rows)
    sha256 = hashlib.sha256(data).hexdigest()
    versions = {'adduct_table': ADDUCT_TABLE_VERSION, 'shift_table': SHIFT_TABLE_VERSION,
                'thresholds': THRESHOLDS_VERSION}
    params = {'test_group': test_group, 'ppm': ppm, 'ppm_bands': bands, 'polarity': polarity,
              'max_shifts': max_shifts, 'coelution_rt_tol': rt_tol}
    if offset is not None:
        params['ppm_offset'] = offset
    payload = {'params': params, 'versions': versions, 'thresholds': dict(THRESHOLDS),
               'registration_id': registration_id, 'index': index, 'summary': summary,
               'columns': SCREEN_COLUMNS,
               'rows': rows, 'output_path': str(target), 'output_sha256': sha256}
    revision = _publish(store, analysis_id, payload, target, data, base_revision)
    return {'screen_revision': revision, 'registration_id': registration_id, 'polarity': polarity,
            'test_group': test_group, 'versions': versions, 'thresholds': dict(THRESHOLDS),
            'index': index, 'params': params, **summary, 'output_path': str(target), 'output_sha256': sha256,
            'next_tools': ['get_xenobiotic_screen']}


def _cursor_scope(revision: str, feature_ids: list[str] | None, level: str | None) -> dict:
    """カーソルが結び付く版と絞り込みの条件。特徴の指定は順に依らない。"""
    digest = None
    if feature_ids is not None:
        digest = hashlib.sha256(canonical_json(sorted(set(feature_ids))).encode()).hexdigest()
    return {'rev': revision, 'level': level, 'feature_ids': digest}


def _encode_cursor(scope: dict, offset: int) -> str:
    raw = canonical_json({**scope, 'offset': offset}).encode()
    return base64.urlsafe_b64encode(raw).decode()


def _cursor_offset(cursor: str, scope: dict) -> int:
    try:
        value = json.loads(base64.urlsafe_b64decode(cursor.encode()))
    except (ValueError, TypeError, UnicodeError, binascii.Error) as exc:
        raise ValueError('invalid_cursor') from exc
    if not isinstance(value, dict) or not isinstance(value.get('rev'), str):
        raise ValueError('invalid_cursor')
    offset = value.get('offset')
    if type(offset) is not int or offset < 0:
        raise ValueError('invalid_cursor')
    if any(value.get(key) != expected for key, expected in scope.items()):
        raise ValueError('stale_cursor: 走査か絞り込みの条件が変わりました')
    return offset


def get_screen_rows(store: EvidenceStore, analysis_id: str, *, feature_ids: list[str] | None = None,
                    level: str | None = None, cursor: str | None = None, limit: int = 50) -> dict:
    """最新の走査の行を引く。カーソルは走査の版と絞り込みに結び付くので、別の行は混ざらない。

    走査の後に登録し直していたら止めずに ``stale_registration`` で知らせる。
    """
    if type(limit) is not int or not 1 <= limit <= MAX_ROWS_LIMIT:
        raise XenobioticError(f"limit は 1〜{MAX_ROWS_LIMIT} の整数: {limit!r}")
    if level is not None and level not in LEVELS:
        raise XenobioticError(f"level は {LEVELS} のいずれか: {level!r}")
    revision = store.head_revision(analysis_id, SNAPSHOT_KIND)
    if revision is None:
        raise XenobioticError('走査の結果がありません。先に走査してください', code='missing_state')
    snapshot = store.get_snapshot(analysis_id, revision)
    rows = snapshot['rows']
    scope = _cursor_scope(revision, feature_ids, level)
    offset = _cursor_offset(cursor, scope) if cursor is not None else 0
    _, current_registration = load_registration(store, analysis_id)

    known = {row[_FEATURE_ID] for row in rows}
    unknown = [fid for fid in feature_ids or [] if fid not in known]
    wanted = None if feature_ids is None else set(feature_ids)
    matched = [row for row in rows
               if (wanted is None or row[_FEATURE_ID] in wanted)
               and (level is None or row[_LEVEL] == level)]
    end = offset + limit
    return {'screen_revision': revision, 'registration_id': snapshot['registration_id'],
            'stale_registration': snapshot['registration_id'] != current_registration,
            'n_matched': len(matched),
            # 列を足す前に保存した走査は、そのときの列名で出す（feature_id と level の位置は変えていない）
            'rows': as_tsv(snapshot.get('columns') or SCREEN_COLUMNS, matched[offset:end]),
            'next_cursor': _encode_cursor(scope, end) if end < len(matched) else None,
            'unknown_feature_ids': unknown}
