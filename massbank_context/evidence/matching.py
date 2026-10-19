"""群間比較のID対応に従って照合し、候補を保存する。"""
import hashlib
import math
import sqlite3
from dataclasses import asdict
from massbank_context.core import config
from massbank_context.core.evidence_ids import stable_id, canonical_json
from massbank_context.core.evidence_models import Candidate
from massbank_context.index.writer import read_meta
from massbank_context.spectral.matcher import match_spectrum_result, fragment_tolerance


def select_ids(all_ids, requested):
    ids = sorted(all_ids if requested is None else set(requested))
    if set(ids) - set(all_ids):
        raise ValueError('未知のIDが指定されています')
    return ids


def match_features(store, analysis_id: str, feature_ids, *, top_k=5,
                   acquisition_class='LC-ESI-MS2', precursor_tolerance=.01):
    if type(top_k) is not int or not 1 <= top_k <= 20:
        raise ValueError('top_kは1〜20です')
    if not math.isfinite(precursor_tolerance) or precursor_tolerance < 0:
        raise ValueError('precursor_toleranceが不正です')
    bundle, base_revision, state = store.read_view(analysis_id)
    ids = select_ids([f.feature_id for f in bundle.features], feature_ids)
    previous = state.get('matches', {})
    statuses = dict(previous.get('features', {}))
    candidates = [c for c in previous.get('candidates', []) if c['feature_id'] not in ids]
    needed = any(s.feature_id in ids for s in bundle.spectrum_links)
    conn, index_hash, meta = None, None, {}
    try:
        if needed:
            if not config.INDEX_PATH.is_file():
                raise ValueError('missing_index: 照合用索引がありません')
            conn = sqlite3.connect(config.INDEX_PATH.resolve().as_uri() + '?mode=ro', uri=True)
            meta = read_meta(conn)
            if not meta or acquisition_class not in meta['acquisition_classes']:
                raise ValueError('incompatible_index: 測定層が索引にありません')
            with config.INDEX_PATH.open('rb') as stream:
                index_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
            if previous.get('index_hash') and previous['index_hash'] != index_hash:
                statuses, candidates = {}, []
        for fid in ids:
            links = [s for s in bundle.spectrum_links if s.feature_id == fid]
            status = dict(status='no_msms' if not links else 'matched', truncated=False,
                          n_total=0, boundary_tie_count=0)
            source_status = bundle.manifest.get('provenance', {}).get('msms_link_statuses', {}).get(fid, {})
            if not links and source_status.get('status') in ('precursor_mismatch', 'metadata_missing', 'out_of_range'):
                status.update(status='msms_mapping_failed', mapping_reason=source_status['status'])
            for link in links:
                spectrum = bundle.spectra[link.spectrum_file + '::' + link.spectrum_id]
                if not spectrum.peaks:
                    status['status'] = 'empty_spectrum'
                    continue
                result = match_spectrum_result(conn, spectrum, acquisition_class,
                                               precursor_tolerance, top_k)
                status['n_total'] += result.n_total
                status['truncated'] |= result.truncated
                status['boundary_tie_count'] += result.boundary_tie_count
                for rank, match in enumerate(result.candidates, 1):
                    provenance = dict(index_hash=index_hash, index_meta=meta,
                                      acquisition_class=acquisition_class,
                                      precursor_tolerance=precursor_tolerance,
                                      fragment_tolerance=fragment_tolerance(meta['bin_width']),
                                      spectrum_file=link.spectrum_file)
                    cid = stable_id('candidate', [analysis_id, fid, link.spectrum_file,
                                                  link.spectrum_id, match.accession,
                                                  canonical_json(provenance)])
                    candidates.append(asdict(Candidate(cid, fid, 'massbank', match.accession,
                                                       spectrum_id=link.spectrum_id, inchikey=match.inchikey,
                                                       rank=rank, scores=asdict(match), provenance=provenance)))
            statuses[fid] = status
        payload = dict(features=statuses, candidates=candidates,
                       index_hash=index_hash or previous.get('index_hash'))
        store.put_snapshot(analysis_id, 'matches', payload, expected_revision=base_revision)
        return dict(n_candidates=len(candidates), features=statuses,
                    evidence_revision=store.latest_revision(analysis_id))
    finally:
        if conn:
            conn.close()
