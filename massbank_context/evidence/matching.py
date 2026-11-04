"""群間比較のID対応に従って照合し、候補を保存する。"""
import math
from dataclasses import asdict
from massbank_context.core import config
from massbank_context.core.evidence_ids import stable_id, canonical_json
from massbank_context.core.evidence_models import Candidate
from massbank_context.evidence.selection import links_by_feature, mapping_failure, select_ids
from massbank_context.index.writer import index_fingerprint, open_index_readonly, read_meta
from massbank_context.spectral.matcher import match_spectrum_result, fragment_tolerance

# 照合を終えた状態。これ以外のリンク付き特徴量は照合待ちとして統計から外す。
TERMINAL_STATUSES = ('matched', 'empty_spectrum')


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
    selected = set(ids)
    candidates = [c for c in previous.get('candidates', []) if c['feature_id'] not in selected]
    links = links_by_feature(bundle)
    conn, fingerprint, meta = None, None, {}
    try:
        if any(fid in links for fid in ids):
            if not config.INDEX_PATH.is_file():
                raise ValueError('missing_index: 照合用索引がありません')
            conn = open_index_readonly(config.INDEX_PATH)
            meta = read_meta(conn)
            if not meta or acquisition_class not in meta['acquisition_classes']:
                raise ValueError('incompatible_index: 測定層が索引にありません')
            fingerprint = index_fingerprint(config.INDEX_PATH, meta)
            # 別の索引で得た候補と混ぜない。
            if previous.get('candidates') and previous.get('index_fingerprint') != fingerprint:
                statuses, candidates = {}, []
        provenance = dict(index_fingerprint=fingerprint, acquisition_class=acquisition_class,
                          precursor_tolerance=precursor_tolerance,
                          fragment_tolerance=fragment_tolerance(meta['bin_width']) if meta else None)
        for fid in ids:
            feature_links = links.get(fid, [])
            status = dict(status='no_msms' if not feature_links else 'matched', truncated=False,
                          n_total=0, boundary_tie_count=0, n_empty_spectra=0)
            reason = mapping_failure(bundle, fid)
            if not feature_links and reason:
                status.update(status='msms_mapping_failed', mapping_reason=reason)
            for link in feature_links:
                spectrum = bundle.spectra[link.spectrum_file + '::' + link.spectrum_id]
                if not spectrum.peaks:
                    status['n_empty_spectra'] += 1
                    continue
                result = match_spectrum_result(conn, spectrum, acquisition_class,
                                               precursor_tolerance, top_k)
                status['n_total'] += result.n_total
                status['truncated'] |= result.truncated
                status['boundary_tie_count'] += result.boundary_tie_count
                source = dict(provenance, spectrum_file=link.spectrum_file)
                for rank, match in enumerate(result.candidates, 1):
                    cid = stable_id('candidate', [analysis_id, fid, link.spectrum_file,
                                                  link.spectrum_id, match.accession,
                                                  canonical_json(source)])
                    candidates.append(asdict(Candidate(cid, fid, 'massbank', match.accession,
                                                       spectrum_id=link.spectrum_id, inchikey=match.inchikey,
                                                       rank=rank, scores=asdict(match), provenance=source)))
            # 一部のスペクトルだけが空なら照合は成立している。全て空のときだけ空扱い。
            if feature_links and status['n_empty_spectra'] == len(feature_links):
                status['status'] = 'empty_spectrum'
            statuses[fid] = status
        payload = dict(features=statuses, candidates=candidates,
                       index_fingerprint=fingerprint or previous.get('index_fingerprint'),
                       index_meta=meta or previous.get('index_meta', {}))
        store.put_snapshot(analysis_id, 'matches', payload, expected_revision=base_revision)
        return dict(n_candidates=len(candidates), features=statuses,
                    evidence_revision=store.latest_revision(analysis_id))
    finally:
        if conn:
            conn.close()
