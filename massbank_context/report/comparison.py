"""全体概要と、版を固定した特徴量別の根拠を提示する。"""
import base64
import json
from collections import Counter
from dataclasses import asdict
from massbank_context.core.evidence_ids import stable_id, canonical_json
from massbank_context.core.serialization import rows_tsv
from massbank_context.evidence.selection import candidate_dicts, select_ids
from massbank_context.report.comparison_population import build_population

MAX_PAGE_BYTES = 128 * 1024
EVIDENCE_PAGE_SIZE = 64
# 公開時にTSVへ落とす行の一覧。
TABLE_KEYS = ('features', 'candidates', 'evidence', 'lookups')


def evidence_rows(candidates, lookup_state) -> list[dict]:
    """候補に結び付いたPubChem根拠と索引メタデータ根拠（evidence_id順）。"""
    candidate_ids = {c['candidate_id'] for c in candidates}
    rows = [e for e in lookup_state.get('evidence', []) if e['candidate_id'] in candidate_ids]
    catalog = {}
    for record in lookup_state.get('catalog', []):
        catalog.setdefault((record['inchikey'], record['accession']), []).append(record)
    provenance = lookup_state.get('catalog_provenance', {})
    for c in candidates:
        if not c['inchikey']:
            continue
        accession = c['source_id'] if c['source'] == 'massbank' else None
        for record in catalog.get((c['inchikey'], accession), []):
            rows.append(dict(record, evidence_id=stable_id('evidence', [c['candidate_id'], canonical_json(record)]),
                             candidate_id=c['candidate_id'], feature_id=c['feature_id'],
                             source='massbank_rdf_index',
                             status='resolved' if record.get('accession') or record.get('compound')
                             or record.get('classes') else 'empty',
                             provenance=provenance))
    rows.sort(key=lambda e: e['evidence_id'])
    return rows


def tabulate(page: dict) -> dict:
    """行の一覧を列名1回のTSVにする（戻り値の規約）。"""
    return {k: rows_tsv(v) if k in TABLE_KEYS else v for k, v in page.items()}


def _encode(value: dict) -> str:
    return base64.urlsafe_b64encode(canonical_json(value).encode()).decode()


def _decode_cursor(cursor: str, expected: dict, n_ids: int, limit: int) -> dict:
    try:
        cur = json.loads(base64.urlsafe_b64decode(cursor.encode()))
    except (ValueError, TypeError, UnicodeError) as exc:
        raise ValueError('invalid_cursor') from exc
    if not isinstance(cur, dict):
        raise ValueError('invalid_cursor')
    if any(cur.get(k) != v for k, v in expected.items()):
        raise ValueError('stale_cursor: 解析・条件または根拠版が変わりました')
    offset = cur.get('offset')
    page_size = cur.get('page_size', limit)
    evidence_offset = cur.get('evidence_offset', 0)
    if (type(offset) is not int or not 0 <= offset <= n_ids
            or type(page_size) is not int or not 1 <= page_size <= limit
            or type(evidence_offset) is not int or evidence_offset < 0):
        raise ValueError('invalid_cursor')
    return dict(offset=offset, page_size=page_size, evidence_offset=evidence_offset,
                detail_page=cur.get('mode') == 'evidence')


def population_projection(bundle, candidates, matches, excluded=frozenset()):
    return build_population(bundle, candidates, matches.get('features', {}), excluded)


def get_evidence(store, analysis_id: str, *, feature_ids=None, cursor=None, limit=50):
    return evidence_page(store.read_view(analysis_id), analysis_id,
                         feature_ids=feature_ids, cursor=cursor, limit=limit)


def evidence_page(view, analysis_id: str, *, feature_ids=None, cursor=None, limit=50):
    """read_viewで固定した版から1ページを組み立てる。"""
    if type(limit) is not int or not 1 <= limit <= 200:
        raise ValueError('limitは1〜200です')
    bundle, revision, state = view
    matches = state.get('matches', {})
    lookup_state = state.get('lookups', {})
    context = state.get('context', {})
    candidates = candidate_dicts(bundle, matches)
    ids = select_ids([f.feature_id for f in bundle.features], feature_ids)
    base = dict(analysis_id=analysis_id, revision=revision,
                filter=stable_id('filter', [canonical_json(ids)]), limit=limit)
    cur = dict(offset=0, page_size=limit, evidence_offset=0, detail_page=False)
    if cursor:
        cur = _decode_cursor(cursor, base, len(ids), limit)

    # ページに依存しない全件ベースの部分は1回だけ計算する。
    all_candidate_ids = {c['candidate_id'] for c in candidates}
    counts = Counter(l['status'] for l in lookup_state.get('lookups', [])
                     if l['candidate_id'] in all_candidate_ids)
    linked = {l.feature_id for l in bundle.spectrum_links}
    annotated = {c['feature_id'] for c in candidates if c['status'] == 'active'}
    population = bundle.manifest['population']
    summary = dict(n_features_total=population['n_features_total'],
                   n_features_available=len(bundle.features), population_complete=bundle.population_complete,
                   n_significant=sum(f.significant for f in bundle.features),
                   n_up=sum(f.log2fc is not None and f.log2fc > 0 for f in bundle.features),
                   n_down=sum(f.log2fc is not None and f.log2fc < 0 for f in bundle.features),
                   n_missing_log2fc=sum(f.log2fc is None for f in bundle.features),
                   n_with_candidates=len(annotated), n_with_msms=len(linked),
                   n_candidates=len(candidates), lookup_statuses=dict(counts),
                   n_not_queried=len(candidates)-sum(counts.values())+counts.get('not_queried', 0))
    projection = population_projection(bundle, candidates, matches)
    eligibility = dict(projection.exclusions)
    for ik, representative in projection.representatives.items():
        for fid in representative['feature_ids']:
            eligibility[fid] = dict(reason='eligible', inchikey=ik,
                                    representative_feature_id=representative['feature_id'])
    lookups_by_id = {l['candidate_id']: l for l in lookup_state.get('lookups', [])}
    common = dict(analysis_id=analysis_id, evidence_revision=revision,
                  context_revision=context.get('context_revision'), context=context,
                  comparison=bundle.manifest['comparison'], summary=summary,
                  population_derivation=dict(
                      n_features_before=population.get('n_features_before'),
                      exclusions=population.get('exclusions', []),
                      caveats=(population.get('preprocessing_report') or {}).get('caveats', []),
                      derivation_status=population.get('derivation_status', 'not_available')),
                  matching_provenance=dict(index_fingerprint=matches.get('index_fingerprint'),
                                           index_meta=matches.get('index_meta')),
                  scope=lookup_state.get('scope', ['massbank_metadata', 'pubchem_pathways']),
                  n_selected=len(ids))

    offset, page_size = cur['offset'], cur['page_size']
    evidence_offset, detail_page = cur['evidence_offset'], cur['detail_page']
    while True:
        shown = set(ids[offset:offset+page_size])
        selected = [c for c in candidates if c['feature_id'] in shown]
        lookups = []
        for cid in sorted(c['candidate_id'] for c in selected):
            raw = lookups_by_id.get(cid, dict(candidate_id=cid, status='not_queried'))
            lookups.append(dict({k: v for k, v in raw.items() if k != 'hits'}, n_hits=len(raw.get('hits', []))))
        evidence = evidence_rows(selected, lookup_state)
        if evidence_offset > len(evidence):
            raise ValueError('invalid_cursor')
        # 詳細ページでも、残りの特徴量ページへのcursorを失わない。
        next_cursor = (_encode(dict(base, offset=offset+page_size, page_size=limit))
                       if offset+page_size < len(ids) else None)
        next_evidence = (_encode(dict(base, offset=offset, page_size=page_size, mode='evidence',
                                      evidence_offset=evidence_offset+EVIDENCE_PAGE_SIZE))
                         if evidence_offset+EVIDENCE_PAGE_SIZE < len(evidence) else None)
        result = dict(common,
                      features=[dict(asdict(f), statistical_eligibility=eligibility.get(f.feature_id))
                                for f in bundle.features if f.feature_id in shown],
                      candidates=selected,
                      evidence=evidence[evidence_offset:evidence_offset+EVIDENCE_PAGE_SIZE],
                      lookups=lookups,
                      matching={fid: matches.get('features', {}).get(
                                    fid, dict(status='not_matched' if fid in linked else 'no_msms'))
                                for fid in sorted(shown)},
                      truncated=next_cursor is not None, next_cursor=next_cursor,
                      detail_page=detail_page, n_evidence=len(evidence),
                      evidence_truncated=next_evidence is not None, next_evidence_cursor=next_evidence)
        # 件数上限より先にバイト上限へ達した場合は件数を減らす。cursorで残りへ到達できる。
        if (detail_page or page_size == 1
                or len(canonical_json(tabulate(result)).encode()) <= MAX_PAGE_BYTES):
            return result
        page_size = max(1, page_size // 2)
