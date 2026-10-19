"""ページ選択で全体の母集団や版を変えない。"""
import pytest
from evidence_helpers import make_bundle


def setup_store(tmp_path):
    from massbank_context.evidence.store import EvidenceStore
    from massbank_context.experiment.bundle import load_bundle
    store = EvidenceStore(tmp_path / 'store')
    aid = store.register(load_bundle(make_bundle(tmp_path / 'input')))
    return store, aid


def test_page_does_not_change_population(tmp_path):
    from massbank_context.report.comparison import get_evidence
    store, aid = setup_store(tmp_path)
    page = get_evidence(store, aid, limit=1)
    assert page['summary']['n_features_total'] == 8
    assert len(page['features']) == 1
    assert page['next_cursor']
    second = get_evidence(store, aid, limit=1, cursor=page['next_cursor'])
    assert page['features'][0]['feature_id'] != second['features'][0]['feature_id']
    assert get_evidence(store, aid, feature_ids=[])['features'] == []


def test_context_change_invalidates_cursor(tmp_path):
    from massbank_context.report.comparison import get_evidence
    from massbank_context.evidence.context import set_context
    store, aid = setup_store(tmp_path)
    page = get_evidence(store, aid, limit=1)
    set_context(store, aid, {}, False)
    with pytest.raises(ValueError, match='stale_cursor'):
        get_evidence(store, aid, limit=1, cursor=page['next_cursor'])


def test_large_rdf_details_are_separately_pageable(tmp_path):
    from massbank_context.report.comparison import get_evidence
    store, aid = setup_store(tmp_path)
    candidate = store.load(aid).annotations[0]
    rows = [dict(evidence_id=f'e{i:03}',candidate_id=candidate.candidate_id,
                 feature_id=candidate.feature_id,pathway_id=f'PW{i}') for i in range(130)]
    store.put_snapshot(aid,'lookups',dict(evidence=rows,lookups=[dict(candidate_id=candidate.candidate_id,
        status='resolved',hits=rows)],catalog=[]))
    page = get_evidence(store,aid)
    assert len(page['evidence']) <= 64
    assert 'hits' not in page['lookups'][0]
    ids = {e['evidence_id'] for e in page['evidence']}
    while page['next_evidence_cursor']:
        page = get_evidence(store,aid,cursor=page['next_evidence_cursor'])
        ids.update(e['evidence_id'] for e in page['evidence'])
    assert len(ids) == 130


def test_explicit_unqueried_lookup_counts_as_unqueried(tmp_path):
    from massbank_context.report.comparison import get_evidence
    store, aid = setup_store(tmp_path)
    candidate = store.load(aid).annotations[0]
    store.put_snapshot(aid,'lookups',dict(evidence=[],lookups=[
        dict(candidate_id=candidate.candidate_id,status='not_queried')]))
    summary = get_evidence(store,aid)['summary']
    assert summary['n_not_queried'] == summary['n_candidates']


def test_large_feature_page_keeps_reachable_next_page(tmp_path):
    import json
    from evidence_helpers import rewrite_table
    from massbank_context.evidence.store import EvidenceStore
    from massbank_context.experiment.bundle import load_bundle
    from massbank_context.report.comparison import get_evidence
    path=make_bundle(tmp_path/'input')
    rewrite_table(path,'annotations',lambda rows:[r.update(name='x'*18000) for r in rows])
    store=EvidenceStore(tmp_path/'store')
    aid=store.register(load_bundle(path))
    page=get_evidence(store,aid,limit=8)
    assert len(page['features']) < 8
    assert len(json.dumps(page,ensure_ascii=False).encode()) < 128*1024
    ids={f['feature_id'] for f in page['features']}
    while page['next_cursor']:
        page=get_evidence(store,aid,limit=8,cursor=page['next_cursor'])
        ids.update(f['feature_id'] for f in page['features'])
    assert len(ids)==8


def test_small_page_does_not_permanently_shrink_next_pages(tmp_path):
    from evidence_helpers import rewrite_table
    from massbank_context.evidence.store import EvidenceStore
    from massbank_context.experiment.bundle import load_bundle
    from massbank_context.report.comparison import get_evidence
    path=make_bundle(tmp_path/'input')
    rewrite_table(path,'annotations',lambda rows:[r.update(name='x'*120000) for r in rows if r['feature_id']=='f0'])
    store=EvidenceStore(tmp_path/'store');aid=store.register(load_bundle(path))
    first=get_evidence(store,aid,limit=8)
    assert len(first['features']) == 1
    second=get_evidence(store,aid,limit=8,cursor=first['next_cursor'])
    assert len(second['features']) > 1
