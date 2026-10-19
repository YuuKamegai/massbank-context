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
