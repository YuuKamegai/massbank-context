"""RDFの取得不能と関連なしは別の状態。"""
from evidence_helpers import KEY_A, make_bundle


def test_endpoint_cache_separation(tmp_path):
    from massbank_context.kg.evidence_cache import EvidenceCache
    cache = EvidenceCache(tmp_path / 'cache.sqlite')
    cache.put(dict(endpoint='https://one.example', query_version='1', inchikey=KEY_A,
                   status='resolved', hits=[{'pathway_id': 'PW1'}]))
    assert cache.lookup('https://two.example', '1', KEY_A) is None


def test_refresh_failure_keeps_old_version_but_not_current_success(tmp_path):
    from massbank_context.evidence.store import EvidenceStore
    from massbank_context.experiment.bundle import load_bundle
    from massbank_context.kg.evidence import resolve_candidates
    from massbank_context.core import config
    store = EvidenceStore(tmp_path / 'store')
    aid = store.register(load_bundle(make_bundle(tmp_path / 'input')))

    class Adapter:
        endpoint = 'https://fixture.example'
        failed = []
        def fetch(self, keys):
            return {key: [{'pathway_id':'PW1', 'pathway_label':'Path', 'taxid':'10090'}] for key in keys}

    good = resolve_candidates(store, aid, None, adapter=Adapter())
    assert good['n_resolved'] > 0

    class Failed(Adapter):
        def fetch(self, keys):
            self.failed = keys
            return {}

    bad = resolve_candidates(store, aid, None, refresh=True, adapter=Failed())
    assert bad['n_failed'] > 0
    now = store.current(aid, 'lookups')
    assert all(x['status'] != 'resolved' for x in now['lookups'])


def test_invalid_sparql_payload_is_failure_not_empty():
    from massbank_context.kg.evidence import EvidenceAdapter, fetch_pathway_evidence
    invalid_optional = {'results':{'bindings':[dict(value_inchikey={'value':KEY_A},
        pathway={'value':'http://example.org/PW1'},pathway_organism={'value':None})]}}
    for payload in ({'error':'query timeout'}, {}, {'results':{'bindings':[{}]}}, invalid_optional):
        adapter = EvidenceAdapter('https://fixture.example', transport=lambda *args:payload)
        result = fetch_pathway_evidence(adapter,[KEY_A])
        assert result[KEY_A]['status'] == 'failed'
