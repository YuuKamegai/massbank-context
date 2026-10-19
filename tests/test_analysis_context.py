"""確認済み背景を暗黙に書き換えない。"""
from evidence_helpers import make_bundle


def test_context_preserves_old_confirmation(tmp_path):
    from massbank_context.evidence.store import EvidenceStore
    from massbank_context.evidence.context import set_context
    from massbank_context.experiment.bundle import load_bundle
    store = EvidenceStore(tmp_path / 'store')
    aid = store.register(load_bundle(make_bundle(tmp_path / 'input')))
    old = set_context(store, aid, {'species': {'value': 'mouse', 'source': '研究者'}}, True)
    new = set_context(store, aid, {'species': {'value': None}}, False)
    assert old['status'] == 'confirmed'
    assert new['status'] == 'proposed'
    assert old['context_revision'] != new['context_revision']
