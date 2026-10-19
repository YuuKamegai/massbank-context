"""入力と根拠版を再起動後にも復元できる。"""
from evidence_helpers import make_bundle


def test_store_reopens_same_analysis(tmp_path):
    from massbank_context.experiment.bundle import load_bundle
    from massbank_context.evidence.store import EvidenceStore
    bundle = load_bundle(make_bundle(tmp_path / 'input'))
    root = tmp_path / 'store'
    store = EvidenceStore(root)
    aid = store.register(bundle)
    rev = store.put_snapshot(aid, 'matches', {'candidates': [], 'status': 'done'})
    assert EvidenceStore(root).load(aid).analysis_id == aid
    assert EvidenceStore(root).get_snapshot(aid, rev)['status'] == 'done'
    assert store.register(bundle) == aid


def test_snapshot_survives_new_revision(tmp_path):
    from massbank_context.experiment.bundle import load_bundle
    from massbank_context.evidence.store import EvidenceStore
    store = EvidenceStore(tmp_path / 'store')
    aid = store.register(load_bundle(make_bundle(tmp_path / 'input')))
    old = store.put_snapshot(aid, 'matches', {'number': 1})
    new = store.put_snapshot(aid, 'matches', {'number': 2})
    assert old != new
    assert store.get_snapshot(aid, old)['number'] == 1
    assert store.current(aid, 'matches')['number'] == 2


def test_outdated_writer_does_not_overwrite_new_state(tmp_path):
    import pytest
    from massbank_context.experiment.bundle import load_bundle
    from massbank_context.evidence.store import EvidenceStore
    store = EvidenceStore(tmp_path / 'store')
    aid = store.register(load_bundle(make_bundle(tmp_path / 'input')))
    revision = store.latest_revision(aid)
    store.put_snapshot(aid, 'matches', {'number': 2})
    with pytest.raises(ValueError, match='stale_revision'):
        store.put_snapshot(aid, 'matches', {'number': 1}, expected_revision=revision)
    assert store.current(aid,'matches')['number'] == 2
