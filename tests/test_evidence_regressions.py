"""レビュー指摘の回帰検証（空スペクトル・識別子・cursor・引用検証・索引の版）。"""
import hashlib
import json
from evidence_helpers import KEY_A, make_bundle, rewrite_table
from test_comparison_evidence import setup_store
from test_matcher import conn


def _bundle_with_spectra(tmp_path, links):
    """f1にMS/MSを対応させたバンドル。linksは(spectrum_id)の並び。"""
    path = make_bundle(tmp_path / 'input')
    msp = path.parent / 'spectra.msp'
    msp.write_text('Name: empty\nPrecursorMZ: 230.0552\nNum Peaks: 0\n\n'
                   'Name: synthetic\nPrecursorMZ: 230.0552\nNum Peaks: 2\n86.0349 20\n171.0069 999\n')
    manifest = json.loads(path.read_text())
    manifest['files'].append(dict(role='spectra', path=msp.name,
                                  sha256=hashlib.sha256(msp.read_bytes()).hexdigest()))
    path.write_text(json.dumps(manifest))
    rewrite_table(path, 'spectrum_links', lambda rows: rows.extend(
        dict(feature_id='f1', spectrum_file=msp.name, spectrum_id=sid,
             source_analysis_id='synthetic-run', mapping_evidence='fixture-explicit-id') for sid in links))
    return path


def _store(tmp_path, path):
    from massbank_context.evidence.store import EvidenceStore
    from massbank_context.experiment.bundle import load_bundle
    store = EvidenceStore(tmp_path / 'store')
    return store, store.register(load_bundle(path))


def _use_index(monkeypatch, conn):
    from pathlib import Path
    from massbank_context.core import config
    monkeypatch.setattr(config, 'INDEX_PATH', Path(conn.execute('PRAGMA database_list').fetchone()[2]))


def test_one_empty_spectrum_does_not_block_matched_feature(tmp_path, monkeypatch, conn):
    from massbank_context.evidence.matching import match_features
    from massbank_context.evidence.selection import candidate_dicts
    from massbank_context.report.comparison_population import build_population
    _use_index(monkeypatch, conn)
    store, aid = _store(tmp_path, _bundle_with_spectra(tmp_path, ['spectra:1', 'spectra:2']))
    status = match_features(store, aid, None)['features']['f1']
    assert status['status'] == 'matched'
    assert status['n_empty_spectra'] == 1
    bundle, _, state = store.read_view(aid)
    projection = build_population(bundle, candidate_dicts(bundle, state['matches']),
                                  state['matches']['features'])
    assert projection.summary['n_pending_matching'] == 0


def test_all_empty_spectra_is_terminal_not_pending(tmp_path, monkeypatch, conn):
    from massbank_context.evidence.matching import match_features
    from massbank_context.evidence.selection import candidate_dicts
    from massbank_context.report.comparison_population import build_population
    _use_index(monkeypatch, conn)
    store, aid = _store(tmp_path, _bundle_with_spectra(tmp_path, ['spectra:1']))
    assert match_features(store, aid, None)['features']['f1']['status'] == 'empty_spectrum'
    bundle, _, state = store.read_view(aid)
    projection = build_population(bundle, candidate_dicts(bundle, state['matches']),
                                  state['matches']['features'])
    assert projection.summary['n_pending_matching'] == 0
    assert projection.exclusions.get('f1', {}).get('reason') != 'matching_incomplete'


def test_candidates_carry_index_fingerprint_not_full_meta(tmp_path, monkeypatch, conn):
    from massbank_context.evidence.matching import match_features
    _use_index(monkeypatch, conn)
    store, aid = _store(tmp_path, _bundle_with_spectra(tmp_path, ['spectra:2']))
    match_features(store, aid, None)
    matches = store.current(aid, 'matches')
    assert matches['index_meta']['bin_width']
    assert all(c['provenance']['index_fingerprint'] == matches['index_fingerprint']
               and 'index_meta' not in c['provenance'] for c in matches['candidates'])


def test_index_fingerprint_follows_meta(tmp_path):
    from massbank_context.index.writer import index_fingerprint
    path = tmp_path / 'index.sqlite'
    path.write_bytes(b'x' * 10)
    meta = dict(built_at='2026-10-19T00:00:00+00:00')
    assert index_fingerprint(path, meta) == index_fingerprint(path, dict(meta))
    assert index_fingerprint(path, meta) != index_fingerprint(path, dict(built_at='2026-10-20T00:00:00+00:00'))


def test_malformed_inchikey_is_unsupported_not_remaining(tmp_path):
    from massbank_context.kg.evidence import resolve_candidates
    store, aid = setup_store(tmp_path)
    store.put_snapshot(aid, 'matches', dict(features={}, candidates=[dict(
        candidate_id='candidate:legacy', feature_id='f4', source='legacy_annotation',
        source_id='x', spectrum_id=None, name='', inchikey='ABCDEF', structure_ids={},
        rank=None, scores={}, status='active', same_structure_as=[], provenance={}, curation={})]))

    class Adapter:
        endpoint = 'https://fixture.example'
        failed = []
        def fetch(self, keys):
            assert 'ABCDEF' not in keys
            return {}

    result = resolve_candidates(store, aid, ['candidate:legacy'], adapter=Adapter())
    assert result['n_remaining'] == 0
    lookup = store.current(aid, 'lookups')['lookups'][0]
    assert lookup['status'] == 'unsupported_identifier'


def test_evidence_detail_page_keeps_feature_cursor(tmp_path):
    from massbank_context.report.comparison import get_evidence
    store, aid = setup_store(tmp_path)
    candidate = store.load(aid).annotations[0]
    rows = [dict(evidence_id=f'e{i:03}', candidate_id=candidate.candidate_id,
                 feature_id=candidate.feature_id, pathway_id=f'PW{i}') for i in range(130)]
    store.put_snapshot(aid, 'lookups', dict(evidence=rows, lookups=[], catalog=[]))
    first = get_evidence(store, aid, limit=1)
    assert first['next_cursor'] and first['next_evidence_cursor']
    detail = get_evidence(store, aid, limit=1, cursor=first['next_evidence_cursor'])
    assert detail['detail_page']
    assert detail['truncated'] is True
    assert detail['next_cursor'] == first['next_cursor']


def test_claims_are_validated_against_all_pages(tmp_path):
    from massbank_context.report.comparison import get_evidence
    from massbank_context.report.comparison_findings import validate_claims
    store, aid = setup_store(tmp_path)
    last = store.load(aid).annotations[-2]
    assert last.feature_id not in {f['feature_id'] for f in get_evidence(store, aid, limit=1)['features']}
    row = dict(evidence_id='e-last', candidate_id=last.candidate_id, feature_id=last.feature_id,
               pathway_id='PW1')
    store.put_snapshot(aid, 'lookups', dict(evidence=[row], lookups=[], catalog=[]))
    ok = validate_claims(store, aid, [dict(claim_id='c', feature_ids=[last.feature_id],
                                           candidate_ids=[last.candidate_id], evidence_ids=['e-last'])])
    assert ok['valid'], ok['issues']
    stale = validate_claims(store, aid, [dict(claim_id='c', evidence_revision='revision:old')])
    assert [i['code'] for i in stale['issues']] == ['evidence_revision_mismatch']


def test_snapshot_write_does_not_decode_bundle(tmp_path, monkeypatch):
    from massbank_context.evidence.store import EvidenceStore
    store, aid = setup_store(tmp_path)

    def fail(*args):
        raise AssertionError('入力全体を復元した')

    monkeypatch.setattr(EvidenceStore, 'load', fail)
    store.put_snapshot(aid, 'matches', {'n': 1}, expected_revision=store.latest_revision(aid))
    assert store.current(aid, 'matches') == {'n': 1}


def test_public_page_lists_are_tsv(tmp_path, monkeypatch):
    from evidence_helpers import read_tsv
    from massbank_context.core import config
    from massbank_context.tools import evidence_tools as tools
    monkeypatch.setattr(config, 'ARTIFACTS_PATH', tmp_path / 'store')
    aid = json.loads(tools.load_evidence_bundle(str(make_bundle(tmp_path / 'input'))))['analysis_id']
    page = json.loads(tools.get_comparison_evidence(aid))
    candidates = read_tsv(page['candidates'])
    # 値中のタブがTSVの列をずらさない（fixtureの名前は「名前<TAB>候補」）。
    assert {c['name'] for c in candidates} == {'名前\t候補'}
    assert all(c['inchikey'] in ('', KEY_A) or len(c['inchikey']) == 27 for c in candidates)
