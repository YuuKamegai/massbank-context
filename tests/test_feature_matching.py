"""順位と切り詰め情報を区別する。"""
from test_matcher import conn
from massbank_context.core.models import QuerySpectrum, Peak


def test_match_result_preserves_existing_rank(conn):
    from massbank_context.spectral.matcher import match_spectrum, match_spectrum_result
    query = QuerySpectrum('q', None, None, (Peak(100, 999),))
    result = match_spectrum_result(conn, query, 'LC-ESI-MS2', top_k=1)
    old = match_spectrum(conn, query, 'LC-ESI-MS2', top_k=1)
    assert result.candidates == old
    assert result.n_total >= len(result.candidates)
    assert result.truncated == (result.n_total > 1)


def test_no_msms_does_not_require_index(tmp_path):
    from evidence_helpers import make_bundle
    from massbank_context.experiment.bundle import load_bundle
    from massbank_context.evidence.store import EvidenceStore
    from massbank_context.evidence.matching import match_features
    store = EvidenceStore(tmp_path / 'store')
    aid = store.register(load_bundle(make_bundle(tmp_path / 'input')))
    result = match_features(store, aid, None)
    assert all(x['status'] == 'no_msms' for x in result['features'].values())


def test_catalog_keeps_license_and_unknown_ttl(conn):
    from massbank_context.evidence.catalog import read_candidate_catalog
    accession = conn.execute('SELECT accession FROM records LIMIT 1').fetchone()[0]
    result = read_candidate_catalog(conn, accessions=[accession], inchikeys=[])
    assert result[0]['accession'] == accession
    assert 'license' in result[0]
    assert result[0]['ttl_source'] is None
