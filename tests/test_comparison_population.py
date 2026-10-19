"""曖昧な候補を有意性で選ばない。"""
from dataclasses import asdict
from evidence_helpers import make_bundle, KEY_CONFLICT


def test_competing_and_opposite_spots_excluded(tmp_path):
    from massbank_context.experiment.bundle import load_bundle
    from massbank_context.report.comparison_population import build_population
    bundle = load_bundle(make_bundle(tmp_path))
    result = build_population(bundle, list(bundle.annotations), {})
    assert result.exclusions['f0']['reason'] == 'candidate_conflict'
    assert result.exclusions['f5']['reason'] == 'direction_conflict'
    assert result.exclusions['f6']['reason'] == 'direction_conflict'
    assert all(c['inchikey'] != KEY_CONFLICT for c in result.compounds)
    assert not result.representatives[next(iter(result.representatives))]['significant'] is None


def test_truncation_and_failed_lookup_suppress_statistics(tmp_path):
    from massbank_context.experiment.bundle import load_bundle
    from massbank_context.report.comparison_population import build_population, compute_pathways
    bundle = load_bundle(make_bundle(tmp_path))
    result = build_population(bundle, list(bundle.annotations), {'f1': {'truncated': True, 'status': 'matched'}})
    assert result.exclusions['f1']['reason'] == 'candidate_truncated'
    missing = compute_pathways(result, [], '10090')
    assert missing['tested'] is False
    assert missing['reason'] == 'lookup_incomplete'
