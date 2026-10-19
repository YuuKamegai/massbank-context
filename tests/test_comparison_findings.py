"""LLMの主張を存在する根拠へ結び付ける。"""
from test_comparison_evidence import setup_store


def test_unknown_evidence_is_invalid(tmp_path):
    from massbank_context.report.comparison import get_evidence
    from massbank_context.report.comparison_findings import validate_claim_references
    store, aid = setup_store(tmp_path)
    result = validate_claim_references([{'claim_id':'x', 'evidence_ids':['missing']}], get_evidence(store, aid))
    assert result[0]['code'] == 'unknown_evidence'


def test_unconfirmed_context_still_returns_facts(tmp_path):
    from massbank_context.report.comparison_findings import build_comparison_findings
    store, aid = setup_store(tmp_path)
    result = build_comparison_findings(store, aid, organism='10090', include_pathways=True)
    assert result['facts']['summary']['n_features_total'] == 8
    assert result['readiness']['context_confirmed'] is False
    assert result['pathways']['reason'] == 'organism_unconfirmed'
    assert result['hypotheses'][0]['claim'] is None
