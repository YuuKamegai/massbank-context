"""LLMの主張を存在する根拠へ結び付ける。"""
import re
import pytest
from evidence_helpers import KEY_A
from test_comparison_evidence import setup_store

# 実際の PubChem RDF の形。生物種 URI の末尾は "TAXID10090" で、数字だけではない。
_PUBCHEM = 'http://rdf.ncbi.nlm.nih.gov/pubchem'


def _pubchem_transport(endpoint, query):
    """KEY_A にマウスの PWID1 とヒトの PWID2 を返す。ほかの InChIKey は空。"""
    bindings = []
    if KEY_A in re.findall(r'"([A-Z]{14}-[A-Z]{10}-[A-Z])"', query):
        for pathway, taxid in (('PWID1', 'TAXID10090'), ('PWID2', 'TAXID9606')):
            bindings.append(dict(value_inchikey={'value': KEY_A},
                                 pathway={'value': f'{_PUBCHEM}/pathway/{pathway}'},
                                 pathway_label={'value': f'Pathway {pathway}'},
                                 pathway_organism={'value': f'{_PUBCHEM}/taxonomy/{taxid}'}))
    return {'results': {'bindings': bindings}}


def _confirmed_store(tmp_path, context_taxid):
    from massbank_context.evidence.context import set_context
    from massbank_context.kg.evidence import EvidenceAdapter, resolve_candidates
    store, aid = setup_store(tmp_path)
    set_context(store, aid, {'species': {'value': {'label': 'mouse', 'taxid': context_taxid},
                                         'source': 'fixture researcher'}}, True)
    resolve_candidates(store, aid, None,
                       adapter=EvidenceAdapter('https://fixture.example', transport=_pubchem_transport))
    return store, aid


@pytest.mark.parametrize('context_taxid, organism', [
    ('10090', '10090'),            # 手順書（two-group-interpretation.md）どおりの書式
    ('10090', 'TAXID10090'),       # 根拠のページに出る書式を organism に写した
    ('TAXID10090', '10090'),
    ('TAXID10090', 'TAXID10090'),
])
def test_pathways_match_pubchem_taxid_in_either_form(tmp_path, context_taxid, organism):
    from massbank_context.report.comparison_findings import build_comparison_findings
    store, aid = _confirmed_store(tmp_path, context_taxid)
    pathways = build_comparison_findings(store, aid, organism=organism, include_pathways=True)['pathways']
    # 母集団に残るのは KEY_A（f1・f2）だけ。マウスの PWID1 だけが数えられ、ヒトの PWID2 は入らない。
    assert pathways['reason'] is None, pathways['reason']
    assert pathways['n_with_species_pathways'] == 1
    assert [row['pathway_id'] for row in pathways['rows']] == ['PWID1']
    assert pathways['rows'][0]['feature_ids'] == ['f1', 'f2']


def test_other_species_is_still_a_context_mismatch(tmp_path):
    from massbank_context.report.comparison_findings import build_comparison_findings
    store, aid = _confirmed_store(tmp_path, '10090')
    pathways = build_comparison_findings(store, aid, organism='TAXID9606', include_pathways=True)['pathways']
    assert pathways['tested'] is False
    assert pathways['reason'] == 'organism_context_mismatch'


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
