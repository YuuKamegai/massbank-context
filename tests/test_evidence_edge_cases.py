"""取得範囲と識別子欠落で誤った結論を返さない。"""
import json
import pytest
from evidence_helpers import make_bundle, rewrite_table
from test_comparison_evidence import setup_store


def test_mapping_failure_is_not_no_msms(tmp_path):
    from massbank_context.experiment.bundle import load_bundle
    from massbank_context.evidence.store import EvidenceStore
    from massbank_context.evidence.matching import match_features
    path = make_bundle(tmp_path/'input')
    data = json.loads(path.read_text())
    data['provenance'] = dict(msms_link_statuses={'f1': {'status':'precursor_mismatch'}})
    path.write_text(json.dumps(data))
    store = EvidenceStore(tmp_path/'store')
    aid = store.register(load_bundle(path))
    assert match_features(store,aid,None)['features']['f1']['status'] == 'msms_mapping_failed'


def test_cursor_shape_error_is_machine_readable(tmp_path):
    from massbank_context.report.comparison import get_evidence
    store, aid = setup_store(tmp_path)
    with pytest.raises(ValueError, match='invalid_cursor'):
        get_evidence(store,aid,cursor='W10=')
