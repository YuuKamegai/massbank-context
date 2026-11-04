"""公開層のエラーと永続化。"""
import json
from evidence_helpers import make_bundle, read_tsv


def test_public_tools_keep_population(tmp_path, monkeypatch):
    from massbank_context.core import config
    from massbank_context.tools import evidence_tools as tools
    monkeypatch.setattr(config, 'ARTIFACTS_PATH', tmp_path / 'store')
    loaded = json.loads(tools.load_evidence_bundle(str(make_bundle(tmp_path / 'input'))))
    aid = loaded['analysis_id']
    page = json.loads(tools.get_comparison_evidence(aid))
    assert page['summary']['n_features_total'] == 8
    assert len(read_tsv(page['features'])) == 8
    error = json.loads(tools.get_comparison_evidence('missing'))
    assert error['error']['code'] == 'invalid_evidence_request'


def test_small_p_values_are_not_rounded_to_zero(tmp_path, monkeypatch):
    from evidence_helpers import rewrite_table
    from massbank_context.core import config
    from massbank_context.tools import evidence_tools as tools
    monkeypatch.setattr(config,'ARTIFACTS_PATH',tmp_path/'store')
    path=make_bundle(tmp_path/'input')
    rewrite_table(path,'features',lambda rows:rows[0].update(p_value='2.1e-13',q_value='3.2e-10'))
    aid=json.loads(tools.load_evidence_bundle(str(path)))['analysis_id']
    feature=read_tsv(json.loads(tools.get_comparison_evidence(aid))['features'])[0]
    assert float(feature['p_value']) == 2.1e-13
    assert float(feature['q_value']) == 3.2e-10
