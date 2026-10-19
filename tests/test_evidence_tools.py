"""公開層のエラーと永続化。"""
import json
from evidence_helpers import make_bundle


def test_public_tools_keep_population(tmp_path, monkeypatch):
    from massbank_context.core import config
    from massbank_context.tools import evidence_tools as tools
    monkeypatch.setattr(config, 'ARTIFACTS_PATH', tmp_path / 'store')
    loaded = json.loads(tools.load_evidence_bundle(str(make_bundle(tmp_path / 'input'))))
    aid = loaded['analysis_id']
    page = json.loads(tools.get_comparison_evidence(aid))
    assert page['summary']['n_features_total'] == 8
    assert len(page['features']) == 8
    error = json.loads(tools.get_comparison_evidence('missing'))
    assert error['error']['code'] == 'invalid_evidence_request'
