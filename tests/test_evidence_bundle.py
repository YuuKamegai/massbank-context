"""候補の有無に依存しない入力契約。"""
import json
import pytest
from evidence_helpers import make_bundle, rewrite_table


def test_mixed_bundle_keeps_all_features(tmp_path):
    from massbank_context.experiment.bundle import load_bundle
    bundle = load_bundle(make_bundle(tmp_path / '日本語'))
    assert len(bundle.features) == 8
    assert bundle.population_complete is True
    assert any(f.log2fc is None for f in bundle.features)
    assert '\t' in bundle.annotations[0].name
    assert load_bundle(bundle.manifest_path).analysis_id == bundle.analysis_id


@pytest.mark.parametrize('column,value', [('p_value', '1.1'), ('q_value', '-.1'),
                                         ('log2fc', 'NaN'), ('significant', 'yes')])
def test_invalid_values_rejected(tmp_path, column, value):
    from massbank_context.experiment.bundle import load_bundle, BundleError
    path = make_bundle(tmp_path)
    rewrite_table(path, 'features', lambda rows: rows[0].update({column: value}))
    with pytest.raises(BundleError):
        load_bundle(path)


def test_duplicate_and_unlinked_features_rejected(tmp_path):
    from massbank_context.experiment.bundle import load_bundle, BundleError
    path = make_bundle(tmp_path)
    rewrite_table(path, 'annotations', lambda rows: rows[0].update(feature_id='other-run'))
    with pytest.raises(BundleError):
        load_bundle(path)


def test_hash_and_path_escape_rejected(tmp_path):
    from massbank_context.experiment.bundle import load_bundle, BundleError
    path = make_bundle(tmp_path)
    data = json.loads(path.read_text())
    data['files'][0]['path'] = '../features.tsv'
    path.write_text(json.dumps(data))
    with pytest.raises(BundleError):
        load_bundle(path)


def test_content_change_and_direction_rejected(tmp_path):
    from massbank_context.experiment.bundle import load_bundle, BundleError
    path = make_bundle(tmp_path)
    (tmp_path / 'features.tsv').write_text('changed')
    with pytest.raises(BundleError, match='ハッシュ'):
        load_bundle(path)
    path = make_bundle(tmp_path)
    data = json.loads(path.read_text())
    data['comparison']['log2fc_sign'] = 'A higher'
    path.write_text(json.dumps(data))
    with pytest.raises(BundleError, match='方向'):
        load_bundle(path)
