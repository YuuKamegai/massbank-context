import json

import pytest

from massbank_context.evidence.store import EvidenceStore
from massbank_context.experiment.bundle import (
    INTENSITY_COLUMNS, BundleError, load_bundle, parse_intensities, parse_samples, role_path)
from tests.evidence_helpers import change_manifest, rewrite_table
from tests.xenobiotic_helpers import make_v3_bundle


def _edit_comparison(path, change):
    data = json.loads(path.read_text(encoding='utf-8'))
    change(data['comparison'])
    path.write_text(json.dumps(data), encoding='utf-8')


def break_bundle(path, breakage):
    if breakage == 'samples_without_intensities':
        data = json.loads(path.read_text(encoding='utf-8'))
        data['files'] = [f for f in data['files'] if f['role'] != 'intensities']
        path.write_text(json.dumps(data), encoding='utf-8')
    elif breakage == 'v2_with_samples':
        change_manifest(path, contract_version=2)
    elif breakage == 'bad_group':
        rewrite_table(path, 'samples', lambda rows: rows[0].update(group='x'))
    elif breakage == 'n_b_mismatch':
        _edit_comparison(path, lambda c: c.update(n_b=5))
    elif breakage == 'missing_intensity_row':
        rewrite_table(path, 'intensities', lambda rows: rows.pop())
    elif breakage == 'duplicate_intensity_row':
        rewrite_table(path, 'intensities', lambda rows: rows.append(dict(rows[0])))
    elif breakage == 'unknown_sample':
        rewrite_table(path, 'intensities', lambda rows: rows[0].update(sample_id='ZZ'))
    elif breakage == 'bad_detected':
        rewrite_table(path, 'intensities', lambda rows: rows[0].update(detected='maybe'))
    elif breakage == 'bad_polarity':
        _edit_comparison(path, lambda c: c.update(polarity='pos'))
    elif breakage == 'missing_polarity_v3':
        _edit_comparison(path, lambda c: c.pop('polarity'))
    elif breakage == 'bad_intensities_stage':
        change_manifest(path, intensities_stage='imputed')
    elif breakage == 'contract_version_4':
        change_manifest(path, contract_version=4)
    elif breakage == 'qc_intensity_row':
        rewrite_table(path, 'intensities', lambda rows: rows.append(dict(rows[0], sample_id='QC1')))
    else:
        raise AssertionError(breakage)


def test_v3_loads_and_keeps_adduct(tmp_path):
    bundle = load_bundle(make_v3_bundle(tmp_path))
    assert bundle.manifest['comparison']['polarity'] == 'positive'
    assert {f.feature_id: f.adduct for f in bundle.features}['k_parent'] is None


def test_v2_still_loads(tmp_path):
    bundle = load_bundle(make_v3_bundle(tmp_path, version=2))
    assert role_path(bundle.manifest, 'intensities') is None


def test_parse_intensities_round_trip(tmp_path):
    path = make_v3_bundle(tmp_path)
    manifest = json.loads(path.read_text(encoding='utf-8'))
    samples = parse_samples((tmp_path / role_path(manifest, 'samples')).read_text(encoding='utf-8-sig'))
    matrix = parse_intensities((tmp_path / role_path(manifest, 'intensities')).read_text(encoding='utf-8-sig'),
                               {'k_parent', 'k_hnk', 'k_frag', 'x_msms', 'mass_only', 'endo', 'no_mz'}, samples)
    assert matrix['k_parent']['T3'] == (300.0, True)
    assert matrix['k_parent']['R1'] == (1.0, False)
    assert 'QC1' not in matrix['k_parent']


def test_unknown_detection_is_none(tmp_path):
    path = make_v3_bundle(tmp_path, detection='unknown')
    manifest = json.loads(path.read_text(encoding='utf-8'))
    assert manifest['population']['detection_unavailable_reason'] == 'sample_axis_mismatch'
    samples = parse_samples((tmp_path / 'samples.tsv').read_text(encoding='utf-8-sig'))
    matrix = parse_intensities((tmp_path / 'intensities.tsv').read_text(encoding='utf-8-sig'),
                               {'k_parent', 'k_hnk', 'k_frag', 'x_msms', 'mass_only', 'endo', 'no_mz'}, samples)
    assert matrix['k_parent']['T3'] == (300.0, None)
    load_bundle(path)


@pytest.mark.parametrize('breakage', [
    'samples_without_intensities', 'v2_with_samples', 'bad_group', 'n_b_mismatch',
    'missing_intensity_row', 'duplicate_intensity_row', 'unknown_sample', 'bad_detected',
    'bad_polarity', 'missing_polarity_v3', 'bad_intensities_stage', 'contract_version_4',
    'qc_intensity_row'])
def test_v3_contract_violations_raise(tmp_path, breakage):
    path = make_v3_bundle(tmp_path)
    break_bundle(path, breakage)
    with pytest.raises(BundleError):
        load_bundle(path)


def test_store_keeps_intensities_source(tmp_path):
    store = EvidenceStore(tmp_path / 'store')
    bundle = load_bundle(make_v3_bundle(tmp_path / 'in'))
    aid = store.register(bundle)
    raw = store.read_source(aid, role_path(bundle.manifest, 'intensities'))
    assert raw.decode('utf-8-sig').splitlines()[0].split('\t') == INTENSITY_COLUMNS


def test_read_source_missing_raises(tmp_path):
    store = EvidenceStore(tmp_path / 'store')
    aid = store.register(load_bundle(make_v3_bundle(tmp_path / 'in', version=2)))
    with pytest.raises(ValueError, match='入力ファイルが保存されていません: intensities.tsv'):
        store.read_source(aid, 'intensities.tsv')
