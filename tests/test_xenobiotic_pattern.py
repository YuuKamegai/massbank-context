"""出方の筋（群の偏り・親との連動）。"""
import csv

import pytest

from massbank_context.experiment.bundle import parse_intensities, parse_samples
from massbank_context.stats.correlation import spearman_rho
from massbank_context.xenobiotic.pattern import (
    detection_counts, group_bias, parent_linked, parent_rho, sample_sides)
from massbank_context.xenobiotic.tables import THRESHOLDS
from xenobiotic_helpers import make_v3_bundle


def _load(root, detection):
    make_v3_bundle(root, detection=detection)

    def read(name):
        return (root / name).read_text(encoding='utf-8-sig')

    with (root / 'features.tsv').open(encoding='utf-8-sig', newline='') as stream:
        ids = {row['feature_id'] for row in csv.DictReader(stream, delimiter='\t')}
    samples = parse_samples(read('samples.tsv'))
    return samples, parse_intensities(read('intensities.tsv'), ids, samples)


@pytest.fixture
def loaded(tmp_path):
    return _load(tmp_path / 'known', 'known')


@pytest.fixture
def samples(loaded):
    return loaded[0]


@pytest.fixture
def matrix(loaded):
    return loaded[1]


@pytest.fixture
def unknown_matrix(tmp_path):
    return _load(tmp_path / 'unknown', 'unknown')[1]


def test_spearman_with_ties():
    assert spearman_rho([1, 2, 3, 4, 5], [5, 6, 7, 8, 7]) == pytest.approx(0.820783, abs=1e-6)
    assert spearman_rho([1, 2, 3], [4, 4, 4]) is None
    assert spearman_rho([1], [2]) is None


def test_bias_for_test_only_feature(matrix, samples):
    sides = sample_sides(samples, 'b')
    assert sides.test == ('T1', 'T2', 'T3', 'T4', 'T5', 'T6') and 'QC1' not in sides.reference
    assert detection_counts(matrix['x_msms'], sides) == (4, 1)
    assert group_bias((4, 1), sides, 'test', THRESHOLDS) is True
    assert group_bias((6, 6), sides, 'test', THRESHOLDS) is False
    assert group_bias((0, 6), sides, 'reference', THRESHOLDS) is True
    assert group_bias((6, 0), sides, 'both', THRESHOLDS) is None


def test_test_side_a_swaps_groups(samples):
    sides = sample_sides(samples, 'a')
    assert sides.test[0] == 'R1' and sides.reference[0] == 'T1'


def test_unknown_detection_is_unavailable_not_zero(unknown_matrix, samples):
    sides = sample_sides(samples, 'b')
    assert detection_counts(unknown_matrix['k_parent'], sides) == (None, None)
    assert group_bias((None, None), sides, 'test', THRESHOLDS) is None


def test_partly_unknown_counts_only_known_detections(samples):
    sides = sample_sides(samples, 'b')
    row = {sid: (1.0, None) for sid in sides.test + sides.reference}
    row['T1'] = (1.0, True)
    assert detection_counts(row, sides) == (1, None)


def test_parent_rho_skips_imputed(matrix, samples):
    sides = sample_sides(samples, 'b')
    rho, n = parent_rho(matrix['k_hnk'], matrix['k_parent'], sides)
    assert rho == pytest.approx(1.0) and n == 6
    assert parent_linked(rho, n, THRESHOLDS)
    assert not parent_linked(1.0, 4, THRESHOLDS)
    assert not parent_linked(None, 6, THRESHOLDS)


def test_parent_rho_excludes_undetected_and_missing(samples):
    sides = sample_sides(samples, 'b')
    child = {s: (float(i), True) for i, s in enumerate(sides.test)}
    parent = {s: (float(i), True) for i, s in enumerate(sides.test)}
    parent['T1'] = (1.0, False)
    child['T2'] = (None, True)
    assert parent_rho(child, parent, sides)[1] == 4
