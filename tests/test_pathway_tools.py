"""パスウェイ系ツール。"""
import json

import pytest

from massbank_context.core import session as session_module
from massbank_context.tools.pathway_tools import load_differential

GOOD = """# contract_version = 1
# exported_at = 2026-08-31T12:00:00+00:00
# group_a = gf_AIN\tn_a = 5
# group_b = gf_HFD\tn_b = 5
# log2fc_sign = positive means group_b is higher
# q_threshold = 0.05\tlog2fc_threshold = 1.0\tlog_transform = true
# n_features_total = 10\tn_with_inchikey = 2\tn_unannotated = 8
spot_id\tname\tname_source\tontology\tinchikey\tinchikey_source\tmsi_level\tmz\trt\tlog2fc\tp_value\tq_value\tmean_a\tmean_b\tsignificant
1\tPC 34:1\tarf2\tPC\tAAAAAAAAAAAAAA-BBBBBBBBBB-C\tarf2\t\t760.5851\t12.34\t1.900000\t0.001\t0.01\t10\t40\ttrue
2\tCitric acid\tarf2\t\tKRKNYBCHXYNGOX-UHFFFAOYSA-N\tarf2\t\t191.0197\t2.10\t-0.500000\t0.4\t0.6\t20\t14\tfalse
"""


@pytest.fixture(autouse=True)
def clean_session():
    """各テストを独立した可変セッションで実行する。"""
    session_module.session.reset()
    yield
    session_module.session.reset()


def test_load_differential_reports_counts(tmp_path):
    """読み込み結果は行データではなく件数と比較の向きを返す。"""
    path = tmp_path / "differential.tsv"
    path.write_text(GOOD, encoding="utf-8")
    payload = json.loads(load_differential(str(path)))
    assert payload["n_rows"] == 2
    assert payload["n_inchikeys"] == 2
    assert payload["n_unannotated"] == 8
    assert payload["group_a"] == "gf_AIN"
    assert payload["group_b"] == "gf_HFD"
    assert "group_b" in payload["log2fc_sign"]
    assert session_module.session.differential is not None


def test_load_differential_rejects_bad_contract(tmp_path):
    """不正な契約はセッションを更新せず契約エラーとして返す。"""
    path = tmp_path / "differential.tsv"
    path.write_text(
        GOOD.replace("contract_version = 1", "contract_version = 99"),
        encoding="utf-8",
    )
    payload = json.loads(load_differential(str(path)))
    assert payload["error"]["code"] == "contract_mismatch"
    assert session_module.session.differential is None


def test_load_differential_missing_file(tmp_path):
    """存在しないファイルは契約エラーとして返す。"""
    payload = json.loads(load_differential(str(tmp_path / "nope.tsv")))
    assert payload["error"]["code"] == "contract_mismatch"
    assert "見つかりません" in payload["error"]["reason"]


def test_payload_does_not_carry_rows(tmp_path):
    """行そのものは payload に載せない（spec §12）。"""
    path = tmp_path / "differential.tsv"
    path.write_text(GOOD, encoding="utf-8")
    payload = json.loads(load_differential(str(path)))
    assert "rows" not in payload
    assert "AAAAAAAAAAAAAA" not in json.dumps(payload)


def test_load_differential_rejects_invalid_meta_without_replacing_session(tmp_path):
    """不正な整数メタは契約エラーにし、既存の結果を置き換えない。"""
    valid_path = tmp_path / "valid.tsv"
    valid_path.write_text(GOOD, encoding="utf-8")
    load_differential(str(valid_path))
    previous = session_module.session.differential
    previous_source = session_module.session.differential_source

    invalid_path = tmp_path / "invalid.tsv"
    invalid_path.write_text(
        GOOD.replace("n_features_total = 10", "n_features_total = not-a-number"),
        encoding="utf-8",
    )
    payload = json.loads(load_differential(str(invalid_path)))

    assert payload["error"]["code"] == "contract_mismatch"
    assert "n_features_total" in payload["error"]["reason"]
    assert session_module.session.differential is previous
    assert session_module.session.differential_source == previous_source
