"""事実と仮説の分離（spec §10）。"""
import json

import pytest

from massbank_context.core import config
from massbank_context.core import session as session_module
from massbank_context.kg import cache
from massbank_context.report.findings import (
    INTERPRETATION_WORDS,
    build_findings,
)
from massbank_context.tools.pathway_tools import load_differential
from massbank_context.tools.report_tools import draft_pathway_findings


def _facts():
    return {
        "comparison": {"group_a": "gf_AIN", "group_b": "gf_HFD",
                       "n_a": "5", "n_b": "5",
                       "log2fc_sign": "正 = gf_HFD が高い（上昇）"},
        "population": {"n_features_total": 20, "n_with_inchikey": 6,
                       "n_unannotated": 14, "n_background": 6,
                       "n_no_pathway": 0},
        "provenance": {"endpoint": "https://rdfportal.org/pubchem/sparql",
                       "fetched_at": ["2026-08-31T00:00:00+00:00", None],
                       "organism": "TAXID10090"},
    }


def _rows():
    return [{"pathway_id": "PWID_UP", "pathway_label": "Up pathway",
             "direction": "UP_DOMINANT", "n_members": 3}]


def test_limits_always_has_five_items():
    findings = build_findings(_facts(), _rows(), coverage=0.5)
    assert len(findings["limits"]) == 5


def test_coverage_is_substituted_into_limits():
    findings = build_findings(_facts(), _rows(), coverage=0.145)
    assert "14.5%" in findings["limits"][0]


def test_hypothesis_slot_is_empty_and_requires_evidence():
    findings = build_findings(_facts(), _rows(), coverage=0.5)
    slot = findings["hypotheses"][0]
    assert slot["claim"] is None
    assert slot["evidence_rows"] == []
    assert "evidence_rows" in slot["required"]


def test_facts_contain_no_interpretation_words():
    """事実側に解釈語が混入しないことを機械的に検査する。

    「活性化」は hypotheses の中にしか置けない。facts に混ざると、
    計算していない主張が計算結果の顔をして出ていく。
    """
    findings = build_findings(_facts(), _rows(), coverage=0.5)
    text = json.dumps(findings["facts"], ensure_ascii=False)
    for word in INTERPRETATION_WORDS:
        assert word not in text


def test_top_pathways_are_carried_verbatim():
    findings = build_findings(_facts(), _rows(), coverage=0.5)
    assert findings["facts"]["top_pathways"] == _rows()


_HEADER = """# contract_version = 1
# group_a = gf_AIN\tn_a = 5
# group_b = gf_HFD\tn_b = 5
# log2fc_sign = positive means group_b is higher
# q_threshold = 0.05\tlog2fc_threshold = 1.0
# n_features_total = 20\tn_with_inchikey = 6\tn_unannotated = 14
spot_id\tname\tname_source\tontology\tinchikey\tinchikey_source\tmsi_level\tmz\trt\tlog2fc\tp_value\tq_value\tmean_a\tmean_b\tsignificant
"""


@pytest.fixture
def prepared(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "PATHWAY_CACHE_PATH",
                        tmp_path / "pathway_cache.sqlite")
    session_module.session.reset()
    lines = [_HEADER.rstrip("\n")]
    spec = [("1", "IK1", "2.0", "true"), ("2", "IK2", "2.0", "true"),
            ("3", "IK3", "2.0", "true"), ("4", "IK4", "0.1", "false"),
            ("5", "IK5", "0.1", "false"), ("6", "IK6", "0.1", "false")]
    for spot, ik, fc, sig in spec:
        lines.append("\t".join([spot, f"C{spot}", "arf2", "", ik, "arf2", "",
                                "100.0", "1.0", fc, "0.001", "0.01",
                                "1", "2", sig]))
    path = tmp_path / "differential.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    load_differential(str(path))
    conn = cache.open_cache(config.PATHWAY_CACHE_PATH)
    for ik in ("IK1", "IK2", "IK3", "IK4", "IK5", "IK6"):
        cache.record_result(conn, ik, config.PUBCHEM_SPARQL_ENDPOINT, [
            {"pathway_id": "PWID_UP", "pathway_label": "Up pathway",
             "taxid": "TAXID10090"}])
    conn.close()
    yield
    session_module.session.reset()


def test_draft_returns_three_blocks(prepared):
    payload = json.loads(draft_pathway_findings(organism="TAXID10090"))
    assert set(payload) == {"facts", "hypotheses", "limits"}
    assert payload["facts"]["comparison"]["group_b"] == "gf_HFD"
    assert payload["facts"]["provenance"]["organism"] == "TAXID10090"
    assert len(payload["limits"]) == 5


def test_draft_propagates_errors(prepared):
    payload = json.loads(draft_pathway_findings())
    assert payload["error"]["code"] == "incompatible_organism"


def test_draft_provenance_carries_endpoint_and_fetched_at(prepared):
    """コントローラ裁定: provenance は pathway_activity の payload をそのまま写す。

    endpoint と fetched_at をローカルで再計算せず、pathway_activity の
    payload の値と一致させる（Task 7 で計算済みの値を二重実装しない）。
    """
    payload = json.loads(draft_pathway_findings(organism="TAXID10090"))
    provenance = payload["facts"]["provenance"]
    assert provenance["endpoint"] == config.PUBCHEM_SPARQL_ENDPOINT
    assert provenance["fetched_at"][0] is not None
    assert provenance["fetched_at"][1] is not None
