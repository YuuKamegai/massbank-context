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
from massbank_context.tools.pathway_tools import load_differential, pathway_activity
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


def test_draft_top_pathways_have_numeric_and_bool_types(prepared):
    """レビュー指摘 1: top_pathways は TSV 再パースで全部文字列になっていた。

    設計書 §10 は top_pathways を「数値のみ、解釈なし」と定める。
    build_findings に手組みの dict を渡すテストでは TSV 再パースの不具合を
    捕まえられないため、実際の draft_pathway_findings の戻り値を検査する。
    """
    payload = json.loads(draft_pathway_findings(organism="TAXID10090"))
    row = payload["facts"]["top_pathways"][0]
    # 文字列のまま残るべき列
    assert isinstance(row["pathway_id"], str)
    assert isinstance(row["pathway_label"], str)
    assert isinstance(row["taxid"], str)
    assert isinstance(row["direction"], str)
    # 数値・bool に復元されるべき列
    assert isinstance(row["n_members"], int)
    assert isinstance(row["n_significant"], int)
    assert isinstance(row["n_up"], int)
    assert isinstance(row["n_down"], int)
    assert isinstance(row["median_log2fc"], float)
    assert isinstance(row["fisher_p"], float)
    assert isinstance(row["fisher_q"], float)
    assert isinstance(row["wilcoxon_p"], float)
    assert row["tested"] is True
    assert row["wilcoxon_tested"] is True


def test_draft_population_truncated_is_carried(prepared):
    """レビュー指摘 2: truncated が黙って落ちていた（打ち切りを必ず出す）。"""
    payload = json.loads(draft_pathway_findings(organism="TAXID10090"))
    activity_payload = json.loads(pathway_activity(organism="TAXID10090"))
    assert "truncated" in payload["facts"]["population"]
    assert payload["facts"]["population"]["truncated"] == activity_payload["truncated"]


def test_draft_population_pathway_coverage_matches_payload(prepared):
    """コントローラ裁定 3 の未達分: population.pathway_coverage が payload と一致すること。"""
    payload = json.loads(draft_pathway_findings(organism="TAXID10090"))
    activity_payload = json.loads(pathway_activity(organism="TAXID10090"))
    assert (payload["facts"]["population"]["pathway_coverage"]
            == activity_payload["pathway_coverage"])


def test_draft_payload_is_strict_rfc8259_json(prepared):
    """NaN / Infinity を含まない JSON を返す（厳密なパーサで落ちないこと）。

    Python の json.loads は既定で NaN / Infinity を通してしまうので、
    parse_constant で明示的に検出する。契約リーダが非有限値を弾いているのが
    一次防御で、これはその出口側の回帰テストである。
    """
    def _reject(name):
        raise AssertionError(f"RFC 8259 に無い定数が出力されました: {name}")

    raw = draft_pathway_findings(organism="TAXID10090")
    json.loads(raw, parse_constant=_reject)
    assert "NaN" not in raw
    assert "Infinity" not in raw


def test_draft_group_sizes_are_integers(prepared):
    """facts.comparison の n_a / n_b は int（設計書 §10 のサンプルどおり）。"""
    payload = json.loads(draft_pathway_findings(organism="TAXID10090"))
    comparison = payload["facts"]["comparison"]
    assert comparison["n_a"] == 5
    assert comparison["n_b"] == 5
    assert isinstance(comparison["n_a"], int)
    assert isinstance(comparison["n_b"], int)


def test_draft_does_not_expose_wilcoxon_q(prepared):
    """wilcoxon_q は top_pathways にも出さない。"""
    payload = json.loads(draft_pathway_findings(organism="TAXID10090"))
    assert "wilcoxon_q" not in payload["facts"]["top_pathways"][0]


def test_draft_reports_collapsed_spots(prepared):
    """併合件数を提示層でも出す（黙って行を減らさない）。"""
    payload = json.loads(draft_pathway_findings(organism="TAXID10090"))
    assert payload["facts"]["population"]["n_collapsed_spots"] == 0


def test_draft_survives_tab_in_pathway_label(tmp_path, monkeypatch):
    """pathway_label にタブが含まれても top_pathways が壊れない。

    TSV へ落として再パースしていたときは as_tsv がエスケープしないため列が
    ずれ、facts.top_pathways が黙って壊れていた。共有した内部関数から
    typed dict を直接受け取るので、往復そのものが無くなっている。
    """
    monkeypatch.setattr(config, "PATHWAY_CACHE_PATH",
                        tmp_path / "pathway_cache.sqlite")
    session_module.session.reset()
    lines = [_HEADER.rstrip("\n")]
    spec = [("1", "IK1", "2.0", "true"), ("2", "IK2", "2.0", "true"),
            ("3", "IK3", "2.0", "true")]
    for spot, ik, fc, sig in spec:
        lines.append("\t".join([spot, f"C{spot}", "arf2", "", ik, "arf2", "",
                                "100.0", "1.0", fc, "0.001", "0.01",
                                "1", "2", sig]))
    path = tmp_path / "differential.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    load_differential(str(path))
    conn = cache.open_cache(config.PATHWAY_CACHE_PATH)
    for ik in ("IK1", "IK2", "IK3"):
        cache.record_result(conn, ik, config.PUBCHEM_SPARQL_ENDPOINT, [
            {"pathway_id": "PWID_TAB", "pathway_label": "Up\tpathway",
             "taxid": "TAXID10090"}])
    conn.close()
    try:
        payload = json.loads(draft_pathway_findings(organism="TAXID10090"))
        row = payload["facts"]["top_pathways"][0]
        assert row["pathway_id"] == "PWID_TAB"
        assert row["pathway_label"] == "Up\tpathway"
        assert row["taxid"] == "TAXID10090"
        assert row["n_members"] == 3
    finally:
        session_module.session.reset()
