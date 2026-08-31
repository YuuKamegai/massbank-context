"""pathway_activity ツール（spec §7）。"""
import json

import pytest

from massbank_context.core import config
from massbank_context.core import session as session_module
from massbank_context.kg import cache
from massbank_context.tools.pathway_tools import load_differential, pathway_activity

HEADER = """# contract_version = 1
# group_a = gf_AIN\tn_a = 5
# group_b = gf_HFD\tn_b = 5
# log2fc_sign = positive means group_b is higher
# q_threshold = 0.05\tlog2fc_threshold = 1.0
# n_features_total = 20\tn_with_inchikey = 6\tn_unannotated = 14
spot_id\tname\tname_source\tontology\tinchikey\tinchikey_source\tmsi_level\tmz\trt\tlog2fc\tp_value\tq_value\tmean_a\tmean_b\tsignificant
"""

ROWS = [
    ("1", "IK1", "2.0", "true"), ("2", "IK2", "2.0", "true"),
    ("3", "IK3", "2.0", "true"), ("4", "IK4", "0.1", "false"),
    ("5", "IK5", "0.1", "false"), ("6", "IK6", "0.1", "false"),
]


def _tsv():
    lines = [HEADER.rstrip("\n")]
    for spot, ik, fc, sig in ROWS:
        lines.append("\t".join([spot, f"C{spot}", "arf2", "", ik, "arf2", "",
                                "100.0", "1.0", fc, "0.001", "0.01",
                                "1", "2", sig]))
    return "\n".join(lines) + "\n"


@pytest.fixture(autouse=True)
def prepared(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "PATHWAY_CACHE_PATH",
                        tmp_path / "pathway_cache.sqlite")
    session_module.session.reset()
    path = tmp_path / "differential.tsv"
    path.write_text(_tsv(), encoding="utf-8")
    load_differential(str(path))
    conn = cache.open_cache(config.PATHWAY_CACHE_PATH)
    endpoint = config.PUBCHEM_SPARQL_ENDPOINT
    for ik in ("IK1", "IK2", "IK3"):
        cache.record_result(conn, ik, endpoint, [
            {"pathway_id": "PWID_UP", "pathway_label": "Up pathway",
             "taxid": "TAXID10090"}])
    for ik in ("IK4", "IK5", "IK6"):
        cache.record_result(conn, ik, endpoint, [
            {"pathway_id": "PWID_FLAT", "pathway_label": "Flat pathway",
             "taxid": "TAXID10090"}])
    conn.close()
    yield
    session_module.session.reset()


def test_requires_organism():
    payload = json.loads(pathway_activity())
    assert payload["error"]["code"] == "incompatible_organism"
    assert payload["error"]["available"][0]["taxid"] == "TAXID10090"


def test_rejects_unknown_organism():
    payload = json.loads(pathway_activity(organism="TAXID9999"))
    assert payload["error"]["code"] == "incompatible_organism"


def test_returns_rows_and_header():
    payload = json.loads(pathway_activity(organism="TAXID10090"))
    assert payload["n_background"] == 6
    assert payload["n_no_pathway"] == 0
    assert payload["n_pathways_total"] == 2
    assert payload["group_a"] == "gf_AIN"
    assert "上昇" in payload["log2fc_sign"]
    lines = payload["candidates"].splitlines()
    assert lines[0].split("\t")[0] == "pathway_id"
    assert lines[1].split("\t")[0] == "PWID_UP"
    assert "UP_DOMINANT" in lines[1]


def test_blocks_when_lookup_incomplete():
    conn = cache.open_cache(config.PATHWAY_CACHE_PATH)
    cache.record_failure(conn, "IK1", config.PUBCHEM_SPARQL_ENDPOINT)
    conn.close()
    payload = json.loads(pathway_activity(organism="TAXID10090"))
    assert payload["error"]["code"] == "missing_state"
    assert "resolve_pathways" in payload["error"]["required_tools"]


def test_reports_truncation():
    payload = json.loads(pathway_activity(organism="TAXID10090", top_k=1))
    assert payload["truncated"] is True
    assert payload["n_pathways_total"] == 2
    assert len(payload["candidates"].splitlines()) == 2   # ヘッダ + 1 行


def test_facts_for_presentation_layer():
    """spec §10 の facts が要求する 3 点: pathway_coverage / endpoint / fetched_at。"""
    payload = json.loads(pathway_activity(organism="TAXID10090"))
    # 6 化合物全てにパスウェイが付いている（n_no_pathway == 0）ので被覆は 1.0。
    assert payload["pathway_coverage"] == 1.0
    assert payload["endpoint"] == config.PUBCHEM_SPARQL_ENDPOINT
    fetched_at = payload["fetched_at"]
    assert len(fetched_at) == 2
    assert fetched_at[0] is not None and fetched_at[1] is not None
    assert fetched_at[0] <= fetched_at[1]


def test_pathway_coverage_excludes_no_pathway_compounds():
    """パスウェイが 1 つも無い化合物がいる場合、被覆はその分だけ下がる。"""
    conn = cache.open_cache(config.PATHWAY_CACHE_PATH)
    # IK4 を「照会したがヒット 0」に上書きする（pending には含めない）。
    cache.record_result(conn, "IK4", config.PUBCHEM_SPARQL_ENDPOINT, [])
    conn.close()
    payload = json.loads(pathway_activity(organism="TAXID10090"))
    assert payload["n_no_pathway"] == 1
    assert payload["n_background"] == 5
    assert payload["pathway_coverage"] == pytest.approx(5 / 6)


def _tsv_from(rows):
    """任意の (spot_id, inchikey, log2fc, significant) から契約 TSV を組む。"""
    lines = [HEADER.rstrip("\n")]
    for spot, ik, fc, sig in rows:
        lines.append("\t".join([spot, f"C{spot}", "arf2", "", ik, "arf2", "",
                                "100.0", "1.0", fc, "0.001", "0.01",
                                "1", "2", sig]))
    return "\n".join(lines) + "\n"


def _candidates(payload):
    """candidates の TSV を pathway_id -> 列 dict に開く（検証用）。"""
    lines = payload["candidates"].splitlines()
    columns = lines[0].split("\t")
    rows = [dict(zip(columns, line.split("\t"))) for line in lines[1:]]
    return {row["pathway_id"]: row for row in rows}


def _load_with_decoy_spots(tmp_path):
    """有意な行と、それより |log2fc| が大きい非有意な行を同じ InChIKey に持たせる。"""
    rows = [
        ("1", "IK1", "2.0", "true"), ("7", "IK1", "-3.0", "false"),
        ("2", "IK2", "2.0", "true"), ("8", "IK2", "-3.5", "false"),
        ("3", "IK3", "2.0", "true"), ("4", "IK4", "0.1", "false"),
        ("5", "IK5", "0.1", "false"), ("6", "IK6", "0.1", "false"),
    ]
    path = tmp_path / "collapsed.tsv"
    path.write_text(_tsv_from(rows), encoding="utf-8")
    load_differential(str(path))


def test_collapse_prefers_significant_spot(tmp_path):
    """同一 InChIKey の併合で、有意な行を非有意な行に置き換えない。

    |log2fc| だけで選ぶと IK1 は -3.0/非有意、IK2 は -3.5/非有意が採られ、
    n_significant が 3 から 1 に落ちる。Fisher の k と n が動いて濃縮の
    p 値まで変わり、direction は UP_DOMINANT から SINGLE_MEMBER に化ける。
    significant は上流 ms-data-parser の閾値判定であり、こちらで取り消さない。
    """
    _load_with_decoy_spots(tmp_path)
    payload = json.loads(pathway_activity(organism="TAXID10090"))
    assert payload["n_significant"] == 3
    up = _candidates(payload)["PWID_UP"]
    assert up["n_significant"] == "3"
    assert up["n_up"] == "3"
    assert up["n_down"] == "0"
    assert up["direction"] == "UP_DOMINANT"
    assert up["median_log2fc"] == "2.0"


def test_collapse_counts_are_reported(tmp_path):
    """併合が起きたこと自体を payload に出す（黙って行を減らさない・spec §12）。"""
    _load_with_decoy_spots(tmp_path)
    payload = json.loads(pathway_activity(organism="TAXID10090"))
    assert payload["n_rows"] == 8
    assert payload["n_inchikeys"] == 6
    assert payload["n_collapsed_spots"] == 2
    assert "significant" in payload["collapse_rule"]
    assert "log2fc" in payload["collapse_rule"]


def test_collapse_counts_are_zero_without_duplicates():
    """重複が無ければ n_collapsed_spots は 0 になる。"""
    payload = json.loads(pathway_activity(organism="TAXID10090"))
    assert payload["n_rows"] == 6
    assert payload["n_inchikeys"] == 6
    assert payload["n_collapsed_spots"] == 0


def test_collapse_keeps_largest_absolute_log2fc_among_significant(tmp_path):
    """有意な行が複数あるときは、その中で |log2fc| が最大の行を採る。"""
    rows = [
        ("1", "IK1", "1.0", "true"), ("7", "IK1", "-4.0", "true"),
        ("2", "IK2", "2.0", "true"), ("3", "IK3", "2.0", "true"),
        ("4", "IK4", "0.1", "false"), ("5", "IK5", "0.1", "false"),
        ("6", "IK6", "0.1", "false"),
    ]
    path = tmp_path / "both_significant.tsv"
    path.write_text(_tsv_from(rows), encoding="utf-8")
    load_differential(str(path))
    payload = json.loads(pathway_activity(organism="TAXID10090"))
    up = _candidates(payload)["PWID_UP"]
    assert up["n_significant"] == "3"
    assert up["n_down"] == "1"          # -4.0 が採られている
    assert up["direction"] == "MIXED"


def test_wilcoxon_q_is_not_exposed():
    """wilcoxon_q は公開しない。fisher_q と並べると合成したくなる（spec §5.3）。"""
    payload = json.loads(pathway_activity(organism="TAXID10090"))
    assert "wilcoxon_q" not in payload["candidates"]
    assert "wilcoxon_q" not in json.dumps(payload, ensure_ascii=False)
