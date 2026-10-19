"""MCP 登録面を固定する（登録漏れの早期検出）。"""
import asyncio

from massbank_context import server

EXPECTED_TOOLS = sorted([
    "load_evidence_bundle",
    "set_analysis_context",
    "match_feature_spectra",
    "resolve_candidate_evidence",
    "get_comparison_evidence",
    "draft_comparison_findings",
    "draft_pathway_findings",
    "index_status",
    "inspect_spectrum",
    "load_differential",
    "load_spectra",
    "match_spectra",
    "pathway_activity",
    "resolve_pathways",
])


def test_registered_tools():
    tools = asyncio.run(server.mcp.list_tools())
    assert sorted(t.name for t in tools) == EXPECTED_TOOLS


def test_no_structured_output():
    """structured_output=False が全ツールに効いていること（content の二重送出を防ぐ）。"""
    tools = asyncio.run(server.mcp.list_tools())
    assert all(getattr(t, "outputSchema", None) is None for t in tools)
