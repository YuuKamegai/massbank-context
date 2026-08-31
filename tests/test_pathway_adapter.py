"""PubChem パスウェイのアダプタ。HTTP 層は差し替える。"""
import pytest

from massbank_context.kg.pathway import (
    PubChemPathwayAdapter,
    UpstreamError,
    build_pathway_query,
)

ENDPOINT = "https://rdfportal.org/pubchem/sparql"


def _binding(inchikey, pathway="PWID1", label="Citric Acid Cycle",
             organism="http://rdf.ncbi.nlm.nih.gov/pubchem/taxonomy/TAXID10090"):
    row = {
        "value_inchikey": {"type": "literal", "value": inchikey},
        "pathway": {"type": "uri",
                    "value": f"http://rdf.ncbi.nlm.nih.gov/pubchem/pathway/{pathway}"},
        "pathway_label": {"type": "literal", "value": label},
    }
    if organism:
        row["pathway_organism"] = {"type": "uri", "value": organism}
    return row


class TestQuery:
    def test_includes_all_inchikeys_and_pathway_graph(self):
        query = build_pathway_query(["A", "B"])
        assert '"A"' in query and '"B"' in query
        assert "rdf.ncbi.nlm.nih.gov/pubchem/pathway" in query
        assert "RO_0000057" in query

    def test_rejects_empty(self):
        with pytest.raises(ValueError):
            build_pathway_query([])


class TestFetch:
    def test_parses_bindings(self):
        def transport(endpoint, query):
            return {"results": {"bindings": [_binding("A")]}}

        adapter = PubChemPathwayAdapter(ENDPOINT, transport=transport)
        got = adapter.fetch(["A"])
        assert got["A"][0]["pathway_id"] == "PWID1"
        assert got["A"][0]["pathway_label"] == "Citric Acid Cycle"
        assert got["A"][0]["taxid"] == "TAXID10090"
        assert adapter.failed == []

    def test_missing_organism_becomes_empty_taxid(self):
        def transport(endpoint, query):
            return {"results": {"bindings": [_binding("A", organism=None)]}}

        adapter = PubChemPathwayAdapter(ENDPOINT, transport=transport)
        assert adapter.fetch(["A"])["A"][0]["taxid"] == ""

    def test_batch_failure_retries_individually(self):
        """バッチが落ちたら 1 件ずつ引き直し、どの InChIKey で落ちたかを特定する。"""
        calls = []

        def transport(endpoint, query):
            calls.append(query)
            if '"A"' in query and '"B"' in query:
                raise UpstreamError("batch failed")
            if '"B"' in query:
                raise UpstreamError("B is bad")
            return {"results": {"bindings": [_binding("A")]}}

        adapter = PubChemPathwayAdapter(ENDPOINT, transport=transport)
        got = adapter.fetch(["A", "B"])
        assert list(got) == ["A"]
        assert adapter.failed == ["B"]
        assert len(calls) == 3   # バッチ 1 回 + 個別 2 回

    def test_single_inchikey_failure_is_not_retried_twice(self):
        def transport(endpoint, query):
            raise UpstreamError("down")

        adapter = PubChemPathwayAdapter(ENDPOINT, transport=transport)
        assert adapter.fetch(["A"]) == {}
        assert adapter.failed == ["A"]
