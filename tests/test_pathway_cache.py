"""パスウェイキャッシュ（spec §6.2）。"""
import pytest

from massbank_context.kg import cache

ENDPOINT = "https://rdfportal.org/pubchem/sparql"


@pytest.fixture
def conn(tmp_path):
    connection = cache.open_cache(tmp_path / "pathway_cache.sqlite")
    yield connection
    connection.close()


def _hit(pathway_id="PWID1", label="Citric Acid Cycle", taxid="TAXID10090"):
    return {"pathway_id": pathway_id, "pathway_label": label, "taxid": taxid}


class TestPending:
    def test_all_pending_at_first(self, conn):
        assert cache.pending_inchikeys(conn, ["A", "B"], ENDPOINT) == ["A", "B"]

    def test_recorded_hit_is_not_pending(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [_hit()])
        assert cache.pending_inchikeys(conn, ["A", "B"], ENDPOINT) == ["B"]

    def test_empty_result_is_not_pending(self, conn):
        """「照会したがヒット 0」は照会済みである。

        これを未照会と混同すると、被覆の分母が壊れる。
        """
        cache.record_result(conn, "A", ENDPOINT, [])
        assert cache.pending_inchikeys(conn, ["A"], ENDPOINT) == []

    def test_failure_stays_pending(self, conn):
        cache.record_failure(conn, "A", ENDPOINT)
        assert cache.pending_inchikeys(conn, ["A"], ENDPOINT) == ["A"]
        assert cache.failed_inchikeys(conn, ["A"]) == ["A"]

    def test_refresh_makes_everything_pending(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [_hit()])
        assert cache.pending_inchikeys(conn, ["A"], ENDPOINT,
                                       refresh=True) == ["A"]

    def test_other_endpoint_is_pending(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [_hit()])
        assert cache.pending_inchikeys(conn, ["A"], "http://other") == ["A"]


class TestMemberships:
    def test_returns_only_requested(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [_hit()])
        cache.record_result(conn, "B", ENDPOINT, [_hit(pathway_id="PWID2")])
        got = cache.memberships(conn, ["A"])
        assert list(got) == ["A"]
        assert got["A"][0]["pathway_id"] == "PWID1"

    def test_filters_by_taxid(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [
            _hit(pathway_id="PWID1", taxid="TAXID10090"),
            _hit(pathway_id="PWID2", taxid="TAXID562"),
        ])
        got = cache.memberships(conn, ["A"], taxid="TAXID10090")
        assert [e["pathway_id"] for e in got["A"]] == ["PWID1"]

    def test_refresh_replaces_previous_hits(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [_hit(pathway_id="PWID1")])
        cache.record_result(conn, "A", ENDPOINT, [_hit(pathway_id="PWID2")])
        got = cache.memberships(conn, ["A"])
        assert [e["pathway_id"] for e in got["A"]] == ["PWID2"]


class TestOrganismCounts:
    def test_counts_compounds_not_rows(self, conn):
        """同じ化合物が同じ生物種の複数パスウェイに属しても 1 と数える。"""
        cache.record_result(conn, "A", ENDPOINT, [
            _hit(pathway_id="PWID1", taxid="TAXID10090"),
            _hit(pathway_id="PWID2", taxid="TAXID10090"),
        ])
        cache.record_result(conn, "B", ENDPOINT, [
            _hit(pathway_id="PWID3", taxid="TAXID562")])
        counts = cache.organism_counts(conn, ["A", "B"])
        by_taxid = {c["taxid"]: c["n_compounds"] for c in counts}
        assert by_taxid == {"TAXID10090": 1, "TAXID562": 1}

    def test_labels_are_filled_when_known(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [_hit(taxid="TAXID10090")])
        cache.put_organism_labels(conn, {"TAXID10090": "Mus musculus"})
        assert cache.organism_counts(conn, ["A"])[0]["label"] == "Mus musculus"


class TestCacheMeta:
    def test_record_result_writes_endpoint_to_cache_meta(self, conn):
        """cache_meta に endpoint も残す（spec §6.2、来歴として必要）。"""
        cache.record_result(conn, "A", ENDPOINT, [_hit()])
        row = conn.execute(
            "SELECT value FROM cache_meta WHERE key = 'endpoint'").fetchone()
        assert row[0] == ENDPOINT

    def test_record_failure_writes_endpoint_to_cache_meta(self, conn):
        cache.record_failure(conn, "A", "http://other")
        row = conn.execute(
            "SELECT value FROM cache_meta WHERE key = 'endpoint'").fetchone()
        assert row[0] == "http://other"

    def test_cache_meta_keeps_last_used_endpoint(self, conn):
        cache.record_result(conn, "A", ENDPOINT, [_hit()])
        cache.record_failure(conn, "B", "http://other")
        row = conn.execute(
            "SELECT value FROM cache_meta WHERE key = 'endpoint'").fetchone()
        assert row[0] == "http://other"
