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
        assert cache.failed_inchikeys(conn, ["A"], ENDPOINT) == ["A"]

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


class TestFailedIsScopedToEndpoint:
    def test_failure_on_other_endpoint_is_not_counted(self, conn):
        """failed も endpoint で絞る（pending_inchikeys と揃える）。

        絞らないと、MASSBANK_CONTEXT_PUBCHEM_SPARQL を切り替えたときに旧
        エンドポイントの failed 行が n_failed に混じり、被覆の分母
        （len(inchikeys) - n_failed）が不当に縮む。
        """
        cache.record_failure(conn, "A", ENDPOINT)
        assert cache.failed_inchikeys(conn, ["A"], ENDPOINT) == ["A"]
        assert cache.failed_inchikeys(conn, ["A"], "http://other") == []


class TestTaxidIsNeverNull:
    def test_missing_taxid_is_stored_as_empty_string(self, conn):
        """taxid の無いヒットは NULL ではなく空文字列で入れる。

        PRIMARY KEY は (inchikey, pathway_id, taxid) だが、SQLite は NULL 同士を
        相異とみなすため、NULL のままだと同じヒットが重複して入る。
        """
        hit = {"pathway_id": "PWID1", "pathway_label": "No organism"}
        cache.record_result(conn, "A", ENDPOINT, [hit, dict(hit)])
        entries = cache.memberships(conn, ["A"])["A"]
        assert len(entries) == 1
        assert entries[0]["taxid"] == ""
        nulls = conn.execute(
            "SELECT COUNT(*) FROM pathway_hit WHERE taxid IS NULL").fetchone()[0]
        assert nulls == 0


class TestCacheVersion:
    def test_existing_version_mismatch_is_refused(self, tmp_path):
        """版違いのキャッシュを黙って開かない。

        CREATE TABLE IF NOT EXISTS は既存テーブルに列を足さないので、
        スキーマを v2 にしても v1 ファイルはそのまま開けてしまう。
        """
        path = tmp_path / "old_cache.sqlite"
        connection = cache.open_cache(path)
        connection.execute(
            "UPDATE cache_meta SET value = '99' WHERE key = 'contract_version'")
        connection.commit()
        connection.close()
        with pytest.raises(cache.CacheVersionError) as excinfo:
            cache.open_cache(path)
        assert "99" in str(excinfo.value)
        # どうすればよいかを書く（ファイルを消して引き直す）。
        assert "削除" in str(excinfo.value)

    def test_same_version_opens(self, tmp_path):
        path = tmp_path / "cache.sqlite"
        cache.open_cache(path).close()
        cache.open_cache(path).close()
