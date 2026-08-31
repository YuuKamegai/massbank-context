"""パスウェイ単位の集計（純関数）。"""
import pytest

from massbank_context.stats.activity import pathway_activity_rows


def _compounds(spec):
    """spec は (inchikey, log2fc, significant) のリスト。"""
    return [{"inchikey": ik, "log2fc": fc, "significant": sig}
            for ik, fc, sig in spec]


def _members(mapping, pathway_id="PWID1", label="Citric Acid Cycle",
             taxid="TAXID10090"):
    return {ik: [{"pathway_id": pathway_id, "pathway_label": label,
                  "taxid": taxid}] for ik in mapping}


class TestDirection:
    def test_up_dominant(self):
        compounds = _compounds([("A", 2.0, True), ("B", 1.5, True),
                                ("C", 0.1, False)])
        rows, _ = pathway_activity_rows(
            compounds, _members(["A", "B", "C"]), min_members=3)
        assert rows[0]["direction"] == "UP_DOMINANT"
        assert rows[0]["n_up"] == 2
        assert rows[0]["n_down"] == 0

    def test_down_dominant(self):
        compounds = _compounds([("A", -2.0, True), ("B", -1.5, True),
                                ("C", 0.1, False)])
        rows, _ = pathway_activity_rows(
            compounds, _members(["A", "B", "C"]), min_members=3)
        assert rows[0]["direction"] == "DOWN_DOMINANT"

    def test_mixed(self):
        compounds = _compounds([("A", 2.0, True), ("B", -1.5, True),
                                ("C", 0.1, False)])
        rows, _ = pathway_activity_rows(
            compounds, _members(["A", "B", "C"]), min_members=3)
        assert rows[0]["direction"] == "MIXED"

    def test_single_member(self):
        compounds = _compounds([("A", 2.0, True), ("B", 0.1, False),
                                ("C", 0.1, False)])
        rows, _ = pathway_activity_rows(
            compounds, _members(["A", "B", "C"]), min_members=3)
        assert rows[0]["direction"] == "SINGLE_MEMBER"

    def test_no_significant_member(self):
        compounds = _compounds([("A", 0.2, False), ("B", 0.1, False),
                                ("C", 0.1, False)])
        rows, _ = pathway_activity_rows(
            compounds, _members(["A", "B", "C"]), min_members=3)
        assert rows[0]["direction"] == "NO_SIGNIFICANT_MEMBER"


class TestMedianUsesAllMembers:
    def test_median_includes_non_significant(self):
        """中央値は有意でないメンバーも含めた全員から取る。"""
        compounds = _compounds([("A", 3.0, True), ("B", 3.0, True),
                                ("C", 0.0, False)])
        rows, _ = pathway_activity_rows(
            compounds, _members(["A", "B", "C"]), min_members=3)
        assert rows[0]["median_log2fc"] == pytest.approx(3.0)
        assert rows[0]["n_members"] == 3


class TestTestingGate:
    def test_small_pathway_is_not_tested_but_kept(self):
        compounds = _compounds([("A", 2.0, True), ("B", 0.1, False),
                                ("C", 0.1, False)])
        memberships = {"A": [{"pathway_id": "PWID_SMALL",
                              "pathway_label": "Small", "taxid": "TAXID10090"}],
                       "B": [], "C": []}
        rows, summary = pathway_activity_rows(compounds, memberships,
                                              min_members=3)
        row = [r for r in rows if r["pathway_id"] == "PWID_SMALL"][0]
        assert row["tested"] is False
        assert row["fisher_p"] is None
        assert row["fisher_q"] is None
        assert summary["n_pathways_total"] == 1
        assert summary["n_pathways_tested"] == 0


class TestOrdering:
    def test_sorted_by_q_then_pathway_id(self):
        compounds = _compounds([("A", 2.0, True), ("B", 2.0, True),
                                ("C", 2.0, True), ("D", 0.0, False),
                                ("E", 0.0, False), ("F", 0.0, False)])
        memberships = {
            "A": [{"pathway_id": "PWID_B", "pathway_label": "B",
                   "taxid": "TAXID10090"}],
            "B": [{"pathway_id": "PWID_B", "pathway_label": "B",
                   "taxid": "TAXID10090"}],
            "C": [{"pathway_id": "PWID_B", "pathway_label": "B",
                   "taxid": "TAXID10090"}],
            "D": [{"pathway_id": "PWID_A", "pathway_label": "A",
                   "taxid": "TAXID10090"}],
            "E": [{"pathway_id": "PWID_A", "pathway_label": "A",
                   "taxid": "TAXID10090"}],
            "F": [{"pathway_id": "PWID_A", "pathway_label": "A",
                   "taxid": "TAXID10090"}],
        }
        rows, _ = pathway_activity_rows(compounds, memberships, min_members=3)
        assert rows[0]["pathway_id"] == "PWID_B"
        assert rows[0]["fisher_q"] <= rows[1]["fisher_q"]


class TestSeparateBenjaminiHochberg:
    def test_wilcoxon_q_is_corrected_separately_from_enrichment(self):
        """Wilcoxon p=[0.0625, 0.125] は別系列で BH 補正して両方 0.125。"""
        compounds = _compounds([
            ("A", 1.0, True), ("B", 2.0, True), ("C", 3.0, True),
            ("D", 4.0, True), ("E", 5.0, True),
            ("F", -1.0, False), ("G", 2.0, False), ("H", 3.0, False),
            ("I", 4.0, False), ("J", 5.0, False),
        ])
        memberships = {
            **_members(["A", "B", "C", "D", "E"], "PWID_A", "A"),
            **_members(["F", "G", "H", "I", "J"], "PWID_B", "B"),
        }
        rows, _ = pathway_activity_rows(compounds, memberships, min_members=5)
        by_id = {row["pathway_id"]: row for row in rows}
        assert by_id["PWID_A"]["wilcoxon_p"] == pytest.approx(0.0625)
        assert by_id["PWID_B"]["wilcoxon_p"] == pytest.approx(0.125)
        assert by_id["PWID_A"]["wilcoxon_q"] == pytest.approx(0.125)
        assert by_id["PWID_B"]["wilcoxon_q"] == pytest.approx(0.125)
        assert by_id["PWID_A"]["fisher_q"] != by_id["PWID_A"]["wilcoxon_q"]
