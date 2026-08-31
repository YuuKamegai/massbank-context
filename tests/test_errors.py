"""エラーエンベロープ（spec §8・§9）。"""
import json

import pytest

from massbank_context.core import errors


def test_contract_mismatch_names_what_was_expected():
    payload = json.loads(errors.contract_mismatch(
        reason="必須列がありません: inchikey",
        expected="contract_version = 1",
        path="C:/tmp/differential.tsv"))
    assert payload["error"]["code"] == "contract_mismatch"
    assert payload["error"]["expected"] == "contract_version = 1"
    assert payload["error"]["path"].endswith("differential.tsv")


def test_incompatible_organism_lists_choices():
    payload = json.loads(errors.incompatible_organism(
        available=[{"taxid": "TAXID9606", "n_compounds": 67}],
        message="organism を指定してください。"))
    assert payload["error"]["code"] == "incompatible_organism"
    assert payload["error"]["available"][0]["taxid"] == "TAXID9606"


def test_upstream_unavailable_names_endpoint():
    payload = json.loads(errors.upstream_unavailable(
        endpoint="https://rdfportal.org/pubchem/sparql",
        message="接続できません。"))
    assert payload["error"]["code"] == "upstream_unavailable"
    assert "rdfportal" in payload["error"]["endpoint"]


@pytest.mark.parametrize("call", [
    lambda: errors.contract_mismatch("", "x", "y"),
    lambda: errors.incompatible_organism([], ""),
    lambda: errors.upstream_unavailable("", "m"),
])
def test_required_arguments(call):
    """手掛かりの無いエンベロープを黙って出さない。"""
    with pytest.raises(ValueError):
        call()
