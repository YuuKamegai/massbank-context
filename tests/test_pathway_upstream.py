"""腐敗防止テスト: 実エンドポイントの SPARQL スキーマ追従を検知する。

既定では走らない（pytest.ini の addopts = -m "not network"）。
実行するときは `-m network` を付ける。合成 fixture では上流のスキーマ変更を
検知できないため、これだけは実物を見る（spec §11）。
"""
import pytest

from massbank_context.kg.pathway import PubChemPathwayAdapter

ENDPOINT = "https://rdfportal.org/pubchem/sparql"
# クエン酸。中心代謝なのでパスウェイが必ず付く。
CITRIC_ACID = "KRKNYBCHXYNGOX-UHFFFAOYSA-N"


@pytest.mark.network
def test_real_endpoint_returns_pathways():
    adapter = PubChemPathwayAdapter(ENDPOINT)
    got = adapter.fetch([CITRIC_ACID])
    assert CITRIC_ACID in got, "クエン酸にパスウェイが付かない = スキーマが変わった疑い"
    entry = got[CITRIC_ACID][0]
    assert entry["pathway_id"].startswith("PWID")
    assert entry["pathway_label"]
