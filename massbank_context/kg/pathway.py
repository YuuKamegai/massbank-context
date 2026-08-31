"""PubChem のパスウェイを SPARQL で引く（spec §5.2 の SparqlAdapter 実装）。

obo:RO_0000057 は **has participant** であり、その化合物が基質か産物か・
反応の向き・量論を持たない。ここで得られるのは「参加している」だけで、
活性そのものではない（spec §3.1）。

HTTP 層は transport で差し替えられるようにしてある。テストで実エンドポイントを
叩かないためと、腐敗防止テスト（別マーク）だけが実物を見るようにするため。
"""
from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request

DEFAULT_TIMEOUT = 120

_QUERY_TEMPLATE = """
PREFIX sio: <http://semanticscience.org/resource/>
PREFIX vocab: <http://rdf.ncbi.nlm.nih.gov/pubchem/vocabulary#>
PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX up: <http://purl.uniprot.org/core/>
SELECT DISTINCT ?value_inchikey ?pathway ?pathway_label ?pathway_organism
FROM <http://rdf.ncbi.nlm.nih.gov/pubchem/inchikey>
FROM <http://rdf.ncbi.nlm.nih.gov/pubchem/compound>
FROM <http://rdf.ncbi.nlm.nih.gov/pubchem/pathway>
WHERE {
  VALUES ?value_inchikey { %(values)s }
  ?inchikey_node sio:SIO_000300 ?value_inchikey ; sio:SIO_000011 ?cid .
  ?cid a vocab:Compound .
  ?pathway a vocab:Pathway ;
           obo:RO_0000057 ?cid ;
           dcterms:title ?pathway_label .
  OPTIONAL { ?pathway up:organism ?pathway_organism . }
}
"""

_ORGANISM_QUERY_TEMPLATE = """
PREFIX dcterms: <http://purl.org/dc/terms/>
SELECT DISTINCT ?taxonomy ?label
FROM <http://rdf.ncbi.nlm.nih.gov/pubchem/taxonomy>
WHERE {
  VALUES ?taxonomy { %(values)s }
  ?taxonomy dcterms:title ?label .
}
"""


class UpstreamError(RuntimeError):
    """外部エンドポイントへの照会が失敗した。"""


def _values(items: list[str], as_uri: bool = False) -> str:
    if as_uri:
        return " ".join(f"<{item}>" for item in items)
    return " ".join(json.dumps(item) for item in items)


def build_pathway_query(inchikeys: list[str]) -> str:
    """InChIKey から CID を経てパスウェイを引くクエリ。"""
    if not inchikeys:
        raise ValueError("inchikeys が空です。")
    return _QUERY_TEMPLATE % {"values": _values(inchikeys)}


def build_organism_label_query(taxids: list[str]) -> str:
    """taxonomy URI からラベルを引くクエリ。"""
    if not taxids:
        raise ValueError("taxids が空です。")
    uris = [f"http://rdf.ncbi.nlm.nih.gov/pubchem/taxonomy/{t}" for t in taxids]
    return _ORGANISM_QUERY_TEMPLATE % {"values": _values(uris, as_uri=True)}


def http_transport(endpoint: str, query: str) -> dict:
    """既定の HTTP 層。SPARQL results JSON を返す。"""
    request = urllib.request.Request(
        endpoint,
        data=urllib.parse.urlencode({"query": query}).encode(),
        headers={"Accept": "application/sparql-results+json"})
    try:
        with urllib.request.urlopen(request, timeout=DEFAULT_TIMEOUT) as response:
            return json.load(response)
    except (urllib.error.URLError, OSError, ValueError) as error:
        raise UpstreamError(str(error)) from error


def _tail(uri: str) -> str:
    return uri.rsplit("/", 1)[-1] if uri else ""


class PubChemPathwayAdapter:
    """バッチで引き、落ちたら 1 件ずつ引き直すアダプタ。"""

    def __init__(self, endpoint: str, transport=None) -> None:
        self.endpoint = endpoint
        self._transport = transport or http_transport
        self.failed: list[str] = []

    def _query(self, inchikeys: list[str]) -> dict[str, list[dict]]:
        payload = self._transport(self.endpoint, build_pathway_query(inchikeys))
        result: dict[str, list[dict]] = {}
        for binding in payload.get("results", {}).get("bindings", []):
            inchikey = binding["value_inchikey"]["value"]
            result.setdefault(inchikey, []).append({
                "pathway_id": _tail(binding["pathway"]["value"]),
                "pathway_label": binding.get("pathway_label", {}).get("value", ""),
                "taxid": _tail(binding.get("pathway_organism", {}).get("value", "")),
            })
        return result

    def fetch(self, inchikeys: list[str]) -> dict[str, list[dict]]:
        """1 バッチ分を引く。落ちたら 1 件ずつ引き直す。

        バッチ全体を諦めると、1 件の悪い InChIKey が 50 件を巻き添えにする。
        どれが悪いのかを特定できないと、再試行の対象も決められない。
        """
        if not inchikeys:
            return {}
        try:
            return self._query(inchikeys)
        except UpstreamError:
            if len(inchikeys) == 1:
                self.failed.extend(inchikeys)
                return {}
        result: dict[str, list[dict]] = {}
        for inchikey in inchikeys:
            try:
                result.update(self._query([inchikey]))
            except UpstreamError:
                self.failed.append(inchikey)
        return result

    def fetch_organism_labels(self, taxids: list[str]) -> dict[str, str]:
        """taxid のラベルを引く（best effort。失敗しても解析は続ける）。"""
        taxids = [t for t in taxids if t]
        if not taxids:
            return {}
        try:
            payload = self._transport(self.endpoint,
                                      build_organism_label_query(taxids))
        except UpstreamError:
            return {}
        return {_tail(b["taxonomy"]["value"]): b["label"]["value"]
                for b in payload.get("results", {}).get("bindings", [])}
