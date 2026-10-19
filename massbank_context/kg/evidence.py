"""候補別にRDF取得結果と来歴を結び付ける。"""
import sqlite3
from datetime import datetime, timezone
from massbank_context.core import config
from massbank_context.core.evidence_ids import stable_id, canonical_json
from massbank_context.experiment.bundle import INCHIKEY_PATTERN
from massbank_context.evidence.selection import candidate_dicts
from massbank_context.evidence.matching import select_ids
from massbank_context.evidence.catalog import read_candidate_catalog
from massbank_context.kg.evidence_cache import EvidenceCache
from massbank_context.kg.pathway import PubChemPathwayAdapter, build_pathway_query

QUERY_VERSION = 'pubchem-participant-v1'
PARTICIPANT = 'http://purl.obolibrary.org/obo/RO_0000057'


class EvidenceAdapter(PubChemPathwayAdapter):
    def _query(self, inchikeys):
        payload = self._transport(self.endpoint, build_pathway_query(inchikeys))
        result = {}
        for binding in payload.get('results', {}).get('bindings', []):
            uri = binding['pathway']['value']
            result.setdefault(binding['value_inchikey']['value'], []).append(dict(
                pathway_id=uri.rsplit('/', 1)[-1], pathway_uri=uri,
                pathway_label=binding.get('pathway_label', {}).get('value', ''),
                organism_uri=binding.get('pathway_organism', {}).get('value'),
                taxid=binding.get('pathway_organism', {}).get('value', '').rsplit('/', 1)[-1],
                predicate=PARTICIPANT))
        return result


def fetch_pathway_evidence(adapter, inchikeys):
    results = {}
    hits = adapter.fetch(inchikeys)
    failed = set(adapter.failed)
    for ik in inchikeys:
        entries = hits.get(ik, [])
        results[ik] = dict(inchikey=ik, endpoint=adapter.endpoint,
                           query_version=QUERY_VERSION,
                           fetched_at=datetime.now(timezone.utc).isoformat(),
                           status='failed' if ik in failed else 'resolved' if entries else 'empty',
                           hits=entries)
    return results


def resolve_candidates(store, analysis_id, candidate_ids, *, refresh=False, adapter=None):
    bundle = store.load(analysis_id)
    candidates = candidate_dicts(bundle, store.current(analysis_id, 'matches'))
    by_id = {c['candidate_id']: c for c in candidates}
    ids = select_ids(by_id, candidate_ids)
    adapter = adapter or EvidenceAdapter(config.PUBCHEM_SPARQL_ENDPOINT)
    cache = EvidenceCache(store.root / 'rdf-cache.sqlite')
    keys = sorted({by_id[cid]['inchikey'] for cid in ids if by_id[cid]['inchikey']
                   and INCHIKEY_PATTERN.fullmatch(by_id[cid]['inchikey'])})
    snapshots, pending = {}, []
    for ik in keys:
        cached = cache.lookup(adapter.endpoint, QUERY_VERSION, ik)
        if not refresh and cached and cached['status'] in ('resolved', 'empty'):
            snapshots[ik] = cached
        else:
            pending.append(ik)
    for offset in range(0, len(pending), 50):
        batch = pending[offset:offset+50]
        got = fetch_pathway_evidence(adapter, batch)
        for ik, snap in got.items():
            cache.put(snap)
            snapshots[ik] = snap
        if all(s['status'] == 'failed' for s in got.values()):
            break
    # 違うendpointや消えた候補の現在結果を引き継がない。
    previous = store.current(analysis_id, 'lookups')
    lookups = [x for x in previous.get('lookups', []) if x['candidate_id'] in by_id
               and x['candidate_id'] not in ids and x.get('endpoint') == adapter.endpoint]
    evidence = [x for x in previous.get('evidence', [])
                if x['candidate_id'] in {l['candidate_id'] for l in lookups}]
    for cid in ids:
        candidate = by_id[cid]
        ik = candidate['inchikey']
        snap = snapshots.get(ik, dict(inchikey=ik, endpoint=adapter.endpoint,
                                      query_version=QUERY_VERSION, fetched_at=None, hits=[],
                                      status='not_queried' if ik else 'unsupported_identifier'))
        lookup_id = stable_id('lookup', [cid, canonical_json(snap)])
        lookup = dict(snap, candidate_id=cid, feature_id=candidate['feature_id'], lookup_id=lookup_id)
        lookups.append(lookup)
        for hit in snap['hits'] if snap['status'] == 'resolved' else []:
            eid = stable_id('evidence', [lookup_id, canonical_json(hit)])
            evidence.append(dict(hit, evidence_id=eid, candidate_id=cid,
                                 feature_id=candidate['feature_id'], lookup_id=lookup_id,
                                 source='pubchem_rdf', endpoint=snap['endpoint'],
                                 query_version=snap['query_version'], fetched_at=snap['fetched_at']))
    catalog = []
    if config.INDEX_PATH.is_file():
        with sqlite3.connect(config.INDEX_PATH.resolve().as_uri() + '?mode=ro', uri=True) as conn:
            catalog = read_candidate_catalog(conn, accessions=[c['source_id'] for c in candidates
                                                              if c['source'] == 'massbank'], inchikeys=keys)
    payload = dict(lookups=lookups, evidence=evidence, catalog=catalog,
                   scope=['massbank_metadata', 'pubchem_pathways'])
    store.put_snapshot(analysis_id, 'lookups', payload)
    return dict(n_resolved=sum(x['status'] == 'resolved' for x in lookups),
                n_failed=sum(x['status'] == 'failed' for x in lookups),
                n_remaining=sum(x['status'] == 'not_queried' for x in lookups),
                evidence_revision=store.latest_revision(analysis_id))
