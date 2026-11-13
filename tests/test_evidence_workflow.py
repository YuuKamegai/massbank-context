"""公開MCPとローカルSPARQL応答を通す合成縦断検証。"""
import asyncio
import hashlib
import json
import os
import re
import sys
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import parse_qs
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from evidence_helpers import make_bundle, read_tsv, rewrite_table
from test_matcher import conn


def test_public_workflow_with_spectra_and_rdf(tmp_path, conn):
    path = make_bundle(tmp_path / 'input')
    msp = path.parent / 'spectra.msp'
    msp.write_text('Name: synthetic\nPrecursorMZ: 230.0552\nIonMode: Positive\nNum Peaks: 2\n86.0349 20\n171.0069 999\n')
    manifest = json.loads(path.read_text())
    manifest['files'].append(dict(role='spectra', path=msp.name, sha256=hashlib.sha256(msp.read_bytes()).hexdigest()))
    path.write_text(json.dumps(manifest))
    rewrite_table(path, 'spectrum_links', lambda rows: rows.append(dict(
        feature_id='f1', spectrum_file=msp.name, spectrum_id='spectra:1',
        source_analysis_id='synthetic-run', mapping_evidence='fixture-explicit-id')))
    received = []

    class Handler(BaseHTTPRequestHandler):
        def do_POST(self):
            query = parse_qs(self.rfile.read(int(self.headers['Content-Length'])).decode())['query'][0]
            received.append(query)
            bindings = []
            for key in sorted(set(re.findall(r'"([A-Z]{14}-[A-Z]{10}-[A-Z])"', query))):
                # URI は実際の PubChem RDF の形にする。生物種の末尾は "TAXID10090" で、
                # "http://example.org/10090" のような偽の形だと背景の "10090" との食い違いが隠れる。
                bindings.append(dict(value_inchikey={'value':key},
                                     pathway={'value':'http://rdf.ncbi.nlm.nih.gov/pubchem/pathway/PWID1'},
                                     pathway_label={'value':'Synthetic pathway'},
                                     pathway_organism={'value':'http://rdf.ncbi.nlm.nih.gov/pubchem/taxonomy/TAXID10090'}))
            raw = json.dumps({'results': {'bindings': bindings}}).encode()
            self.send_response(200)
            self.send_header('Content-Type', 'application/sparql-results+json')
            self.end_headers()
            self.wfile.write(raw)
        def log_message(self, *args):
            pass

    http = HTTPServer(('127.0.0.1', 0), Handler)
    thread = threading.Thread(target=http.serve_forever, daemon=True)
    thread.start()
    index = conn.execute('PRAGMA database_list').fetchone()[2]

    async def run():
        params = StdioServerParameters(command=sys.executable, args=['-m','massbank_context.server'],
                                       env={**os.environ, 'MASSBANK_CONTEXT_ARTIFACTS':str(tmp_path/'store'),
                                            'MASSBANK_CONTEXT_INDEX':index,
                                            'MASSBANK_CONTEXT_PUBCHEM_SPARQL':f'http://127.0.0.1:{http.server_port}'})
        async with stdio_client(params) as (read, write):
            async with ClientSession(read, write) as session:
                await session.initialize()
                async def call(name, **args):
                    result = await session.call_tool(name, args)
                    payload = json.loads(result.content[0].text)
                    assert not result.isError and 'error' not in payload, payload
                    return payload
                aid = (await call('load_evidence_bundle', path=str(path)))['analysis_id']
                await call('set_analysis_context', analysis_id=aid, confirmed=True,
                           context={'species': {'value': {'label':'synthetic mouse','taxid':'10090'},
                                                'source':'fixture researcher'}})
                matched = await call('match_feature_spectra', analysis_id=aid)
                assert matched['n_candidates'] > 0
                await call('resolve_candidate_evidence', analysis_id=aid)
                page = await call('get_comparison_evidence', analysis_id=aid, limit=1)
                ids = set()
                while True:
                    ids.update(f['feature_id'] for f in read_tsv(page['features']))
                    assert page['summary']['n_features_total'] == 8
                    if not page['next_cursor']:
                        break
                    page = await call('get_comparison_evidence', analysis_id=aid, limit=1, cursor=page['next_cursor'])
                assert len(ids) == 8
                result = await call('draft_comparison_findings', analysis_id=aid,
                                    organism='10090', include_pathways=True)
                assert result['readiness']['context_confirmed']
                evidence = read_tsv(result['facts']['evidence'])
                assert evidence
                assert result['hypotheses'][0]['claim'] is None
                assert result['pathways']['n_features_total'] == 8
                # 母集団に残るのは KEY_A の 1 化合物（f1・f2）で、マウスの PWID1 が付く。
                assert result['pathways']['n_with_species_pathways'] == 1, result['pathways']
                assert [r['pathway_id'] for r in result['pathways']['rows']] == ['PWID1']
                candidates = read_tsv(result['facts']['candidates'])
                row = evidence[0]
                claim = dict(claim_id='c1', feature_ids=[row['feature_id']],
                             candidate_ids=[row['candidate_id']], evidence_ids=[row['evidence_id']],
                             evidence_revision=result['evidence_revision'])
                checked = await call('validate_comparison_claims', analysis_id=aid, claims=[claim])
                assert checked['valid'], checked
                assert any(c['source'] == 'massbank' for c in candidates)

    try:
        asyncio.run(run())
    finally:
        http.shutdown()
        http.server_close()
        thread.join()
    assert received
    assert all('fixture researcher' not in q and 'synthetic-run' not in q for q in received)
