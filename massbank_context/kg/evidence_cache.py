"""エンドポイントとクエリ版で分離したRDF取得履歴。"""
import json
import sqlite3
from pathlib import Path
from massbank_context.core.evidence_ids import canonical_json, stable_id


class EvidenceCache:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with sqlite3.connect(self.path) as conn:
            conn.execute('CREATE TABLE IF NOT EXISTS rdf_snapshots '
                         '(id TEXT PRIMARY KEY, endpoint TEXT, version TEXT, inchikey TEXT, payload TEXT)')

    def lookup(self, endpoint, query_version, inchikey):
        with sqlite3.connect(self.path) as conn:
            row = conn.execute('SELECT payload FROM rdf_snapshots WHERE endpoint=? AND version=? '
                               'AND inchikey=? ORDER BY rowid DESC LIMIT 1',
                               (endpoint, query_version, inchikey)).fetchone()
        return json.loads(row[0]) if row else None

    def put(self, snapshot: dict):
        raw = canonical_json(snapshot)
        with sqlite3.connect(self.path) as conn:
            conn.execute('INSERT OR IGNORE INTO rdf_snapshots VALUES (?,?,?,?,?)',
                         (stable_id('lookup', [raw]), snapshot['endpoint'], snapshot['query_version'],
                          snapshot['inchikey'], raw))
