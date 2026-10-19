"""解析入力と不変の根拠版をSQLiteへ保存する。"""
import json
import sqlite3
from dataclasses import asdict
from pathlib import Path
from massbank_context.core.evidence_ids import canonical_json, stable_id
from massbank_context.core.evidence_models import Candidate, EvidenceBundle, Feature, SpectrumLink
from massbank_context.core.models import Peak, QuerySpectrum


class EvidenceStore:
    def __init__(self, root: Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.path = self.root / 'evidence.sqlite'
        with self.connect() as conn:
            conn.executescript('''
                CREATE TABLE IF NOT EXISTS analyses (id TEXT PRIMARY KEY, payload TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS snapshots (
                    id TEXT PRIMARY KEY, analysis_id TEXT NOT NULL, kind TEXT NOT NULL, payload TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS heads (
                    analysis_id TEXT NOT NULL, kind TEXT NOT NULL, revision TEXT NOT NULL,
                    PRIMARY KEY(analysis_id, kind));
                CREATE TABLE IF NOT EXISTS sources (
                    analysis_id TEXT NOT NULL, path TEXT NOT NULL, payload BLOB NOT NULL,
                    PRIMARY KEY(analysis_id, path));
            ''')

    def connect(self):
        return sqlite3.connect(self.path, timeout=30)

    def register(self, bundle: EvidenceBundle) -> str:
        data = asdict(bundle)
        data['manifest_path'] = str(bundle.manifest_path)
        with self.connect() as conn:
            conn.execute('INSERT OR IGNORE INTO analyses VALUES (?,?)',
                         (bundle.analysis_id, canonical_json(data)))
            for file in bundle.manifest.get('files', []):
                raw = (bundle.manifest_path.parent / file['path']).read_bytes()
                import hashlib
                if hashlib.sha256(raw).hexdigest() != file['sha256']:
                    raise ValueError('登録中に入力が変更されました')
                conn.execute('INSERT OR IGNORE INTO sources VALUES (?,?,?)',
                             (bundle.analysis_id, file['path'], raw))
        return bundle.analysis_id

    def load(self, analysis_id: str) -> EvidenceBundle:
        with self.connect() as conn:
            row = conn.execute('SELECT payload FROM analyses WHERE id=?', (analysis_id,)).fetchone()
        if row is None:
            raise ValueError('解析IDが見つかりません')
        data = json.loads(row[0])
        spectra = {}
        for key, spec in data['spectra'].items():
            spec['peaks'] = tuple(Peak(**p) for p in spec['peaks'])
            spectra[key] = QuerySpectrum(**spec)
        return EvidenceBundle(data['analysis_id'], data['manifest'],
                              tuple(Feature(**f) for f in data['features']),
                              tuple(Candidate(**c) for c in data['annotations']),
                              tuple(SpectrumLink(**s) for s in data['spectrum_links']),
                              spectra, data['population_complete'], Path(data['manifest_path']))

    def put_snapshot(self, analysis_id: str, kind: str, payload: dict) -> str:
        self.load(analysis_id)
        raw = canonical_json(payload)
        revision = stable_id('snapshot', [analysis_id, kind, raw])
        with self.connect() as conn:
            conn.execute('INSERT OR IGNORE INTO snapshots VALUES (?,?,?,?)',
                         (revision, analysis_id, kind, raw))
            conn.execute('INSERT OR REPLACE INTO heads VALUES (?,?,?)',
                         (analysis_id, kind, revision))
        return revision

    def get_snapshot(self, analysis_id: str, revision: str) -> dict:
        with self.connect() as conn:
            row = conn.execute('SELECT payload FROM snapshots WHERE id=? AND analysis_id=?',
                               (revision, analysis_id)).fetchone()
        if row is None:
            raise ValueError('根拠版が見つかりません')
        return json.loads(row[0])

    def current(self, analysis_id: str, kind: str) -> dict:
        with self.connect() as conn:
            row = conn.execute('SELECT revision FROM heads WHERE analysis_id=? AND kind=?',
                               (analysis_id, kind)).fetchone()
        return self.get_snapshot(analysis_id, row[0]) if row else {}

    def latest_revision(self, analysis_id: str) -> str:
        self.load(analysis_id)
        with self.connect() as conn:
            rows = conn.execute('SELECT kind, revision FROM heads WHERE analysis_id=? ORDER BY kind',
                                (analysis_id,)).fetchall()
        return stable_id('revision', [analysis_id, canonical_json(rows)])

    def read_view(self, analysis_id: str) -> tuple[EvidenceBundle, str, dict]:
        """単一トランザクションの読取でページの版を固定する。"""
        with self.connect() as conn:
            conn.execute('BEGIN')
            row = conn.execute('SELECT payload FROM analyses WHERE id=?', (analysis_id,)).fetchone()
            if row is None:
                raise ValueError('解析IDが見つかりません')
            rows = conn.execute('SELECT h.kind, h.revision, s.payload FROM heads h JOIN snapshots s '
                                'ON s.id=h.revision WHERE h.analysis_id=? ORDER BY h.kind', (analysis_id,)).fetchall()
        revision = stable_id('revision', [analysis_id, canonical_json([(r[0], r[1]) for r in rows])])
        return self.load(analysis_id), revision, {r[0]: json.loads(r[2]) for r in rows}
