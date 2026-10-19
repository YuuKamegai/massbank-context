"""索引に実在するメタデータだけを返す。"""


def read_candidate_catalog(conn, *, accessions: list[str], inchikeys: list[str]) -> list[dict]:
    result = []
    old_factory = conn.row_factory
    try:
        import sqlite3
        conn.row_factory = sqlite3.Row
        for accession in sorted(set(accessions)):
            row = conn.execute('SELECT * FROM records WHERE accession=?', (accession,)).fetchone()
            if row:
                result.append(dict(row))
        keys = set(inchikeys) | {r['inchikey'] for r in result}
        compounds = {}
        for ik in sorted(keys):
            row = conn.execute('SELECT * FROM compounds WHERE inchikey=?', (ik,)).fetchone()
            classes = [dict(r) for r in conn.execute('SELECT source,label FROM compound_class '
                                                    'WHERE inchikey=? ORDER BY source,label', (ik,))]
            compounds[ik] = dict(compound=dict(row) if row else None, classes=classes)
        for row in result:
            row.update(compounds[row['inchikey']], ttl_source=None)
        for ik in sorted(set(inchikeys)):
            result.append(dict(inchikey=ik, accession=None, ttl_source=None, **compounds[ik]))
        return result
    finally:
        conn.row_factory = old_factory
