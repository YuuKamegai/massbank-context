"""パスウェイのローカルキャッシュ（spec §6.2）。

lookup_log が要点である。「照会したがヒット 0」と「まだ照会していない」を
区別できないと被覆の分母が壊れる。pathway_hit に行が無いことだけでは
両者が同じに見えるため、照会した事実そのものを別テーブルに記録する。

日数ベースの TTL は持たない。PubChem のパスウェイは滅多に変わらない一方、
日付でキャッシュが勝手に切れると同じ解析が日によって違う結果を出す。

cache_meta には contract_version に加えて endpoint も残す（spec §6.2）。
最後に使ったエンドポイントで上書きする形でよい — 来歴として「このキャッシュは
どのエンドポイントに照会したものか」が分かればよく、endpoint ごとの履歴は
lookup_log の各行が持っている。
"""
from __future__ import annotations

import sqlite3
from datetime import datetime, timezone
from pathlib import Path

CACHE_CONTRACT_VERSION = "1"


class CacheVersionError(RuntimeError):
    """キャッシュファイルの版が対応版と違う。呼び出し側が contract_mismatch に変換する。"""


_SCHEMA = """
CREATE TABLE IF NOT EXISTS pathway_hit (
    inchikey TEXT NOT NULL,
    pathway_id TEXT NOT NULL,
    pathway_label TEXT,
    taxid TEXT,
    PRIMARY KEY (inchikey, pathway_id, taxid)
);
CREATE INDEX IF NOT EXISTS idx_pathway_hit_inchikey ON pathway_hit(inchikey);
CREATE TABLE IF NOT EXISTS lookup_log (
    inchikey TEXT NOT NULL,
    endpoint TEXT NOT NULL,
    fetched_at TEXT NOT NULL,
    status TEXT NOT NULL,
    PRIMARY KEY (inchikey, endpoint)
);
CREATE TABLE IF NOT EXISTS organism (
    taxid TEXT PRIMARY KEY,
    label TEXT
);
CREATE TABLE IF NOT EXISTS cache_meta (
    key TEXT PRIMARY KEY,
    value TEXT
);
"""


def open_cache(path) -> sqlite3.Connection:
    """キャッシュを開く（無ければ作る）。版が違えば開かずに落とす。

    `CREATE TABLE IF NOT EXISTS` は既存テーブルに列を足さないので、スキーマを
    v2 にしても既存の v1 ファイルは**黙って開いてしまう**。そのまま読むと、
    欠けた列が既定値の顔をして統計に入る。書き込んだ版を読み返して照合する。
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(path))
    conn.executescript(_SCHEMA)
    row = conn.execute(
        "SELECT value FROM cache_meta WHERE key = 'contract_version'").fetchone()
    if row is not None and row[0] != CACHE_CONTRACT_VERSION:
        conn.close()
        raise CacheVersionError(
            f"パスウェイキャッシュの版が違います: {row[0]}"
            f"（対応: {CACHE_CONTRACT_VERSION}）。"
            f"{path} を削除してから resolve_pathways を実行し、引き直してください。")
    conn.execute(
        "INSERT OR IGNORE INTO cache_meta(key, value) VALUES ('contract_version', ?)",
        (CACHE_CONTRACT_VERSION,))
    conn.commit()
    return conn


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _put_meta_endpoint(conn: sqlite3.Connection, endpoint: str) -> None:
    """cache_meta の 'endpoint' を最後に使ったものへ更新する（spec §6.2）。"""
    conn.execute(
        "INSERT OR REPLACE INTO cache_meta(key, value) VALUES ('endpoint', ?)",
        (endpoint,))


def pending_inchikeys(conn: sqlite3.Connection, inchikeys: list[str],
                      endpoint: str, refresh: bool = False) -> list[str]:
    """まだ照会が成功していない InChIKey を、入力の並びのまま返す。

    status が failed のものは pending に含める（次回の再試行対象）。
    """
    if refresh:
        return list(inchikeys)
    done = {row[0] for row in conn.execute(
        "SELECT inchikey FROM lookup_log WHERE endpoint = ? "
        "AND status IN ('ok', 'empty')", (endpoint,))}
    return [ik for ik in inchikeys if ik not in done]


def failed_inchikeys(conn: sqlite3.Connection, inchikeys: list[str],
                     endpoint: str) -> list[str]:
    """直近の照会が、**そのエンドポイントで**失敗している InChIKey。

    endpoint で絞るのは pending_inchikeys と揃えるため。絞らないと、
    MASSBANK_CONTEXT_PUBCHEM_SPARQL を切り替えたときに旧エンドポイントの
    failed 行が n_failed に混じり、被覆の分母（len - n_failed）が不当に縮む。
    """
    failed = {row[0] for row in conn.execute(
        "SELECT inchikey FROM lookup_log WHERE status = 'failed' "
        "AND endpoint = ?", (endpoint,))}
    return [ik for ik in inchikeys if ik in failed]


def record_result(conn: sqlite3.Connection, inchikey: str, endpoint: str,
                  hits: list[dict]) -> None:
    """照会結果を記録する。hits が空でも「照会した」ことを必ず残す。

    taxid が無いヒットは NULL ではなく空文字列で入れる。PRIMARY KEY は
    (inchikey, pathway_id, taxid) だが、SQLite は NULL 同士を相異とみなすため、
    NULL のままだと同じヒットが何度でも重複して入る。
    """
    conn.execute("DELETE FROM pathway_hit WHERE inchikey = ?", (inchikey,))
    conn.executemany(
        "INSERT OR REPLACE INTO pathway_hit"
        "(inchikey, pathway_id, pathway_label, taxid) VALUES (?, ?, ?, ?)",
        [(inchikey, h["pathway_id"], h.get("pathway_label"),
          h.get("taxid") or "")
         for h in hits])
    conn.execute(
        "INSERT OR REPLACE INTO lookup_log"
        "(inchikey, endpoint, fetched_at, status) VALUES (?, ?, ?, ?)",
        (inchikey, endpoint, _now(), "ok" if hits else "empty"))
    _put_meta_endpoint(conn, endpoint)
    conn.commit()


def record_failure(conn: sqlite3.Connection, inchikey: str,
                   endpoint: str) -> None:
    """照会に失敗したことを記録する（次回の再試行対象として残す）。"""
    conn.execute(
        "INSERT OR REPLACE INTO lookup_log"
        "(inchikey, endpoint, fetched_at, status) VALUES (?, ?, ?, 'failed')",
        (inchikey, endpoint, _now()))
    _put_meta_endpoint(conn, endpoint)
    conn.commit()


def memberships(conn: sqlite3.Connection, inchikeys: list[str],
                taxid: str | None = None) -> dict[str, list[dict]]:
    """InChIKey ごとのパスウェイ所属。taxid を渡すとその生物種に絞る。"""
    if not inchikeys:
        return {}
    placeholders = ",".join("?" * len(inchikeys))
    sql = ("SELECT inchikey, pathway_id, pathway_label, taxid FROM pathway_hit "
           f"WHERE inchikey IN ({placeholders})")
    params = list(inchikeys)
    if taxid:
        sql += " AND taxid = ?"
        params.append(taxid)
    sql += " ORDER BY inchikey, pathway_id"
    result: dict[str, list[dict]] = {}
    for row in conn.execute(sql, params):
        result.setdefault(row[0], []).append({
            "pathway_id": row[1], "pathway_label": row[2], "taxid": row[3]})
    return result


def organism_counts(conn: sqlite3.Connection,
                    inchikeys: list[str]) -> list[dict]:
    """生物種ごとの**化合物数**（パスウェイ行数ではない）。件数の降順。"""
    if not inchikeys:
        return []
    placeholders = ",".join("?" * len(inchikeys))
    rows = conn.execute(
        "SELECT h.taxid, o.label, COUNT(DISTINCT h.inchikey) "
        "FROM pathway_hit h LEFT JOIN organism o ON o.taxid = h.taxid "
        f"WHERE h.inchikey IN ({placeholders}) "
        "GROUP BY h.taxid, o.label ORDER BY 3 DESC, 1 ASC",
        list(inchikeys)).fetchall()
    return [{"taxid": r[0], "label": r[1] or "", "n_compounds": r[2]}
            for r in rows]


def put_organism_labels(conn: sqlite3.Connection,
                        labels: dict[str, str]) -> None:
    """taxid のラベルを記録する（best effort。引けなくても解析は続く）。"""
    conn.executemany(
        "INSERT OR REPLACE INTO organism(taxid, label) VALUES (?, ?)",
        list(labels.items()))
    conn.commit()


def fetched_at_range(conn: sqlite3.Connection,
                     inchikeys: list[str]) -> tuple[str | None, str | None]:
    """照会日時の最小と最大（来歴として payload に出す）。"""
    if not inchikeys:
        return None, None
    placeholders = ",".join("?" * len(inchikeys))
    row = conn.execute(
        "SELECT MIN(fetched_at), MAX(fetched_at) FROM lookup_log "
        f"WHERE inchikey IN ({placeholders})", list(inchikeys)).fetchone()
    return row[0], row[1]
