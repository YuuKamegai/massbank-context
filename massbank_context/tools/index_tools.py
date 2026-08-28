"""索引の状態を返すツール。"""
from __future__ import annotations

import sqlite3

from massbank_context.core import config
from massbank_context.core.errors import missing_index
from massbank_context.core.serialization import json_payload, round_floats
from massbank_context.index.writer import read_meta

__all__ = ["index_status"]


def index_status() -> str:
    """索引の有無と、その層・件数・ビン幅を返す。"""
    if not config.INDEX_PATH.exists():
        return missing_index(config.BUILD_COMMAND)
    conn = sqlite3.connect(str(config.INDEX_PATH))
    try:
        meta = read_meta(conn)
        if meta is None:
            return missing_index(config.BUILD_COMMAND)
        polarities = sorted(
            row[0] for row in conn.execute(
                "SELECT DISTINCT polarity FROM records WHERE polarity IS NOT NULL"))
    finally:
        conn.close()
    return json_payload(round_floats({**meta, "polarities": polarities}))
