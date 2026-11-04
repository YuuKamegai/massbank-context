"""MS/MS の筋: 索引の読み出し（このタスクでは open_index だけ）。"""
from __future__ import annotations

import sqlite3

from massbank_context.core import config
from massbank_context.index.writer import open_index_readonly, read_meta

__all__ = ["ACQUISITION_CLASS", "open_index"]

#: 外来の分子の走査が照合する層。測定原理の異なる層は跨がない（spec §6.1）。
ACQUISITION_CLASS = "LC-ESI-MS2"


def open_index() -> tuple[sqlite3.Connection, dict]:
    """索引を読み取り専用で開き、(接続, index_meta) を返す。

    索引が無い・読めない場合は ``missing_index``、LC-ESI-MS2 の層が無い場合は
    ``incompatible_index`` を ValueError で返す。別の層で黙って代用しない。
    """
    path = config.INDEX_PATH
    if not path.exists():
        raise ValueError(f"missing_index: 索引がありません: {path}。{config.BUILD_COMMAND}")
    conn = open_index_readonly(path)
    meta = read_meta(conn)
    if meta is None:
        conn.close()
        raise ValueError(f"missing_index: 索引が空または壊れています: {path}。{config.BUILD_COMMAND}")
    if ACQUISITION_CLASS not in meta.get("acquisition_classes", []):
        conn.close()
        raise ValueError(f"incompatible_index: 索引に {ACQUISITION_CLASS} の層がありません"
                         f"（acquisition_classes={meta.get('acquisition_classes')}）")
    return conn, meta
