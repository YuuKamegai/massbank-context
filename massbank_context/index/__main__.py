"""索引ビルド CLI（spec §7.1）。

TTL 全走査は分〜時間を要する。MCP ツールにすると LLM が会話の途中で起動して
長時間固まるため、CLI に出す。MCP 側は index_status() で状態を読むだけ。

使い方:
    python -m massbank_context.index build --ttl-dir <dir> --output <sqlite>

ビン幅は全レコードの質量誤差を測らないと決まらないため 2 パス構成にする。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from massbank_context.core.serialization import json_payload
from massbank_context.index.background import build_background
from massbank_context.index.mass_error import choose_bin_width, measure_mass_error
from massbank_context.index.writer import IndexWriter, open_index
from massbank_context.ttl.reader import read_ttl_file


def build_index(ttl_dir: str | Path, output: str | Path,
                resolution_filter: str = "high",
                limit: int | None = None) -> dict:
    """TTL ディレクトリを走査して索引を作る。"""
    ttl_dir = Path(ttl_dir)
    paths = sorted(ttl_dir.rglob("*.ttl"))
    if limit is not None:
        paths = paths[:limit]

    # --- パス 1: 全部読む ---
    records = []
    failed = 0
    for path in paths:
        try:
            records.extend(read_ttl_file(path))
        except Exception:           # 1 件の異常でビルド全体を止めない
            failed += 1

    # --- ビン幅を実測から決める（spec §6.4）---
    profile = measure_mass_error(records)
    bin_width = choose_bin_width(profile)

    # --- パス 2: ビン化して書く ---
    conn = open_index(output, create=True)
    writer = IndexWriter(conn, bin_width=bin_width,
                         resolution_filter=resolution_filter)
    accepted = 0
    for record in records:
        if writer.add_record(record):
            accepted += 1
    writer.write_mass_error_profile(profile)
    writer.finalize({"ttl_source": str(ttl_dir)})
    background_rows = build_background(conn)
    conn.close()

    return {
        "scanned": len(paths),
        "parsed": len(records),
        "failed": failed,
        "accepted": accepted,
        "rejected": len(records) - accepted,
        "bin_width": bin_width,
        "background_rows": background_rows,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m massbank_context.index")
    sub = parser.add_subparsers(dest="command", required=True)

    build = sub.add_parser("build", help="TTL から SQLite 索引を作る")
    build.add_argument("--ttl-dir", required=True)
    build.add_argument("--output", required=True)
    build.add_argument("--resolution-filter", default="high",
                       help="'high' / 'low' / '' (無効化)")
    build.add_argument("--limit", type=int, default=None,
                       help="先頭 N ファイルだけ処理する（動作確認用）")

    args = parser.parse_args(argv)
    stats = build_index(args.ttl_dir, args.output,
                        resolution_filter=args.resolution_filter,
                        limit=args.limit)
    print(json_payload(stats))
    return 0


if __name__ == "__main__":
    sys.exit(main())
