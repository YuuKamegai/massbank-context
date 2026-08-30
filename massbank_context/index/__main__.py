"""索引ビルド CLI（spec §7.1）。

TTL 全走査は分〜時間を要する。MCP ツールにすると LLM が会話の途中で起動して
長時間固まるため、CLI に出す。MCP 側は index_status() で状態を読むだけ。

使い方:
    python -m massbank_context.index build --ttl-dir <dir> --output <sqlite>

ビン幅は全レコードの質量誤差を測らないと決まらないため 2 パス構成にする。
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from massbank_context.core.serialization import json_payload
from massbank_context.index.background import build_background
from massbank_context.index.classyfire import read_classyfire_tsv
from massbank_context.index.mass_error import choose_bin_width, measure_mass_error
from massbank_context.index.writer import IndexWriter, is_indexable, open_index
from massbank_context.ttl.reader import read_ttl_file_with_generation

logger = logging.getLogger(__name__)


def build_index(ttl_dir: str | Path, output: str | Path,
                resolution_filter: str = "high",
                limit: int | None = None,
                classyfire_tsv: str | Path | None = None) -> dict:
    """TTL ディレクトリを走査して索引を作る。"""
    ttl_dir = Path(ttl_dir)
    paths = sorted(ttl_dir.rglob("*.ttl"))
    if limit is not None:
        paths = paths[:limit]

    # --- パス 1: 全部読む ---
    records = []
    failed = 0
    # 世代ごとの内訳。世代B の TTL 実物は入手できておらず（spec §11 確認事項 1）、
    # 同梱 fixture は SPARQL のパターンから合成したもの。実物の SIO ノードの
    # 入れ子が想定と違えば parse_gen_b は「メタデータは揃っているがピーク 0 本」の
    # レコードを返し、accepted だけを見ていると成功に見えてしまう。
    # `gen_b: 12000 files, 0 peaks` と一目で分かる形で出す。
    generations: dict[str, dict[str, int]] = {}
    for path in paths:
        try:
            generation, parsed = read_ttl_file_with_generation(path)
        except Exception as exc:    # 1 件の異常でビルド全体を止めない
            # どのファイルがなぜ失敗したかを残す。54,183 ファイルの走査で
            # 件数だけ報告されても運用者は追跡できない。
            logger.warning("TTL の読み込みに失敗しました: %s (%r)", path, exc)
            failed += 1
            continue
        counts = generations.setdefault(
            generation, {"files": 0, "records": 0, "peaks": 0})
        counts["files"] += 1
        counts["records"] += len(parsed)
        counts["peaks"] += sum(len(record.peaks) for record in parsed)
        records.extend(parsed)

    for generation, counts in sorted(generations.items()):
        logger.info("TTL 世代 %s: %d files, %d records, %d peaks", generation,
                    counts["files"], counts["records"], counts["peaks"])

    # --- ビン幅を実測から決める（spec §6.4）---
    # 測るのは索引に入るレコードだけ。choose_bin_width は装置クラス横断の
    # 最大値を採るため、ここに層の条件で弾かれるレコードを混ぜると、
    # 索引に 1 件も入らない装置クラスが索引の照合窓を決めてしまう。
    indexable = [r for r in records if is_indexable(r, resolution_filter)]
    logger.info("質量誤差の母集団: %d / %d レコード（層の条件を通ったもの）",
                len(indexable), len(records))
    profile = measure_mass_error(indexable)
    bin_width = choose_bin_width(profile)
    # 「ビン幅は決め打ちではなく実測」という主張は、この 1 関数の出力だけに
    # 乗っている。1 装置クラスの汚染が最大値経由で索引全体のビンを太らせるため、
    # 採用値と内訳の両方を残して運用者が妥当性を目視できるようにする。
    logger.info("採用ビン幅: %.6f（実測できた装置クラス %d 件）",
                bin_width, len(profile))
    for instrument_class, entries in sorted(profile.items()):
        logger.info("質量誤差プロファイル %s: %s", instrument_class,
                    " ".join(f"p{pct:g}={delta:.6f}" for pct, delta in entries))

    # --- パス 2: ビン化して書く ---
    conn = open_index(output, create=True)
    writer = IndexWriter(conn, bin_width=bin_width,
                         resolution_filter=resolution_filter)
    accepted = 0
    for record in records:
        if writer.add_record(record):
            accepted += 1
    writer.write_mass_error_profile(profile)

    # ClassyFire の取り込みはレコードを書いたあとに行う。add_compound_classes は
    # 索引に入った InChIKey だけを通すため、順序が逆だと全件が弾かれる。
    classyfire_stats = None
    if classyfire_tsv is not None:
        classyfire_stats = writer.add_compound_classes(
            read_classyfire_tsv(classyfire_tsv))
        logger.info("ClassyFire: %d 行を書き、索引外の InChIKey を %d 行弾きました",
                    classyfire_stats["written"],
                    classyfire_stats["skipped_unknown_inchikey"])
        if classyfire_stats["written"] == 0:
            logger.warning(
                "ClassyFire のクラスを 1 行も書けませんでした。TSV の InChIKey 表記が"
                "索引側と食い違っている可能性があります。")

    meta = {"ttl_source": str(ttl_dir)}
    if classyfire_stats is not None:
        # どの索引が ClassyFire を持っているかは層別背景（spec §6.6）の
        # 再現に必要な provenance であり、DB を見ただけでは
        # 「取り込んで 0 件だった」と「取り込んでいない」が区別できない。
        meta["classyfire_source"] = str(classyfire_tsv)
        meta["classyfire_rows"] = str(classyfire_stats["written"])
    if generations:
        # spec §6.7 が index_meta に挙げているのに書かれていなかった項目。
        # 世代B の対応は本プロジェクト唯一の未確認事項なので、「この索引は
        # どの世代から作られたか」は最も価値のある provenance である。
        # 複数世代が混ざったビルドでは全部を記録する。
        meta["ttl_generation"] = ",".join(sorted(generations))
    if limit is not None:
        # --limit で作った索引は部分コーパスなので、その事実を残す。
        # 残さないと全コーパス版と見分けがつかなくなる。
        meta["limit"] = str(limit)
    writer.finalize(meta)
    background_rows = build_background(conn)
    conn.close()

    if writer.n_no_peaks:
        logger.warning(
            "ピークが 1 本も無いレコードを %d 件弾きました。パーサが"
            "ピークを拾えていない可能性があります（世代ごとの内訳を確認）。",
            writer.n_no_peaks)

    return {
        "scanned": len(paths),
        "parsed": len(records),
        "failed": failed,
        "accepted": accepted,
        "rejected": len(records) - accepted,
        "no_peaks": writer.n_no_peaks,
        "generations": generations,
        "bin_width": bin_width,
        "background_rows": background_rows,
        "classyfire": classyfire_stats,
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
    build.add_argument("--classyfire-tsv", default=None,
                       help="ClassyFire クラスの TSV（MassBank-RDF の"
                            " script/massbank_inchikey_classyfire.tsv）。"
                            "指定すると compound_class に source 付きで取り込む")

    args = parser.parse_args(argv)
    # ビン幅と世代ごとの内訳は INFO で出す。CLI が唯一の実行経路なので、
    # ここでハンドラを付けないと運用者の目には一切触れない。
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")
    stats = build_index(args.ttl_dir, args.output,
                        resolution_filter=args.resolution_filter,
                        limit=args.limit,
                        classyfire_tsv=args.classyfire_tsv)
    print(json_payload(stats))
    return 0


if __name__ == "__main__":
    sys.exit(main())
