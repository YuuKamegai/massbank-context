"""パスウェイ系のツール。Phase 2〜3（spec §8）。"""
from __future__ import annotations

from pathlib import Path

from massbank_context.core import config
from massbank_context.core.errors import (
    contract_mismatch,
    incompatible_organism,
    missing_state,
    upstream_unavailable,
)
from massbank_context.core.serialization import as_tsv, json_payload, round_floats
from massbank_context.core.session import session
from massbank_context.experiment.contract import ContractError, parse_differential
from massbank_context.kg import cache
from massbank_context.kg.pathway import PubChemPathwayAdapter
from massbank_context.stats.activity import pathway_activity_rows

__all__ = ["load_differential", "resolve_pathways", "pathway_activity"]

#: 併合規則の説明。payload に毎回載せる（黙って行を減らさないため）。
_COLLAPSE_RULE = (
    "同一 InChIKey の複数スポットは 1 件に併合しています。"
    "採用するのは significant=true の行を優先し、その中で |log2fc| が"
    "最大のものです。"
)

_SIGN_NOTE = "log2fc は正なら group_b が高い（上昇）。massbank-context 側で反転しない。"
_BATCH_SIZE = 50

_ACTIVITY_COLUMNS = [
    "pathway_id", "pathway_label", "taxid", "n_members", "n_significant",
    "n_up", "n_down", "median_log2fc", "direction", "fisher_p", "fisher_q",
    "tested", "wilcoxon_p", "wilcoxon_tested",
]


def load_differential(path: str) -> str:
    """ms-data-parser が書き出した差次的結果を読み込む（spec §6.1）。

    行そのものはセッションに保持し、payload には件数と来歴だけを返す。
    n_unannotated は「注釈が付かず書き出されなかった行」であり、
    **「変化が無かった」ではなく「調べていない」**行である。
    """
    file_path = Path(path)
    if not file_path.is_file():
        return contract_mismatch(
            reason=f"ファイルが見つかりません: {path}",
            expected="ms-data-parser の arf_export_differential が書き出した TSV",
            path=path,
        )
    try:
        result = parse_differential(file_path.read_text(encoding="utf-8"))
    except ContractError as error:
        return contract_mismatch(
            reason=str(error),
            expected="contract_version = 1 / log2fc_sign = "
            "positive means group_b is higher / 必須列一式",
            path=path,
        )

    session.differential = result
    session.differential_source = str(file_path)
    return json_payload({
        "source": str(file_path),
        "group_a": result.group_a,
        "group_b": result.group_b,
        "n_a": result.integer_meta("n_a"),
        "n_b": result.integer_meta("n_b"),
        "n_rows": len(result.rows),
        "n_inchikeys": len(result.inchikeys),
        "n_features_total": result.n_features_total,
        "n_unannotated": result.n_unannotated,
        "n_significant": sum(1 for row in result.rows if row["significant"]),
        "log2fc_sign": _SIGN_NOTE,
        "note": (
            "n_unannotated は注釈が付かず書き出されなかった行数です。"
            "「変化が無かった」ではなく「調べていない」行です。"
        ),
    })


def _make_adapter():
    """アダプタの生成点を 1 か所にする（テストで差し替えられるようにするため）。"""
    return PubChemPathwayAdapter(config.PUBCHEM_SPARQL_ENDPOINT)


def _open_cache():
    """キャッシュを開く。版が違えば (None, エラー本文) を返す。

    版違いは呼び出し側のツールでは復旧できない（ファイルを消して引き直す）ので、
    どうすればよいかを contract_mismatch の reason に書いて返す。
    """
    try:
        return cache.open_cache(config.PATHWAY_CACHE_PATH), None
    except cache.CacheVersionError as error:
        return None, contract_mismatch(
            reason=str(error),
            expected=("パスウェイキャッシュの cache_meta.contract_version = "
                      f"{cache.CACHE_CONTRACT_VERSION}"),
            path=str(config.PATHWAY_CACHE_PATH))


def resolve_pathways(refresh: bool = False) -> str:
    """読み込んだ差次的結果の InChIKey にパスウェイを引く（spec §6.2）。

    キャッシュに無いものだけを 50 件ずつ照会する。refresh=True で全件引き直す
    （日数ベースの TTL は持たない。日付でキャッシュが切れると同じ解析が
    日によって違う結果を出すため）。

    照会に失敗した InChIKey が残っている間は pathway_activity を走らせない。
    部分的な母集団で統計を回すと被覆が見かけ上変わり、濃縮の p 値が動く。
    """
    if session.differential is None:
        return missing_state(
            "differential", ["load_differential"],
            "差次的結果が読み込まれていません。先に load_differential を実行してください。")

    inchikeys = session.differential.inchikeys
    endpoint = config.PUBCHEM_SPARQL_ENDPOINT
    conn, error = _open_cache()
    if error:
        return error
    try:
        pending = cache.pending_inchikeys(conn, inchikeys, endpoint,
                                          refresh=refresh)
        adapter = None
        if pending:
            adapter = _make_adapter()
            for start in range(0, len(pending), _BATCH_SIZE):
                batch = pending[start:start + _BATCH_SIZE]
                hits = adapter.fetch(batch)
                failed = set(adapter.failed)
                for inchikey in batch:
                    if inchikey in failed:
                        cache.record_failure(conn, inchikey, endpoint)
                    else:
                        cache.record_result(conn, inchikey, endpoint,
                                            hits.get(inchikey, []))
                n_processed = start + len(batch)
                if len(adapter.failed) >= n_processed:
                    # ここまでの全件が落ちた = 個々の InChIKey の問題ではなく上流が
                    # 死んでいる。残りのバッチを回しても同じで、1 件あたり
                    # DEFAULT_TIMEOUT 秒待つので、900 件なら 30 時間ブロックしてから
                    # エラーを返すことになる（ツール呼び出しとしては事実上ハング）。
                    # 「0 件見つかった」と同じ形で返すと、被覆 0% を結果として読まれる。
                    return upstream_unavailable(
                        endpoint=endpoint,
                        message=(f"照会した {n_processed} 件すべてが失敗したため、"
                                 f"残り {len(pending) - n_processed} 件を照会せずに"
                                 "打ち切りました。エンドポイントに到達できていない"
                                 "可能性があります。復旧後に resolve_pathways を"
                                 "再実行してください。"),
                        n_processed=n_processed,
                        n_remaining=len(pending) - n_processed)

        n_failed = len(cache.failed_inchikeys(conn, inchikeys, endpoint))
        # ラベル解決は pending の有無に依らず走らせる。全件キャッシュヒットの回に
        # 走らないと、未ラベルの taxid が refresh=True 以外では**永久に**埋まらない。
        # 埋まらない先は incompatible_organism.available の label であり、
        # 設計書 §6.3 が「候補を示して呼び出し側に選ばせる」ための唯一の手掛かりである。
        taxids = [o["taxid"] for o in cache.organism_counts(conn, inchikeys)
                  if o["taxid"] and not o["label"]]
        if taxids:
            if adapter is None:
                adapter = _make_adapter()
            cache.put_organism_labels(conn,
                                      adapter.fetch_organism_labels(taxids))

        found = cache.memberships(conn, inchikeys)
        organisms = cache.organism_counts(conn, inchikeys)
        first, last = cache.fetched_at_range(conn, inchikeys)
    finally:
        conn.close()

    n_with_pathway = len(found)
    n_resolved = len(inchikeys) - n_failed
    coverage = (n_with_pathway / n_resolved) if n_resolved else 0.0
    note = ("パスウェイが付かなかった化合物は「変化が無かった」のではなく"
            "「調べていない」ものです。")
    if n_failed:
        note += (f" 照会に失敗した InChIKey が {n_failed} 件あります。"
                 "resolve_pathways を再実行してください"
                 "（未完了のまま pathway_activity は走りません）。")
    return json_payload(round_floats({
        "endpoint": endpoint,
        "n_inchikeys": len(inchikeys),
        "n_queried": len(pending),
        "n_cached": len(inchikeys) - len(pending),
        "n_with_pathway": n_with_pathway,
        "n_failed": n_failed,
        "pathway_coverage": coverage,
        "organisms": organisms,
        "fetched_at": [first, last],
        "note": note,
    }))


def _collapse_by_inchikey(rows: list[dict]) -> dict[str, dict]:
    """同一 InChIKey の複数スポットを 1 件に併合する。

    リピドミクスでは同じ InChIKey が異性体・アダクト・in-source フラグメントで
    複数スポットに割れるのが常態である。パスウェイ単位の Fisher は化合物を
    単位に数えるので、ここで 1 件に落とす必要がある。

    採用順位は **significant=true を優先 → その中で |log2fc| 最大**。
    |log2fc| だけで選ぶと、有意な行が非有意な行に置き換わり得る
    （例: +1.2/significant と -3.0/非 significant なら後者が勝つ）。
    そうなると n_significant が減り、Fisher の k と n が動いて濃縮の p 値まで
    変わる。significant は上流 ms-data-parser の閾値による判定であり、
    こちら側で取り消してはならない（spec §7.2）。
    """
    def _rank(row: dict) -> tuple[bool, float]:
        return (bool(row["significant"]), abs(row["log2fc"] or 0.0))

    collapsed: dict[str, dict] = {}
    for row in rows:
        inchikey = row["inchikey"]
        current = collapsed.get(inchikey)
        if current is None or _rank(row) > _rank(current):
            collapsed[inchikey] = {
                "inchikey": inchikey,
                "log2fc": row["log2fc"],
                "significant": row["significant"],
            }
    return collapsed


def _activity_payload(organism: str, min_members: int, top_k: int):
    """pathway_activity の中身。`(error, rows, payload)` を返す。

    rows は **typed dict のまま**の上位 top_k 行で、payload には candidates
    （TSV）が入っていない。TSV 化は公開ツールの出口だけで行う。
    report_tools はここから rows を直接受け取る — TSV に落として再パースすると、
    型の推測をやり直すことになる（as_tsv はタブ等をエスケープするが、
    再パース側で戻す手間と型の喪失は残る）。

    rows は _ACTIVITY_COLUMNS へ射影する。wilcoxon_q は BH の内部量であり、
    fisher_q と並べると合成したくなる（spec §5.3）ので外へは出さない。
    """
    if session.differential is None:
        return missing_state(
            "differential", ["load_differential"],
            "差次的結果が読み込まれていません。先に load_differential を実行してください。"
        ), [], {}

    differential = session.differential
    inchikeys = differential.inchikeys
    endpoint = config.PUBCHEM_SPARQL_ENDPOINT
    conn, error = _open_cache()
    if error:
        return error, [], {}
    try:
        pending = cache.pending_inchikeys(conn, inchikeys, endpoint)
        if pending:
            return missing_state(
                "resolved_pathways", ["resolve_pathways"],
                f"照会が完了していない InChIKey が {len(pending)} 件あります。"
                "部分的な母集団で統計を回すと被覆が変わって濃縮の p 値が動くため、"
                "resolve_pathways を実行してから呼んでください。"), [], {}
        available = cache.organism_counts(conn, inchikeys)
        if not available:
            return missing_state(
                "resolved_pathways", ["resolve_pathways"],
                "パスウェイが 1 件も引けていません。resolve_pathways を実行してください。"
            ), [], {}
        if not organism or organism not in {o["taxid"] for o in available}:
            return incompatible_organism(
                available=available,
                message=("organism を指定してください。生物種を跨いだパスウェイを"
                         "混ぜると解釈できません（既定で human に寄せることも"
                         "しません）。")), [], {}
        found = cache.memberships(conn, inchikeys, taxid=organism)
        fetched_at = list(cache.fetched_at_range(conn, inchikeys))
    finally:
        conn.close()

    by_inchikey = _collapse_by_inchikey(differential.rows)
    background = [c for ik, c in by_inchikey.items() if ik in found]
    n_no_pathway = len(by_inchikey) - len(background)
    # 分母は設計書 §7.3 の「照会が完了した InChIKey 数」である。ここでは
    # len(background) + n_no_pathway = len(by_inchikey) と書いているが、両者が
    # 一致するのは**上の pending ガードに依存した不変条件**である
    # （pending が 1 件でも残っていればこの関数はここまで来ない）。
    # ガードを緩めるなら、この分母も lookup_log から数え直すこと。
    denominator = len(background) + n_no_pathway
    pathway_coverage = (len(background) / denominator) if denominator else 0.0

    rows, summary = pathway_activity_rows(background, found, min_members)
    truncated = len(rows) > top_k
    shown = [{column: row.get(column) for column in _ACTIVITY_COLUMNS}
             for row in rows[:top_k]]
    payload = {
        "organism": organism,
        "group_a": differential.group_a,
        "group_b": differential.group_b,
        "log2fc_sign": f"上昇 = {differential.group_b} が高い。",
        "min_members": min_members,
        "n_rows": len(differential.rows),
        "n_inchikeys": len(by_inchikey),
        "n_collapsed_spots": len(differential.rows) - len(by_inchikey),
        "n_background": summary["n_background"],
        "n_no_pathway": n_no_pathway,
        "n_significant": summary["n_significant"],
        "n_pathways_total": summary["n_pathways_total"],
        "n_pathways_tested": summary["n_pathways_tested"],
        "pathway_coverage": pathway_coverage,
        "endpoint": endpoint,
        "fetched_at": fetched_at,
        "truncated": truncated,
        "collapse_rule": _COLLAPSE_RULE,
        "note": ("direction は検定ではなく記述です。fisher_q と wilcoxon_p は"
                 "別々の量であり、合成しないでください。"),
    }
    return None, shown, payload


def pathway_activity(organism: str = "", min_members: int = 3,
                     top_k: int = 30) -> str:
    """パスウェイ単位の濃縮と方向を返す（spec §7）。

    organism は必須である。生物種を跨いだパスウェイを黙って混ぜず、既定で
    human に寄せることもしない（実測で非哺乳類が 4 割混ざる）。

    背景は「InChIKey が付き、照会が完了し、指定 organism のパスウェイを
    1 つ以上持つ化合物」に限る。パスウェイを持たない化合物を背景に含めると
    全パスウェイが一律に有意寄りになる。除外件数は n_no_pathway で返す。

    並び順は濃縮の q 値だけで決まる。方向の p 値では並べ替えない
    （スコアを合成しない・spec §5.3）。

    facts（spec §10）: pathway_coverage は「パスウェイが付いた化合物の割合」の
    実測値であり、endpoint / fetched_at は提示層が来歴として出すための値。
    fetched_at は conn を閉じる前に取得する。

    同一 InChIKey が複数スポットに割れている場合は 1 件に併合する。
    併合が起きた件数は n_collapsed_spots として必ず出す（spec §12）。
    """
    error, rows, payload = _activity_payload(organism, min_members, top_k)
    if error:
        return error
    payload["candidates"] = as_tsv(
        _ACTIVITY_COLUMNS,
        [[round_floats(row.get(column)) for column in _ACTIVITY_COLUMNS]
         for row in rows])
    return json_payload(round_floats(payload))
