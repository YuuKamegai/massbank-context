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
        "n_a": result.meta.get("n_a"),
        "n_b": result.meta.get("n_b"),
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
    conn = cache.open_cache(config.PATHWAY_CACHE_PATH)
    try:
        pending = cache.pending_inchikeys(conn, inchikeys, endpoint,
                                          refresh=refresh)
        n_failed = 0
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
            if len(adapter.failed) >= len(pending):
                # 照会した全件が落ちた = 個々の InChIKey の問題ではなく上流が死んでいる。
                # 「0 件見つかった」と同じ形で返すと、被覆 0% を結果として読まれる。
                return upstream_unavailable(
                    endpoint=endpoint,
                    message=(f"照会した {len(pending)} 件すべてが失敗しました。"
                             "エンドポイントに到達できていない可能性があります。"
                             "復旧後に resolve_pathways を再実行してください。"))
            n_failed = len(cache.failed_inchikeys(conn, inchikeys))
            taxids = [o["taxid"] for o in cache.organism_counts(conn, inchikeys)
                      if o["taxid"] and not o["label"]]
            if taxids:
                cache.put_organism_labels(conn,
                                          adapter.fetch_organism_labels(taxids))
        else:
            n_failed = len(cache.failed_inchikeys(conn, inchikeys))

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
    """
    if session.differential is None:
        return missing_state(
            "differential", ["load_differential"],
            "差次的結果が読み込まれていません。先に load_differential を実行してください。")

    differential = session.differential
    inchikeys = differential.inchikeys
    endpoint = config.PUBCHEM_SPARQL_ENDPOINT
    conn = cache.open_cache(config.PATHWAY_CACHE_PATH)
    try:
        pending = cache.pending_inchikeys(conn, inchikeys, endpoint)
        if pending:
            return missing_state(
                "resolved_pathways", ["resolve_pathways"],
                f"照会が完了していない InChIKey が {len(pending)} 件あります。"
                "部分的な母集団で統計を回すと被覆が変わって濃縮の p 値が動くため、"
                "resolve_pathways を実行してから呼んでください。")
        available = cache.organism_counts(conn, inchikeys)
        if not available:
            return missing_state(
                "resolved_pathways", ["resolve_pathways"],
                "パスウェイが 1 件も引けていません。resolve_pathways を実行してください。")
        if not organism or organism not in {o["taxid"] for o in available}:
            return incompatible_organism(
                available=available,
                message=("organism を指定してください。生物種を跨いだパスウェイを"
                         "混ぜると解釈できません（既定で human に寄せることもしません）。"))
        found = cache.memberships(conn, inchikeys, taxid=organism)
        fetched_at = list(cache.fetched_at_range(conn, inchikeys))
    finally:
        conn.close()

    by_inchikey: dict[str, dict] = {}
    for row in differential.rows:
        inchikey = row["inchikey"]
        current = by_inchikey.get(inchikey)
        # 同じ InChIKey が複数スポットに割れることがある。|log2fc| の大きい方を採る。
        if current is None or abs(row["log2fc"] or 0.0) > abs(current["log2fc"] or 0.0):
            by_inchikey[inchikey] = {
                "inchikey": inchikey,
                "log2fc": row["log2fc"],
                "significant": row["significant"],
            }
    background = [c for ik, c in by_inchikey.items() if ik in found]
    n_no_pathway = len(by_inchikey) - len(background)
    denominator = len(background) + n_no_pathway
    pathway_coverage = (len(background) / denominator) if denominator else 0.0

    rows, summary = pathway_activity_rows(background, found, min_members)
    truncated = len(rows) > top_k
    shown = rows[:top_k]
    table = as_tsv(_ACTIVITY_COLUMNS,
                   [[round_floats(row.get(column)) for column in _ACTIVITY_COLUMNS]
                    for row in shown])
    return json_payload(round_floats({
        "organism": organism,
        "group_a": differential.group_a,
        "group_b": differential.group_b,
        "log2fc_sign": f"上昇 = {differential.group_b} が高い。",
        "min_members": min_members,
        "n_background": summary["n_background"],
        "n_no_pathway": n_no_pathway,
        "n_significant": summary["n_significant"],
        "n_pathways_total": summary["n_pathways_total"],
        "n_pathways_tested": summary["n_pathways_tested"],
        "pathway_coverage": pathway_coverage,
        "endpoint": endpoint,
        "fetched_at": fetched_at,
        "truncated": truncated,
        "candidates": table,
        "note": ("direction は検定ではなく記述です。fisher_q と wilcoxon_p は"
                 "別々の量であり、合成しないでください。"),
    }))
