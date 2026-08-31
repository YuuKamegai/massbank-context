"""パスウェイ系のツール。Phase 2〜3（spec §8）。"""
from __future__ import annotations

from pathlib import Path

from massbank_context.core import config
from massbank_context.core.errors import (
    contract_mismatch,
    missing_state,
    upstream_unavailable,
)
from massbank_context.core.serialization import json_payload, round_floats
from massbank_context.core.session import session
from massbank_context.experiment.contract import ContractError, parse_differential
from massbank_context.kg import cache
from massbank_context.kg.pathway import PubChemPathwayAdapter

__all__ = ["load_differential", "resolve_pathways"]

_SIGN_NOTE = "log2fc は正なら group_b が高い（上昇）。massbank-context 側で反転しない。"
_BATCH_SIZE = 50


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
