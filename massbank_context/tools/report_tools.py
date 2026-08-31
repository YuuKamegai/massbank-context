"""提示層のツール（spec §10）。"""
from __future__ import annotations

from massbank_context.core.serialization import json_payload, round_floats
from massbank_context.core.session import session
from massbank_context.report.findings import build_findings
from massbank_context.tools.pathway_tools import _activity_payload

__all__ = ["draft_pathway_findings"]


def draft_pathway_findings(organism: str = "", min_members: int = 3,
                           top_k: int = 10) -> str:
    """事実・仮説スロット・制約の骨組みを返す（spec §10）。

    pathway_activity と同じ計算を使い、結果を事実として写す。解釈は書かない。
    「活性化」「不活性化」の語は hypotheses の中にしか置けない。

    計算は `_activity_payload` を **pathway_activity と共有する**（計画 Task 8）。
    公開ツールの TSV を再パースして型を推測し直すことはしない — 往復すると
    型を推測で復元する羽目になり（"nan" が float('nan') に戻る）、
    pathway_label にタブが含まれれば列がずれて top_pathways が黙って壊れる。

    provenance の endpoint / fetched_at と population の pathway_coverage は
    payload にある値をそのまま写す（ここで再計算しない）。truncated も同様に
    写す（打ち切りを黙って落とさない）。
    """
    error, rows, payload = _activity_payload(organism, min_members, top_k)
    if error:
        return error

    differential = session.differential
    facts = {
        "comparison": {
            "group_a": payload["group_a"],
            "group_b": payload["group_b"],
            "n_a": differential.integer_meta("n_a"),
            "n_b": differential.integer_meta("n_b"),
            "log2fc_sign": payload["log2fc_sign"],
        },
        "population": {
            "n_features_total": differential.n_features_total,
            "n_with_inchikey": len(differential.inchikeys),
            "n_unannotated": differential.n_unannotated,
            "n_background": payload["n_background"],
            "n_no_pathway": payload["n_no_pathway"],
            "n_collapsed_spots": payload["n_collapsed_spots"],
            "n_pathways_total": payload["n_pathways_total"],
            "n_pathways_tested": payload["n_pathways_tested"],
            "pathway_coverage": payload["pathway_coverage"],
            "truncated": payload["truncated"],
        },
        "provenance": {
            "source": session.differential_source,
            "organism": payload["organism"],
            "endpoint": payload["endpoint"],
            "fetched_at": payload["fetched_at"],
        },
    }
    return json_payload(round_floats(
        build_findings(facts, rows, payload["pathway_coverage"])))
