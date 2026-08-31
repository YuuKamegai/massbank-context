"""提示層のツール（spec §10）。"""
from __future__ import annotations

import json

from massbank_context.core.serialization import json_payload, round_floats
from massbank_context.core.session import session
from massbank_context.report.findings import build_findings
from massbank_context.tools.pathway_tools import pathway_activity

__all__ = ["draft_pathway_findings"]


def draft_pathway_findings(organism: str = "", min_members: int = 3,
                           top_k: int = 10) -> str:
    """事実・仮説スロット・制約の骨組みを返す（spec §10）。

    pathway_activity と同じ計算を使い、結果を事実として写す。解釈は書かない。
    「活性化」「不活性化」の語は hypotheses の中にしか置けない。

    provenance の endpoint / fetched_at と population の pathway_coverage は
    pathway_activity の payload にある値をそのまま写す（ここで再計算しない。
    同じ計算を二重に実装しないため）。
    """
    raw = pathway_activity(organism=organism, min_members=min_members,
                           top_k=top_k)
    payload = json.loads(raw)
    if "error" in payload:
        return raw

    differential = session.differential
    columns = payload["candidates"].splitlines()[0].split("\t")
    rows = [dict(zip(columns, line.split("\t")))
            for line in payload["candidates"].splitlines()[1:]]

    facts = {
        "comparison": {
            "group_a": payload["group_a"],
            "group_b": payload["group_b"],
            "n_a": differential.meta.get("n_a"),
            "n_b": differential.meta.get("n_b"),
            "log2fc_sign": payload["log2fc_sign"],
        },
        "population": {
            "n_features_total": differential.n_features_total,
            "n_with_inchikey": len(differential.inchikeys),
            "n_unannotated": differential.n_unannotated,
            "n_background": payload["n_background"],
            "n_no_pathway": payload["n_no_pathway"],
            "n_pathways_total": payload["n_pathways_total"],
            "n_pathways_tested": payload["n_pathways_tested"],
            "pathway_coverage": payload["pathway_coverage"],
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
