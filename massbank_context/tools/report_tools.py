"""提示層のツール（spec §10）。"""
from __future__ import annotations

import json

from massbank_context.core.serialization import json_payload, round_floats
from massbank_context.core.session import session
from massbank_context.report.findings import build_findings
from massbank_context.tools.pathway_tools import pathway_activity

__all__ = ["draft_pathway_findings"]

#: TSV に戻しても文字列のまま残す列（コントローラ裁定）。
_STRING_COLUMNS = {"pathway_id", "pathway_label", "taxid", "direction"}


def _coerce_cell(value: str):
    """TSV セルの文字列を元の型に復元する（コントローラ裁定）。

    as_tsv() は全セルを str() に通すため、そのまま dict に詰めると
    facts.top_pathways が全部文字列になり、設計書 §10 の「数値のみ」に反する。
    空文字列は None、"True"/"False" は bool、それ以外は int → float の順に
    変換を試み、どちらも失敗すれば文字列のまま返す。
    """
    if value == "":
        return None
    if value == "True":
        return True
    if value == "False":
        return False
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    return value


def _rows_from_tsv(table: str) -> list[dict]:
    """as_tsv() で作った TSV を行の dict のリストに戻す。数値・bool 列は型を復元する。"""
    lines = table.splitlines()
    columns = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        cells = line.split("\t")
        row = {}
        for column, cell in zip(columns, cells):
            row[column] = cell if column in _STRING_COLUMNS else _coerce_cell(cell)
        rows.append(row)
    return rows


def draft_pathway_findings(organism: str = "", min_members: int = 3,
                           top_k: int = 10) -> str:
    """事実・仮説スロット・制約の骨組みを返す（spec §10）。

    pathway_activity と同じ計算を使い、結果を事実として写す。解釈は書かない。
    「活性化」「不活性化」の語は hypotheses の中にしか置けない。

    provenance の endpoint / fetched_at と population の pathway_coverage は
    pathway_activity の payload にある値をそのまま写す（ここで再計算しない。
    同じ計算を二重に実装しないため）。truncated も同様に写す（打ち切りを
    黙って落とさない）。
    """
    raw = pathway_activity(organism=organism, min_members=min_members,
                           top_k=top_k)
    payload = json.loads(raw)
    if "error" in payload:
        return raw

    differential = session.differential
    rows = _rows_from_tsv(payload["candidates"])

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
