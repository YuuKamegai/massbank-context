"""パスウェイ単位の集計（純関数）。

plain な dict と list だけを受け取り、plain な dict を返す。experiment/ や kg/ の
データモデルを知らないので、stats/ は leaf のままでいられる。

濃縮と方向は**別々の量**であり、合成しない（spec §5.3）。並び順を決めるのは
濃縮の q 値だけで、方向の p 値では並べ替えない。
"""
from __future__ import annotations

from massbank_context.stats.fdr import benjamini_hochberg
from massbank_context.stats.fisher import fisher_exact_greater
from massbank_context.stats.wilcoxon import wilcoxon_signed_rank


def _median(values: list[float]) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2 == 1:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def _direction(n_up: int, n_down: int) -> str:
    """有意メンバーの上昇・低下の内訳から方向ラベルを決める。

    これは検定ではなく**記述**である。恣意的な閾値を置かないため、
    「片方が 0 かどうか」だけで判定する。
    """
    n_significant = n_up + n_down
    if n_significant == 0:
        return "NO_SIGNIFICANT_MEMBER"
    if n_significant == 1:
        return "SINGLE_MEMBER"
    if n_down == 0:
        return "UP_DOMINANT"
    if n_up == 0:
        return "DOWN_DOMINANT"
    return "MIXED"


def pathway_activity_rows(compounds: list[dict],
                          memberships: dict[str, list[dict]],
                          min_members: int) -> tuple[list[dict], dict]:
    """背景の化合物とパスウェイ所属から、パスウェイ単位の行を作る。

    引数:
        compounds: [{"inchikey", "log2fc", "significant"}]。**背景の化合物だけ**を渡す
            （パスウェイを 1 つも持たない化合物は呼び出し側で除外しておく）。
        memberships: inchikey -> [{"pathway_id", "pathway_label", "taxid"}]。
        min_members: この数未満のメンバーしか持たないパスウェイは検定しない。

    戻り値: (rows, summary)。rows は濃縮の q 値昇順、同値は pathway_id 昇順。
    """
    by_inchikey = {compound["inchikey"]: compound for compound in compounds}
    n_background = len(by_inchikey)
    n_significant_total = sum(
        1 for compound in by_inchikey.values() if compound.get("significant"))

    grouped: dict[str, dict] = {}
    for inchikey, entries in memberships.items():
        if inchikey not in by_inchikey:
            continue
        for entry in entries:
            pathway_id = entry["pathway_id"]
            bucket = grouped.setdefault(pathway_id, {
                "pathway_id": pathway_id,
                "pathway_label": entry.get("pathway_label", ""),
                "taxid": entry.get("taxid", ""),
                "inchikeys": set(),
            })
            bucket["inchikeys"].add(inchikey)

    rows: list[dict] = []
    for bucket in grouped.values():
        members = [by_inchikey[inchikey] for inchikey in bucket["inchikeys"]]
        log2fcs = [
            member["log2fc"] for member in members
            if member.get("log2fc") is not None
        ]
        significant = [member for member in members if member.get("significant")]
        n_up = sum(
            1 for member in significant if (member.get("log2fc") or 0.0) > 0)
        n_down = sum(
            1 for member in significant if (member.get("log2fc") or 0.0) < 0)
        n_members = len(members)
        tested = n_members >= min_members
        fisher_p = None
        if tested:
            fisher_p = fisher_exact_greater(
                k=len(significant), n=n_significant_total,
                K=n_members, N=n_background)
        wilcoxon_p, method = (
            wilcoxon_signed_rank(log2fcs) if tested else (None, "insufficient"))
        rows.append({
            "pathway_id": bucket["pathway_id"],
            "pathway_label": bucket["pathway_label"],
            "taxid": bucket["taxid"],
            "n_members": n_members,
            "n_significant": len(significant),
            "n_up": n_up,
            "n_down": n_down,
            "median_log2fc": _median(log2fcs),
            "direction": _direction(n_up, n_down),
            "fisher_p": fisher_p,
            "fisher_q": None,
            "tested": tested,
            "wilcoxon_p": wilcoxon_p,
            "wilcoxon_tested": method != "insufficient",
        })

    for row, q_value in zip(
            rows, benjamini_hochberg([row["fisher_p"] for row in rows])):
        row["fisher_q"] = q_value
    # Wilcoxon の BH は濃縮とは別に掛ける（合成しないため分けて補正する）。
    for row, q_value in zip(
            rows, benjamini_hochberg([row["wilcoxon_p"] for row in rows])):
        row["wilcoxon_q"] = q_value

    rows.sort(key=lambda row: (
        row["fisher_q"] if row["fisher_q"] is not None else 2.0,
        row["pathway_id"],
    ))
    summary = {
        "n_background": n_background,
        "n_significant": n_significant_total,
        "n_pathways_total": len(rows),
        "n_pathways_tested": sum(1 for row in rows if row["tested"]),
    }
    return rows, summary
