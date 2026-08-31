"""事実と仮説を分けた骨組み（spec §10）。

LLM を呼ばない。空のスロットを返すだけで、埋めるのは LLM である。
「活性化」「不活性化」の語は hypotheses の中にしか置けない。facts 側の語彙は
「上昇（group_b が高い）」「低下」で、群名を必ず伴わせる。
"""
from __future__ import annotations

#: facts 側に現れてはならない語。test で機械的に検査する。
INTERPRETATION_WORDS = ("活性化", "不活性化", "活性が上がった", "活性が下がった")


def _limits(coverage: float) -> list[str]:
    """毎回・固定で出す制約。省略可能にすると、いずれ省略される。"""
    return [
        f"パスウェイ被覆は {coverage * 100:.1f}% です。"
        "パスウェイが付かなかった化合物は「変化が無かった」のではなく「調べていない」ものです。",
        "PubChem のパスウェイは has participant であり、基質か産物か・反応の向き・"
        "量論を持ちません。代謝物の増減からパスウェイの活性は直接は導けません。",
        "化合物同定は MS-DIAL のライブラリ注釈由来で、MS/MS による独立検証は"
        "していません。",
        "メンバー数が閾値未満のパスウェイは検定していません（tested=false）。",
        "異性体は区別していません（MSI Level 2）。",
    ]


def build_findings(facts: dict, rows: list[dict], coverage: float) -> dict:
    """事実・仮説スロット・制約の 3 ブロックを組む。"""
    body = dict(facts)
    body["top_pathways"] = rows
    return {
        "facts": body,
        "hypotheses": [{
            "slot": "pathway_activity",
            "pathway_id": None,
            "claim": None,
            "evidence_rows": [],
            "required": ["pathway_id", "claim", "evidence_rows"],
            "note": ("活性化・不活性化の解釈はここにしか書けません。"
                     "claim には必ず evidence_rows（pathway_id の列）を伴わせてください。"),
        }],
        "limits": _limits(coverage),
    }
