"""パスウェイ系のツール。Phase 2〜3（spec §8）。"""
from __future__ import annotations

from pathlib import Path

from massbank_context.core.errors import contract_mismatch
from massbank_context.core.serialization import json_payload
from massbank_context.core.session import session
from massbank_context.experiment.contract import ContractError, parse_differential

__all__ = ["load_differential"]

_SIGN_NOTE = "log2fc は正なら group_b が高い（上昇）。massbank-context 側で反転しない。"


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
