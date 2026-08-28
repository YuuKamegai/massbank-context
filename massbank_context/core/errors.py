"""クライアント非依存のエラーエンベロープ（spec §8）。

MCP クライアントが「次に何をすればよいか」をサーバ固有の日本語文面を解釈せずに
機械的に判断できるようにする。本文の JSON に載せるのは、使用中の MCP SDK では
isError=true の経路が structuredContent を捨てるため。

依存グラフの leaf（serialization のみ）。
"""
from __future__ import annotations

from massbank_context.core.serialization import json_payload


def missing_index(required_action: str) -> str:
    """索引が無い。ツールでは直せないので CLI を案内する。"""
    if not required_action:
        raise ValueError("required_action は必須です。")
    return json_payload({
        "error": {
            "code": "missing_index",
            "required_action": required_action,
        }
    })


def missing_state(state: str, required_tools: list[str], message: str) -> str:
    """前提となるセッション状態が無い。ツールで復旧できる。"""
    if not state or not required_tools or not message:
        raise ValueError("state / required_tools / message はいずれも必須です。")
    return json_payload({
        "error": {
            "code": "missing_state",
            "state": state,
            "required_tools": list(required_tools),
            "message": message,
        }
    })


def incompatible_index(query_acquisition_class: str,
                       available: list[str],
                       message: str) -> str:
    """索引の層が要求と合わない（spec §6.1）。

    測定原理の異なるスペクトルを混ぜて照合することは無意味なので、黙って
    別の層で代用せず明示的に拒否する。
    """
    if not query_acquisition_class or not message:
        raise ValueError("query_acquisition_class / message は必須です。")
    return json_payload({
        "error": {
            "code": "incompatible_index",
            "query_acquisition_class": query_acquisition_class,
            "available": list(available),
            "message": message,
        }
    })
