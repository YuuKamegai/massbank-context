"""入力順に依存しない根拠識別子。"""
import hashlib
import json


def canonical_json(value) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True,
                      separators=(',', ':'), allow_nan=False)


def stable_id(kind: str, parts: list[str]) -> str:
    return kind + ':' + hashlib.sha256(canonical_json(parts).encode()).hexdigest()
