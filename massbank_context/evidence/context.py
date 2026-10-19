"""研究者が確認した背景と提案段階の背景を区別する。"""
from massbank_context.core.evidence_ids import stable_id, canonical_json

FIELDS = {'species', 'tissue', 'treatment', 'groups', 'measurement', 'question'}


def set_context(store, analysis_id: str, context: dict, confirmed: bool) -> dict:
    if type(confirmed) is not bool or not isinstance(context, dict) or set(context) - FIELDS:
        raise ValueError('背景の項目または確認状態が不正です')
    status = 'confirmed' if confirmed else 'proposed'
    normalized = {}
    for key in sorted(FIELDS):
        item = context.get(key, {'value': None, 'source': ''})
        if not isinstance(item, dict) or 'value' not in item:
            raise ValueError('背景はvalue/sourceを持つオブジェクトで指定してください')
        if confirmed and item['value'] is not None and not item.get('source'):
            raise ValueError('確認済み背景には研究者由来の説明が必要です')
        normalized[key] = dict(value=item['value'], source=item.get('source', ''),
                               status=status if item['value'] is not None else 'unknown')
    result = dict(status=status, context=normalized)
    result['context_revision'] = stable_id('context', [analysis_id, canonical_json(result)])
    store.put_snapshot(analysis_id, 'context', result)
    return result
