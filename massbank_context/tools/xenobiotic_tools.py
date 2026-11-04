"""外来の分子の走査の公開層。水準は判定ではなく材料で、解釈は会話側LLMが行う。"""
from massbank_context.tools.evidence_tools import _invoke
from massbank_context.xenobiotic.registry import register
from massbank_context.xenobiotic.screen import get_screen_rows, run_screen

__all__ = ['register_xenobiotics', 'screen_xenobiotic_features', 'get_xenobiotic_screen']


def register_xenobiotics(analysis_id: str, entries: list[dict]) -> str:
    """外来の分子（薬物・代謝物など）を解析に登録する。entries の欄に従い formula は必須で、名前は照合に使わない。登録し直しは全体の置き換え。"""
    return _invoke(lambda store: register(store, analysis_id, entries))


def screen_xenobiotic_features(analysis_id: str, test_group: str, ppm: float,
                               ppm_bands: list[dict] | None = None, polarity: str | None = None,
                               max_shifts: int = 2, output_path: str | None = None) -> str:
    """登録した外来の分子を全特徴量に対して走査する。test_group は comparison の group_a/group_b の名前で、v2 は polarity が必須。全行は TSV に書き戻り値は要約だけ。水準は判定ではなく材料。"""
    return _invoke(lambda store: run_screen(
        store, analysis_id, test_group=test_group, ppm=ppm, ppm_bands=ppm_bands,
        polarity=polarity, max_shifts=max_shifts, output_path=output_path))


def get_xenobiotic_screen(analysis_id: str, feature_ids: list[str] | None = None,
                          level: str | None = None, cursor: str | None = None,
                          limit: int = 50) -> str:
    """最新の走査の行を最大200行取得する。next_cursor を cursor へ渡す。走査し直すと古いカーソルは stale_cursor になる。"""
    return _invoke(lambda store: get_screen_rows(
        store, analysis_id, feature_ids=feature_ids, level=level, cursor=cursor, limit=limit))
