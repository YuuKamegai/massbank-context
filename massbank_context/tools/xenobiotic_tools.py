"""外来の分子の走査の公開層。水準は判定ではなく材料で、解釈は会話側LLMが行う。"""
from massbank_context.tools.evidence_tools import _invoke
from massbank_context.xenobiotic.registry import register
from massbank_context.xenobiotic.screen import get_screen_rows, run_screen

__all__ = ['register_xenobiotics', 'screen_xenobiotic_features', 'get_xenobiotic_screen']


def register_xenobiotics(analysis_id: str, entries: list[dict]) -> str:
    """外来の分子（薬物・代謝物など）を解析に登録する。entries は辞書の配列で、欄は entry_id（必須、一意）、role（必須: drug / drug_metabolite / vehicle / anesthetic / co_medication / other）、name（必須、表示用で照合に使わない）、formula（必須、中性分子の組成式）、received_by（必須: test / reference / both。test は走査の test_group の側で、溶媒・麻酔薬は both）、inchikey（任意、無ければ MS/MS の筋は unavailable）、parent_entry_id（任意、代謝物の親の entry_id）。未知の欄は誤り。登録し直しは全体の置き換え。"""
    return _invoke(lambda store: register(store, analysis_id, entries))


def screen_xenobiotic_features(analysis_id: str, test_group: str, ppm: float,
                               ppm_bands: list[dict] | None = None, polarity: str | None = None,
                               max_shifts: int = 2, output_path: str | None = None) -> str:
    """登録した外来の分子を全特徴量に対して走査する。test_group は comparison の group_a/group_b の名前で、v2 は polarity が必須。ppm は既定の許容幅。ppm_bands は [{"min_mz": 数か null, "max_mz": 数か null, "ppm": 数}] の配列で、半開区間 [min_mz, max_mz) を理論 m/z で選び、どの帯にも入らなければ ppm を使う（未知のキーと帯の重なりは誤り）。全行は TSV に書き戻り値は要約だけ。水準は判定ではなく材料。spec §8 の較正は改訂した期待に合格したが手順書には入っていない（要約の calibration）。"""
    return _invoke(lambda store: run_screen(
        store, analysis_id, test_group=test_group, ppm=ppm, ppm_bands=ppm_bands,
        polarity=polarity, max_shifts=max_shifts, output_path=output_path))


def get_xenobiotic_screen(analysis_id: str, feature_ids: list[str] | None = None,
                          level: str | None = None, cursor: str | None = None,
                          limit: int = 50) -> str:
    """最新の走査の行を最大200行取得する。next_cursor を cursor へ渡す。走査し直すか level / feature_ids を変えると古いカーソルは stale_cursor になる。走査の後に登録し直していれば stale_registration が true。"""
    return _invoke(lambda store: get_screen_rows(
        store, analysis_id, feature_ids=feature_ids, level=level, cursor=cursor, limit=limit))
