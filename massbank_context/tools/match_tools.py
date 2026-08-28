"""照合ツール。Phase 1 ではバイアスを適用しない。"""
from __future__ import annotations

import sqlite3

from massbank_context.core import config
from massbank_context.core.errors import (
    incompatible_index,
    missing_index,
    missing_state,
)
from massbank_context.core.serialization import as_tsv, json_payload
from massbank_context.core.session import session
from massbank_context.index.writer import read_meta
from massbank_context.spectral.matcher import (
    DEFAULT_PRECURSOR_TOLERANCE,
    fragment_tolerance,
    match_spectrum,
)

__all__ = ["match_spectra"]

_COLUMNS = ["spectrum_id", "accession", "inchikey", "matched_fragments",
            "total_query_peaks", "weighted_score", "cosine",
            "unexplained_peak_ratio", "band", "precursor_delta"]

_DEFAULT_ACQUISITION_CLASS = "LC-ESI-MS2"

#: 1 回の戻り値に載せる候補行の上限。targets の既定は「読み込んだ全スペクトル」で、
#: 実物の .dcl は MSDecResult を数千本持つ。1 本あたり top_k 行が積まれるので、
#: 既定呼び出しがそのまま数万行になりうる — spec §12 が防ごうとしている
#: payload 肥大が、利便のための既定値から入ってくる形になる。
#: 打ち切ったときは必ず truncated で明示する（黙って切らない）。
#: 判定はスペクトル単位で行う。1 本分の候補を途中で切ると
#: 「上位 top_k」という行の意味が壊れるため、上限に達したかどうかは次の
#: スペクトルに進む前に見る。したがって実際の行数は上限を最大 top_k - 1 行
#: だけ超えうる（ソフトな上限）。
_MAX_CANDIDATE_ROWS = 200


def match_spectra(spectrum_ids: list[str] | None = None,
                  acquisition_class: str = _DEFAULT_ACQUISITION_CLASS,
                  precursor_tolerance: float = DEFAULT_PRECURSOR_TOLERANCE,
                  top_k: int = 10) -> str:
    """読み込んだスペクトルを索引と照合し、候補を TSV で返す。

    Phase 1 ではバイアス層が無いため bias.applied は常に false。
    objective 未設定はエラーにしない（spec §8）。
    """
    if not config.INDEX_PATH.exists():
        return missing_index(config.BUILD_COMMAND)
    if not session.spectra:
        return missing_state(
            "spectra", ["load_spectra"],
            "実験スペクトルが読み込まれていません。先に load_spectra を実行してください。")

    conn = sqlite3.connect(str(config.INDEX_PATH))
    try:
        meta = read_meta(conn)
        if meta is None:
            return missing_index(config.BUILD_COMMAND)
        available = meta.get("acquisition_classes", [])
        if acquisition_class not in available:
            return incompatible_index(
                acquisition_class, available,
                f"索引には {available} の層しかありません。測定原理の異なる層を"
                f"跨いだ照合は結果が無意味になるため行いません。")

        targets = (spectrum_ids if spectrum_ids is not None
                   else sorted(session.spectra))
        rows = []
        unknown_ids: list[str] = []
        reported = 0
        truncated = False
        for spectrum_id in targets:
            spectrum = session.spectra.get(spectrum_id)
            if spectrum is None:
                # 読み込まれていない id を黙って飛ばしたうえで n_spectra に
                # 数えると、打ち間違えた id が「照合したが何も当たらなかった」
                # と読めてしまう。数えず、id をそのまま返して区別させる。
                unknown_ids.append(spectrum_id)
                continue
            if len(rows) >= _MAX_CANDIDATE_ROWS:
                truncated = True
                break
            for candidate in match_spectrum(
                    conn, spectrum, acquisition_class,
                    precursor_tolerance=precursor_tolerance, top_k=top_k):
                rows.append([
                    spectrum_id, candidate.accession, candidate.inchikey,
                    candidate.matched_fragments, candidate.total_query_peaks,
                    round(candidate.weighted_score, 4),
                    round(candidate.cosine, 4),
                    round(candidate.unexplained_peak_ratio, 4),
                    candidate.band,
                    None if candidate.precursor_delta is None
                    else round(candidate.precursor_delta, 4),
                ])
            reported += 1
        bin_width = meta["bin_width"]
    finally:
        conn.close()

    payload = {
        "acquisition_class": acquisition_class,
        # 実際に照合したスペクトル数。要求した数ではない（finding 9）。
        "n_spectra": reported,
        "n_spectra_requested": len(targets),
        "n_candidates": len(rows),
        # この 2 つが無いと数字を再現・検算できない。fragment_tolerance は
        # ±1 ビン探索を Da に直したもので、他のどこにも出てこない（finding 13）。
        "precursor_tolerance": precursor_tolerance,
        "fragment_tolerance": round(fragment_tolerance(bin_width), 6),
        "candidates": as_tsv(_COLUMNS, rows),
        "bias": {
            "applied": False,
            "reason": "Phase 1 ではバイアス層が未実装。並び順は weighted_score のみで決まる。",
        },
        "note": "band が PRECURSOR_ONLY の候補は精密質量しか合っていない。"
                "PASS と同一視しないこと。",
    }
    if unknown_ids:
        payload["unknown_spectrum_ids"] = unknown_ids
    if truncated:
        payload["truncated"] = True
        payload["truncated_reason"] = (
            f"候補行が上限 {_MAX_CANDIDATE_ROWS} 行に達したため、"
            f"{len(targets)} 本中 {reported} 本までで打ち切りました"
            "（打ち切りはスペクトル単位なので n_candidates は上限を"
            "わずかに超えることがあります）。"
            "spectrum_ids で対象を絞るか top_k を下げてください。")
    return json_payload(payload)
