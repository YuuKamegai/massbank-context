"""ms-data-parser が書き出す差次的結果を読む（spec §6.1 の契約）。

契約に合わないファイルは**読まずに落とす**。特に log2fc の向きは、書かれて
いなければ読まない。向きが分からないまま読むと、上昇と低下が静かに入れ替わる。
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

SUPPORTED_CONTRACT_VERSION = "1"
EXPECTED_SIGN = "positive means group_b is higher"

REQUIRED_COLUMNS = [
    "spot_id", "name", "ontology", "inchikey", "mz", "rt",
    "log2fc", "p_value", "q_value", "mean_a", "mean_b", "significant",
]
REQUIRED_META = [
    "contract_version", "group_a", "group_b", "log2fc_sign",
    "n_features_total", "n_unannotated",
]
REQUIRED_INTEGER_META = ["n_features_total", "n_unannotated"]
OPTIONAL_INTEGER_META = ["n_a", "n_b", "n_with_inchikey"]


class ContractError(ValueError):
    """契約に合わない入力。呼び出し側が contract_mismatch に変換する。"""


@dataclass
class DifferentialSet:
    """検証済みの差次的結果と、その来歴を保持する。"""

    meta: dict[str, str] = field(default_factory=dict)
    rows: list[dict] = field(default_factory=list)

    @property
    def group_a(self) -> str:
        """比較元グループ名を返す。"""
        return self.meta["group_a"]

    @property
    def group_b(self) -> str:
        """比較先グループ名を返す。"""
        return self.meta["group_b"]

    @property
    def n_features_total(self) -> int:
        """入力特徴量の総数を返す。"""
        return int(self.meta["n_features_total"])

    @property
    def n_unannotated(self) -> int:
        """注釈されず書き出されなかった特徴量数を返す。"""
        return int(self.meta["n_unannotated"])

    def integer_meta(self, key: str) -> int | None:
        """任意の整数メタを int で返す。キーが無ければ None を返す。

        値は parse_differential が OPTIONAL_INTEGER_META で検証済みなので、
        ここでの int() は失敗しない。payload に生文字列を出すと、設計書 §10 の
        サンプル（int）と型が食い違う。**無い値は 0 ではなく None** である
        （「群サイズが 0」と「群サイズが書かれていない」は別物）。
        """
        value = self.meta.get(key)
        return None if value is None else int(value)

    @property
    def inchikeys(self) -> list[str]:
        """重複を除いた InChIKey を決定的な順序で返す。"""
        return sorted({row["inchikey"] for row in self.rows})


def _parse_meta_line(line: str, meta: dict[str, str]) -> None:
    """`# key = value` を、タブ区切りで複数個まで拾う。"""
    for chunk in line.lstrip("#").split("\t"):
        if "=" not in chunk:
            continue
        key, _, value = chunk.partition("=")
        meta[key.strip()] = value.strip()


def _to_float(value: str, column: str) -> float | None:
    """空欄を None として数値列を変換し、不正な値は拒否する。

    `nan` / `inf` は float() が黙って受理するが、ここでは**契約違反として落とす**。
    上流の BH 補正や t 検定は定数入力で普通に nan を吐き、TSV には `nan` と書かれる。
    これを通すと (a) median_log2fc が nan になり（NaN の比較は全て False なので
    sorted の中央に居座る）、(b) `NaN > 0` も `NaN < 0` も False なので n_up / n_down が
    n_significant と食い違い、(c) 最終的に payload が `{"x":NaN}` になる。
    これは RFC 8259 の JSON ではなく、厳密なパーサを持つクライアントが落ちる。
    「欠測」を表したいなら空欄で書くのが契約である。
    """
    value = (value or "").strip()
    if not value:
        return None
    try:
        number = float(value)
    except ValueError as error:
        raise ContractError(f"{column} が数値ではありません: {value!r}") from error
    if not math.isfinite(number):
        raise ContractError(
            f"{column} が有限の数値ではありません: {value!r}"
            "（欠測は空欄で書いてください）"
        )
    return number


def _validate_nonnegative_integer_meta(meta: dict[str, str], key: str) -> None:
    """集計用メタデータを非負整数として検証する。"""
    value = meta[key]
    try:
        integer = int(value)
    except ValueError as error:
        raise ContractError(f"{key} が整数ではありません: {value!r}") from error
    if integer < 0:
        raise ContractError(f"{key} は 0 以上の整数でなければなりません: {value!r}")


def parse_differential(text: str) -> DifferentialSet:
    """契約 TSV を読む。契約に合わなければ ContractError を送出する。"""
    meta: dict[str, str] = {}
    body: list[str] = []
    for line in text.splitlines():
        if line.startswith("#"):
            _parse_meta_line(line, meta)
        elif line.strip():
            body.append(line)

    missing_meta = [key for key in REQUIRED_META if key not in meta]
    if missing_meta:
        raise ContractError(f"メタ行が足りません: {', '.join(missing_meta)}")
    if meta["contract_version"] != SUPPORTED_CONTRACT_VERSION:
        raise ContractError(
            f"未対応の contract_version です: {meta['contract_version']}"
            f"（対応: {SUPPORTED_CONTRACT_VERSION}）"
        )
    if meta["log2fc_sign"] != EXPECTED_SIGN:
        raise ContractError(
            f"log2fc_sign が想定と違います: {meta['log2fc_sign']!r}"
            f"（想定: {EXPECTED_SIGN!r}）。向きが違うファイルは読みません。"
        )
    for key in REQUIRED_INTEGER_META:
        _validate_nonnegative_integer_meta(meta, key)
    for key in OPTIONAL_INTEGER_META:
        if key in meta:
            _validate_nonnegative_integer_meta(meta, key)
    if len(body) < 2:
        raise ContractError("ヘッダ行かデータ行がありません。")

    columns = body[0].split("\t")
    missing_columns = [column for column in REQUIRED_COLUMNS if column not in columns]
    if missing_columns:
        raise ContractError(f"必須列がありません: {', '.join(missing_columns)}")

    index = {name: position for position, name in enumerate(columns)}
    rows: list[dict] = []
    for line in body[1:]:
        cells = line.split("\t")
        if len(cells) != len(columns):
            raise ContractError(
                f"列数が合いません（期待 {len(columns)}、実際 {len(cells)}）: {line[:60]}"
            )
        inchikey = cells[index["inchikey"]].strip()
        if not inchikey:
            raise ContractError(
                "inchikey が空の行があります。書き出し側で落とすべき行です。"
            )
        significant = cells[index["significant"]].strip().lower()
        if significant not in {"true", "false"}:
            raise ContractError(
                "significant は true または false でなければなりません: "
                f"{significant!r}"
            )
        rows.append({
            "spot_id": cells[index["spot_id"]].strip(),
            "name": cells[index["name"]].strip(),
            "ontology": cells[index["ontology"]].strip(),
            "inchikey": inchikey,
            "mz": _to_float(cells[index["mz"]], "mz"),
            "rt": _to_float(cells[index["rt"]], "rt"),
            "log2fc": _to_float(cells[index["log2fc"]], "log2fc"),
            "p_value": _to_float(cells[index["p_value"]], "p_value"),
            "q_value": _to_float(cells[index["q_value"]], "q_value"),
            "mean_a": _to_float(cells[index["mean_a"]], "mean_a"),
            "mean_b": _to_float(cells[index["mean_b"]], "mean_b"),
            "significant": significant == "true",
        })
    return DifferentialSet(meta=meta, rows=rows)
