"""外来の分子の登録: 検証・索引の情報・保存。"""
from __future__ import annotations

import sqlite3
from dataclasses import asdict, dataclass
from typing import Sequence

from massbank_context.core.evidence_ids import canonical_json, stable_id
from massbank_context.core.formula import FormulaError, hill_formula, monoisotopic_mass, parse_formula
from massbank_context.core.serialization import as_tsv
from massbank_context.evidence.store import EvidenceStore
from massbank_context.experiment.bundle import INCHIKEY_PATTERN
from massbank_context.xenobiotic import msms

__all__ = ["XenobioticError", "ROLES", "RECEIVED_BY", "ENTRY_KEYS", "XenoEntry",
           "validate_entries", "describe_entries", "register", "load_registration"]

ROLES = ("drug", "drug_metabolite", "vehicle", "anesthetic", "co_medication", "other")
RECEIVED_BY = ("test", "reference", "both")
ENTRY_KEYS = {"entry_id", "role", "name", "inchikey", "formula", "parent_entry_id", "received_by"}

SNAPSHOT_KIND = "xenobiotics"
#: compounds.exact_mass と計算値の許容差（Da）
MASS_TOLERANCE_DA = 0.002


class XenobioticError(ValueError):
    """外来の分子の要求が不正。code は機械可読（_invoke が拾う）。"""

    def __init__(self, message: str, code: str = "invalid_xenobiotic_request"):
        super().__init__(message)
        self.code = code


@dataclass(frozen=True)
class XenoEntry:
    entry_id: str
    role: str
    name: str
    inchikey: str | None
    formula: str  # Hill 表記に正規化済み
    composition: dict
    monoisotopic_mass: float
    parent_entry_id: str | None
    received_by: str


def _text(raw: dict, key: str, *, required: bool = False) -> str | None:
    """空文字・None は未指定として扱う。文字列以外は拒否する。"""
    value = raw.get(key)
    if value is not None and not isinstance(value, str):
        raise XenobioticError(f"{key} は文字列でなければなりません: {value!r}")
    value = (value or "").strip() or None
    if value is None and required:
        raise XenobioticError(f"{key} は必須です: {raw!r}")
    return value


def validate_entries(entries: list[dict]) -> list[XenoEntry]:
    """登録要求を検証し XenoEntry にする。不正なら XenobioticError。"""
    if not isinstance(entries, list) or not entries:
        raise XenobioticError("entries は 1 件以上の一覧でなければなりません")
    result: list[XenoEntry] = []
    for raw in entries:
        if not isinstance(raw, dict):
            raise XenobioticError(f"entry は辞書でなければなりません: {raw!r}")
        extra = set(raw) - ENTRY_KEYS
        if extra:
            raise XenobioticError(f"未知のキー: {sorted(extra)}（使えるのは {sorted(ENTRY_KEYS)}）")
        entry_id = _text(raw, "entry_id", required=True)
        role = _text(raw, "role", required=True)
        if role not in ROLES:
            raise XenobioticError(f"role は {ROLES} のいずれか: {role!r}")
        received_by = _text(raw, "received_by") or "test"
        if received_by not in RECEIVED_BY:
            raise XenobioticError(f"received_by は {RECEIVED_BY} のいずれか: {received_by!r}")
        inchikey = _text(raw, "inchikey")
        if inchikey is not None and not INCHIKEY_PATTERN.match(inchikey):
            raise XenobioticError(f"InChIKey の形式が不正: {inchikey!r}")
        formula = _text(raw, "formula", required=True)
        try:
            composition = parse_formula(formula)
        except FormulaError as exc:
            raise XenobioticError(f"組成式が不正（{entry_id}）: {exc}") from exc
        result.append(XenoEntry(
            entry_id=entry_id, role=role, name=_text(raw, "name") or "", inchikey=inchikey,
            formula=hill_formula(composition), composition=composition,
            monoisotopic_mass=monoisotopic_mass(composition),
            parent_entry_id=_text(raw, "parent_entry_id"), received_by=received_by))
    ids = [e.entry_id for e in result]
    if len(set(ids)) != len(ids):
        raise XenobioticError(f"entry_id が重複しています: {ids}")
    for entry in result:
        if entry.parent_entry_id is None:
            continue
        if entry.parent_entry_id == entry.entry_id:
            raise XenobioticError(f"parent_entry_id が自分自身です: {entry.entry_id}")
        if entry.parent_entry_id not in ids:
            raise XenobioticError(
                f"parent_entry_id が登録に無い: {entry.entry_id} -> {entry.parent_entry_id}")
    return result


def _counts(rows) -> str:
    return ",".join(f"{value}:{n}" for value, n in rows)


def describe_entries(conn: sqlite3.Connection, entries: Sequence[XenoEntry]) -> list[dict]:
    """各 entry に索引上の情報（LC-ESI-MS2 の件数・極性・付加体）と警告を付ける。"""
    out = []
    for entry in entries:
        info = {"entry_id": entry.entry_id, "n_records": 0, "polarities": "",
                "precursor_types": "", "warnings": ""}
        warnings: list[str] = []
        if entry.inchikey is None:
            warnings.append("no_inchikey")
        else:
            compound = conn.execute(
                "SELECT molecular_formula, exact_mass FROM compounds WHERE inchikey=?",
                (entry.inchikey,)).fetchone()
            base = "FROM records WHERE inchikey=? AND acquisition_class=?"
            args = (entry.inchikey, msms.ACQUISITION_CLASS)
            info["n_records"] = conn.execute(f"SELECT COUNT(*) {base}", args).fetchone()[0]
            any_record = conn.execute(
                "SELECT 1 FROM records WHERE inchikey=? LIMIT 1", (entry.inchikey,)).fetchone()
            if compound is None and any_record is None:
                warnings.append("not_in_index")
            elif info["n_records"] == 0:
                warnings.append("no_lc_esi_ms2_records")
            info["polarities"] = _counts(conn.execute(
                f"SELECT polarity, COUNT(*) {base} GROUP BY polarity ORDER BY polarity", args))
            info["precursor_types"] = _counts(conn.execute(
                f"SELECT precursor_type, COUNT(*) {base} GROUP BY precursor_type"
                " ORDER BY COUNT(*) DESC, precursor_type", args))
            if compound is not None:
                index_formula, exact_mass = compound
                if index_formula:
                    try:
                        index_hill = hill_formula(parse_formula(index_formula))
                    except FormulaError:
                        index_hill = index_formula
                    if index_hill != entry.formula:
                        warnings.append(f"formula_mismatch_index:{index_hill}")
                if exact_mass is not None:
                    diff = abs(exact_mass - entry.monoisotopic_mass)
                    if diff > MASS_TOLERANCE_DA:
                        warnings.append(f"mass_mismatch_index:{diff:.4f}")
        info["warnings"] = ";".join(warnings)
        out.append(info)
    return out


_COLUMNS = ["entry_id", "role", "name", "inchikey", "formula", "monoisotopic_mass",
            "parent_entry_id", "received_by", "n_records", "polarities", "precursor_types",
            "warnings"]


def register(store: EvidenceStore, analysis_id: str, entries: list[dict]) -> dict:
    """外来の分子を解析に登録する。登録し直しは全体の置き換え。"""
    _, base_revision, _ = store.read_view(analysis_id)  # 解析が無ければ ValueError
    validated = validate_entries(entries)
    conn, _meta = msms.open_index()
    try:
        described = describe_entries(conn, validated)
    finally:
        conn.close()
    as_dicts = [asdict(e) for e in validated]
    registration_id = stable_id("xenobiotics", [analysis_id, canonical_json(as_dicts)])
    store.put_snapshot(analysis_id, SNAPSHOT_KIND,
                       {"registration_id": registration_id, "entries": as_dicts},
                       expected_revision=base_revision)
    rows = []
    for entry, info in zip(validated, described):
        values = {**asdict(entry), **info, "monoisotopic_mass": round(entry.monoisotopic_mass, 4)}
        rows.append([values[c] for c in _COLUMNS])
    return {"registration_id": registration_id, "entries": as_tsv(_COLUMNS, rows),
            "n_entries": len(validated)}


def load_registration(store: EvidenceStore, analysis_id: str) -> tuple[list[XenoEntry], str]:
    """登録済みの (entries, registration_id)。無ければ missing_state。"""
    snapshot = store.current(analysis_id, SNAPSHOT_KIND)
    if not snapshot:
        raise XenobioticError("外来の分子が登録されていません。先に登録してください",
                              code="missing_state")
    return [XenoEntry(**e) for e in snapshot["entries"]], snapshot["registration_id"]
