"""パスウェイ系ツール。"""
import json

import pytest

from massbank_context.core import config
from massbank_context.core import session as session_module
from massbank_context.kg import cache
from massbank_context.tools.pathway_tools import load_differential, resolve_pathways

GOOD = """# contract_version = 1
# exported_at = 2026-08-31T12:00:00+00:00
# group_a = gf_AIN\tn_a = 5
# group_b = gf_HFD\tn_b = 5
# log2fc_sign = positive means group_b is higher
# q_threshold = 0.05\tlog2fc_threshold = 1.0\tlog_transform = true
# n_features_total = 10\tn_with_inchikey = 2\tn_unannotated = 8
spot_id\tname\tname_source\tontology\tinchikey\tinchikey_source\tmsi_level\tmz\trt\tlog2fc\tp_value\tq_value\tmean_a\tmean_b\tsignificant
1\tPC 34:1\tarf2\tPC\tAAAAAAAAAAAAAA-BBBBBBBBBB-C\tarf2\t\t760.5851\t12.34\t1.900000\t0.001\t0.01\t10\t40\ttrue
2\tCitric acid\tarf2\t\tKRKNYBCHXYNGOX-UHFFFAOYSA-N\tarf2\t\t191.0197\t2.10\t-0.500000\t0.4\t0.6\t20\t14\tfalse
"""


@pytest.fixture(autouse=True)
def clean_session(tmp_path, monkeypatch):
    """各テストを独立した可変セッションで実行する（キャッシュも毎回作り直す）。"""
    monkeypatch.setattr(config, "PATHWAY_CACHE_PATH",
                        tmp_path / "pathway_cache.sqlite")
    session_module.session.reset()
    yield
    session_module.session.reset()


def test_load_differential_reports_counts(tmp_path):
    """読み込み結果は行データではなく件数と比較の向きを返す。"""
    path = tmp_path / "differential.tsv"
    path.write_text(GOOD, encoding="utf-8")
    payload = json.loads(load_differential(str(path)))
    assert payload["n_rows"] == 2
    assert payload["n_inchikeys"] == 2
    assert payload["n_unannotated"] == 8
    assert payload["group_a"] == "gf_AIN"
    assert payload["group_b"] == "gf_HFD"
    assert "group_b" in payload["log2fc_sign"]
    assert session_module.session.differential is not None


def test_load_differential_rejects_bad_contract(tmp_path):
    """不正な契約はセッションを更新せず契約エラーとして返す。"""
    path = tmp_path / "differential.tsv"
    path.write_text(
        GOOD.replace("contract_version = 1", "contract_version = 99"),
        encoding="utf-8",
    )
    payload = json.loads(load_differential(str(path)))
    assert payload["error"]["code"] == "contract_mismatch"
    assert session_module.session.differential is None


def test_load_differential_missing_file(tmp_path):
    """存在しないファイルは契約エラーとして返す。"""
    payload = json.loads(load_differential(str(tmp_path / "nope.tsv")))
    assert payload["error"]["code"] == "contract_mismatch"
    assert "見つかりません" in payload["error"]["reason"]


def test_payload_does_not_carry_rows(tmp_path):
    """行そのものは payload に載せない（spec §12）。"""
    path = tmp_path / "differential.tsv"
    path.write_text(GOOD, encoding="utf-8")
    payload = json.loads(load_differential(str(path)))
    assert "rows" not in payload
    assert "AAAAAAAAAAAAAA" not in json.dumps(payload)


def test_load_differential_rejects_invalid_meta_without_replacing_session(tmp_path):
    """不正な整数メタは契約エラーにし、既存の結果を置き換えない。"""
    valid_path = tmp_path / "valid.tsv"
    valid_path.write_text(GOOD, encoding="utf-8")
    load_differential(str(valid_path))
    previous = session_module.session.differential
    previous_source = session_module.session.differential_source

    invalid_path = tmp_path / "invalid.tsv"
    invalid_path.write_text(
        GOOD.replace("n_features_total = 10", "n_features_total = not-a-number"),
        encoding="utf-8",
    )
    payload = json.loads(load_differential(str(invalid_path)))

    assert payload["error"]["code"] == "contract_mismatch"
    assert "n_features_total" in payload["error"]["reason"]
    assert session_module.session.differential is previous
    assert session_module.session.differential_source == previous_source


def _load(tmp_path, text=GOOD):
    path = tmp_path / "differential.tsv"
    path.write_text(text, encoding="utf-8")
    load_differential(str(path))


def _stub_adapter(hits, failed=None):
    class _Stub:
        def __init__(self):
            self.failed = list(failed or [])

        def fetch(self, inchikeys):
            return {ik: hits[ik] for ik in inchikeys if ik in hits}

        def fetch_organism_labels(self, taxids):
            return {"TAXID10090": "Mus musculus"}

    return _Stub()


def test_resolve_pathways_requires_differential():
    payload = json.loads(resolve_pathways())
    assert payload["error"]["code"] == "missing_state"
    assert "load_differential" in payload["error"]["required_tools"]


def test_resolve_pathways_reports_coverage(tmp_path, monkeypatch):
    _load(tmp_path)
    hits = {"KRKNYBCHXYNGOX-UHFFFAOYSA-N": [
        {"pathway_id": "PWID1", "pathway_label": "Citric Acid Cycle",
         "taxid": "TAXID10090"}]}
    monkeypatch.setattr("massbank_context.tools.pathway_tools._make_adapter",
                        lambda: _stub_adapter(hits))
    payload = json.loads(resolve_pathways())
    assert payload["n_queried"] == 2
    assert payload["n_with_pathway"] == 1
    assert payload["pathway_coverage"] == pytest.approx(0.5)
    assert payload["organisms"][0]["taxid"] == "TAXID10090"
    assert payload["organisms"][0]["label"] == "Mus musculus"
    assert payload["n_failed"] == 0


def test_resolve_pathways_reports_failures(tmp_path, monkeypatch):
    _load(tmp_path)
    monkeypatch.setattr(
        "massbank_context.tools.pathway_tools._make_adapter",
        lambda: _stub_adapter({}, failed=["KRKNYBCHXYNGOX-UHFFFAOYSA-N"]))
    payload = json.loads(resolve_pathways())
    assert payload["n_failed"] == 1
    assert "再実行" in payload["note"]


def test_resolve_pathways_reports_upstream_down(tmp_path, monkeypatch):
    """照会した全件が落ちたら、部分的な成功として返さず upstream_unavailable を返す。

    「0 件見つかった」と「エンドポイントが落ちていた」を同じ形で返すと、
    被覆 0% を結果として読まれる。
    """
    _load(tmp_path)
    monkeypatch.setattr(
        "massbank_context.tools.pathway_tools._make_adapter",
        lambda: _stub_adapter({}, failed=["AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                                          "KRKNYBCHXYNGOX-UHFFFAOYSA-N"]))
    payload = json.loads(resolve_pathways())
    assert payload["error"]["code"] == "upstream_unavailable"
    assert "rdfportal" in payload["error"]["endpoint"]


def test_resolve_pathways_is_cached(tmp_path, monkeypatch):
    _load(tmp_path)
    hits = {"KRKNYBCHXYNGOX-UHFFFAOYSA-N": [
        {"pathway_id": "PWID1", "pathway_label": "Citric Acid Cycle",
         "taxid": "TAXID10090"}]}
    monkeypatch.setattr("massbank_context.tools.pathway_tools._make_adapter",
                        lambda: _stub_adapter(hits))
    resolve_pathways()

    def _explode():
        raise AssertionError("2 回目はキャッシュで済むはず")

    monkeypatch.setattr("massbank_context.tools.pathway_tools._make_adapter",
                        _explode)
    payload = json.loads(resolve_pathways())
    assert payload["n_queried"] == 0
    assert payload["n_with_pathway"] == 1


def test_load_differential_returns_group_sizes_as_integers(tmp_path):
    """n_a / n_b は生文字列ではなく int で返す（設計書 §10 のサンプルは int）。"""
    path = tmp_path / "differential.tsv"
    path.write_text(GOOD, encoding="utf-8")
    payload = json.loads(load_differential(str(path)))
    assert payload["n_a"] == 5
    assert payload["n_b"] == 5
    assert isinstance(payload["n_a"], int)
    assert isinstance(payload["n_b"], int)


def test_load_differential_group_sizes_stay_none_when_absent(tmp_path):
    """書かれていない群サイズは 0 に丸めず None のままにする。"""
    text = GOOD.replace("# group_a = gf_AIN\tn_a = 5", "# group_a = gf_AIN")
    text = text.replace("# group_b = gf_HFD\tn_b = 5", "# group_b = gf_HFD")
    path = tmp_path / "differential.tsv"
    path.write_text(text, encoding="utf-8")
    payload = json.loads(load_differential(str(path)))
    assert payload["n_a"] is None
    assert payload["n_b"] is None


def test_load_differential_rejects_non_finite_number(tmp_path):
    """nan を含む TSV は読み込まず契約エラーにする（不正 JSON の入口を塞ぐ）。"""
    path = tmp_path / "differential.tsv"
    path.write_text(GOOD.replace("\t1.900000\t", "\tnan\t", 1), encoding="utf-8")
    payload = json.loads(load_differential(str(path)))
    assert payload["error"]["code"] == "contract_mismatch"
    assert "log2fc" in payload["error"]["reason"]
    assert session_module.session.differential is None


def _load_many(tmp_path, count):
    """InChIKey を count 件持つ差次的結果を読み込ませる。"""
    lines = GOOD.splitlines()[:8]
    for i in range(count):
        lines.append("\t".join([
            str(i), f"C{i}", "arf2", "", f"IK{i:04d}AAAAAA-BBBBBBBBBB-C",
            "arf2", "", "100.0", "1.0", "1.0", "0.001", "0.01", "1", "2",
            "true"]))
    path = tmp_path / "many.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    load_differential(str(path))


class _DeadAdapter:
    """常に全件失敗するアダプタ（上流が落ちている状態）。"""

    def __init__(self):
        self.failed = []
        self.batches = []

    def fetch(self, inchikeys):
        self.batches.append(list(inchikeys))
        self.failed.extend(inchikeys)
        return {}

    def fetch_organism_labels(self, taxids):
        raise AssertionError("上流が落ちているのにラベルを引くべきではない")


def test_resolve_pathways_breaks_out_on_dead_upstream(tmp_path, monkeypatch):
    """最初のバッチが全滅したら残りを照会せずに打ち切る。

    判定がループの後にあると、DEFAULT_TIMEOUT = 120 秒 × 残件を待ってから
    エラーを返すことになり（900 件なら 30 時間）、MCP ツール呼び出しとしては
    事実上ハングする。
    """
    _load_many(tmp_path, 120)
    dead = _DeadAdapter()
    monkeypatch.setattr("massbank_context.tools.pathway_tools._make_adapter",
                        lambda: dead)
    payload = json.loads(resolve_pathways())
    assert payload["error"]["code"] == "upstream_unavailable"
    # 50 件（1 バッチ）で打ち切る。残り 70 件は照会しない。
    assert len(dead.batches) == 1
    assert payload["error"]["n_processed"] == 50
    assert payload["error"]["n_remaining"] == 70


def test_resolve_pathways_continues_when_batch_partially_succeeds(
        tmp_path, monkeypatch):
    """一部だけ失敗するのは「上流が死んでいる」ではない。最後まで回す。"""
    _load_many(tmp_path, 120)

    class _Flaky:
        def __init__(self):
            self.failed = []
            self.batches = []

        def fetch(self, inchikeys):
            self.batches.append(list(inchikeys))
            self.failed.extend(inchikeys[1:])   # 先頭 1 件だけ成功する
            return {inchikeys[0]: [{"pathway_id": "PWID1",
                                    "pathway_label": "P",
                                    "taxid": "TAXID10090"}]}

        def fetch_organism_labels(self, taxids):
            return {"TAXID10090": "Mus musculus"}

    flaky = _Flaky()
    monkeypatch.setattr("massbank_context.tools.pathway_tools._make_adapter",
                        lambda: flaky)
    payload = json.loads(resolve_pathways())
    assert "error" not in payload
    assert len(flaky.batches) == 3          # 50 + 50 + 20
    assert payload["n_with_pathway"] == 3
    assert payload["n_failed"] == 117


class _LabelsOnlyAdapter:
    """パスウェイは引かず、ラベルだけ引けるアダプタ。"""

    def __init__(self):
        self.failed = []
        self.label_calls = []

    def fetch(self, inchikeys):
        raise AssertionError("全件キャッシュヒットなので照会は起きないはず")

    def fetch_organism_labels(self, taxids):
        self.label_calls.append(list(taxids))
        return {"TAXID10090": "Mus musculus"}


def test_resolve_pathways_fills_organism_label_on_full_cache_hit(
        tmp_path, monkeypatch):
    """全件キャッシュヒットの回でも未ラベルの taxid を埋める。

    ラベル解決が `if pending:` の中にしか無いと、未ラベルの taxid は
    refresh=True 以外では永久に埋まらない。埋まらない先は
    incompatible_organism.available の label で、設計書 §6.3 が
    「候補を示して呼び出し側に選ばせる」ための唯一の手掛かりである。
    """
    _load(tmp_path)
    conn = cache.open_cache(config.PATHWAY_CACHE_PATH)
    for inchikey in ("AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                     "KRKNYBCHXYNGOX-UHFFFAOYSA-N"):
        cache.record_result(conn, inchikey, config.PUBCHEM_SPARQL_ENDPOINT, [
            {"pathway_id": "PWID1", "pathway_label": "Citric Acid Cycle",
             "taxid": "TAXID10090"}])
    conn.close()

    adapter = _LabelsOnlyAdapter()
    monkeypatch.setattr("massbank_context.tools.pathway_tools._make_adapter",
                        lambda: adapter)
    payload = json.loads(resolve_pathways())
    assert payload["n_queried"] == 0
    assert adapter.label_calls == [["TAXID10090"]]
    assert payload["organisms"][0]["label"] == "Mus musculus"


def test_resolve_pathways_refuses_cache_with_other_version(tmp_path, monkeypatch):
    """版違いのキャッシュは開かず contract_mismatch を返す。"""
    _load(tmp_path)
    conn = cache.open_cache(config.PATHWAY_CACHE_PATH)
    conn.execute(
        "UPDATE cache_meta SET value = '99' WHERE key = 'contract_version'")
    conn.commit()
    conn.close()
    payload = json.loads(resolve_pathways())
    assert payload["error"]["code"] == "contract_mismatch"
    assert "99" in payload["error"]["reason"]
    assert str(config.PATHWAY_CACHE_PATH) == payload["error"]["path"]
