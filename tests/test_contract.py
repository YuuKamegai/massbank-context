"""契約 TSV のリーダ（spec §6.1）。"""
import pytest

from massbank_context.experiment.contract import ContractError, parse_differential

GOOD = """# contract_version = 1
# exported_at = 2026-08-31T12:00:00+00:00
# source_arf = C:/x/AlignmentResult_2026_01_01_00_00_00.arf
# source_arf2 = C:/x/AlignmentResult_2026_01_01_00_00_00.arf2
# group_a = gf_AIN\tn_a = 5
# group_b = gf_HFD\tn_b = 5
# log2fc_sign = positive means group_b is higher
# q_threshold = 0.05\tlog2fc_threshold = 1.0\tlog_transform = true
# n_features_total = 10\tn_with_inchikey = 2\tn_unannotated = 8
spot_id\tname\tname_source\tontology\tinchikey\tinchikey_source\tmsi_level\tmz\trt\tlog2fc\tp_value\tq_value\tmean_a\tmean_b\tsignificant
1\tPC 34:1\tarf2\tPC\tAAAAAAAAAAAAAA-BBBBBBBBBB-C\tarf2\t\t760.5851\t12.34\t1.900000\t0.001\t0.01\t10\t40\ttrue
2\tCitric acid\tarf2\t\tKRKNYBCHXYNGOX-UHFFFAOYSA-N\tarf2\t\t191.0197\t2.10\t-0.500000\t0.4\t0.6\t20\t14\tfalse
"""


def test_parses_meta_and_rows():
    """メタ情報と型変換済みのデータ行を返す。"""
    result = parse_differential(GOOD)
    assert result.group_a == "gf_AIN"
    assert result.group_b == "gf_HFD"
    assert result.n_features_total == 10
    assert result.n_unannotated == 8
    assert len(result.rows) == 2
    assert result.rows[0]["inchikey"] == "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
    assert result.rows[0]["log2fc"] == pytest.approx(1.9)
    assert result.rows[0]["significant"] is True
    assert result.rows[1]["significant"] is False


# ms-data-parser の統計エンジン（dataset_compare）の結果を dataset_export_differential が書くと、
# msi_note の後ろに足される 6 行（export_contract.stats_meta_lines。契約の版は 1 のまま）。
# 値は ms-data-parser main 4a4eee3 の build_meta + stats_meta_lines の出力の形をそのまま写した。
STATS_META_LINES = (
    "# msi_level_note = from mzTab-M\n"
    "# test = moderated_t\teffect_size_definition = log2_mean_difference\tci_level = 0.95\n"
    "# verdict_q_threshold = 0.05\tequivalence_fc = 1.2\tbh_population = matrix\tbh_n = 10"
    "\tq_interpretation = discovery_fdr\n"
    "# unit = subject\tn_units_a = 4\tn_units_b = 5\n"
    "# mean_definition = linear arithmetic mean over units after preprocessing\n"
    "# matrix_id = mx_51f8258026b883d1\trule_hash = 7309dc29\tcontrast_id = gf_HFD_vs_gf_AIN\n"
    "# significant_rule = q<=q_threshold & |log2fc|>=log2fc_threshold & verdict=changed\n"
)


def test_stats_engine_meta_lines_do_not_disturb_contract():
    """統計エンジンのメタ行が足されても、既存のメタ・行の読みは変わらず、新しいキーは素のまま残る。"""
    marker = "# n_features_total = 10\tn_with_inchikey = 2\tn_unannotated = 8\n"
    text = GOOD.replace(marker, marker + STATS_META_LINES)
    result = parse_differential(text)
    base = parse_differential(GOOD)
    assert result.rows == base.rows
    assert {k: result.meta[k] for k in base.meta} == base.meta
    # n_a / n_b は注入の数のまま。比較の単位の数は別のキー
    assert result.integer_meta("n_a") == 5
    assert result.meta["n_units_a"] == "4"
    assert result.meta["test"] == "moderated_t"
    # 値の中の `<=` `>=` `=` で key が切れない（先頭の `=` だけで分ける）
    assert result.meta["significant_rule"] == (
        "q<=q_threshold & |log2fc|>=log2fc_threshold & verdict=changed")


def test_rejects_unknown_contract_version():
    """未対応のバージョンを黙って解釈しない。"""
    text = GOOD.replace("contract_version = 1", "contract_version = 99")
    with pytest.raises(ContractError) as excinfo:
        parse_differential(text)
    assert "contract_version" in str(excinfo.value)


def test_rejects_missing_sign_meta():
    """符号の向きが書かれていないファイルは読まない。"""
    text = "\n".join(
        line for line in GOOD.splitlines() if not line.startswith("# log2fc_sign")
    )
    with pytest.raises(ContractError) as excinfo:
        parse_differential(text)
    assert "log2fc_sign" in str(excinfo.value)


def test_rejects_unexpected_sign_direction():
    """group_a を正とする異なる符号規約を拒否する。"""
    text = GOOD.replace(
        "positive means group_b is higher", "positive means group_a is higher"
    )
    with pytest.raises(ContractError):
        parse_differential(text)


def test_rejects_missing_column():
    """下流処理に必要な InChIKey 列の欠落を拒否する。"""
    text = GOOD.replace("\tinchikey\t", "\t")
    with pytest.raises(ContractError) as excinfo:
        parse_differential(text)
    assert "inchikey" in str(excinfo.value)


def test_rejects_empty_body():
    """ヘッダだけ、またはデータ行がないファイルを拒否する。"""
    text = "\n".join(GOOD.splitlines()[:10])
    with pytest.raises(ContractError):
        parse_differential(text)


@pytest.mark.parametrize(
    ("original", "replacement", "column"),
    [
        ("\t760.5851\t", "\tinvalid\t", "mz"),
        ("\t12.34\t", "\tinvalid\t", "rt"),
        ("\t1.900000\t", "\tinvalid\t", "log2fc"),
        ("\t0.001\t", "\tinvalid\t", "p_value"),
        ("\t0.01\t", "\tinvalid\t", "q_value"),
        ("\t10\t40\ttrue", "\tinvalid\t40\ttrue", "mean_a"),
        ("\t20\t14\tfalse", "\t20\tinvalid\tfalse", "mean_b"),
    ],
)
def test_rejects_non_numeric_required_numeric_cell(original, replacement, column):
    """非空の必須数値セルを欠損値として扱わない。"""
    text = GOOD.replace(original, replacement, 1)
    with pytest.raises(ContractError) as excinfo:
        parse_differential(text)
    assert column in str(excinfo.value)


def test_keeps_empty_numeric_cell_as_none():
    """空の必須数値セルは既存契約どおり None として読み込む。"""
    result = parse_differential(GOOD.replace("\t760.5851\t", "\t\t", 1))
    assert result.rows[0]["mz"] is None


def test_rejects_unknown_significant_value():
    """真偽値でない significant を false に丸め込まない。"""
    text = GOOD.replace("\ttrue\n", "\tunknown\n", 1)
    with pytest.raises(ContractError) as excinfo:
        parse_differential(text)
    assert "significant" in str(excinfo.value)


@pytest.mark.parametrize(
    ("original", "replacement", "key"),
    [
        ("n_features_total = 10", "n_features_total = not-a-number", "n_features_total"),
        ("n_unannotated = 8", "n_unannotated = -1", "n_unannotated"),
        ("n_a = 5", "n_a = 1.5", "n_a"),
        ("n_b = 5", "n_b = -1", "n_b"),
        ("n_with_inchikey = 2", "n_with_inchikey = invalid", "n_with_inchikey"),
    ],
)
def test_rejects_invalid_integer_meta(original, replacement, key):
    """集計用メタデータを parse 中に非負整数として検証する。"""
    with pytest.raises(ContractError) as excinfo:
        parse_differential(GOOD.replace(original, replacement, 1))
    assert key in str(excinfo.value)


@pytest.mark.parametrize("literal", ["nan", "NaN", "inf", "-inf", "Infinity"])
def test_rejects_non_finite_numeric_cell(literal):
    """nan / inf を契約エラーにする（float() は黙って受理してしまう）。

    上流の BH 補正や t 検定は定数入力で普通に nan を吐き、TSV には `nan` と
    書かれる。通すと median_log2fc が nan になり、`NaN > 0` も `NaN < 0` も
    False なので n_up / n_down が n_significant と食い違い、最終的に
    RFC 8259 でない JSON（`{"x":NaN}`）が出ていく。
    """
    text = GOOD.replace("\t1.900000\t", f"\t{literal}\t", 1)
    with pytest.raises(ContractError) as excinfo:
        parse_differential(text)
    assert "log2fc" in str(excinfo.value)


@pytest.mark.parametrize("column", ["mz", "rt", "p_value", "q_value"])
def test_non_finite_is_rejected_in_every_numeric_column(column):
    """log2fc だけでなく全ての数値列で非有限値を弾く。"""
    replacements = {
        "mz": ("\t760.5851\t", "\tnan\t"),
        "rt": ("\t12.34\t", "\tinf\t"),
        "p_value": ("\t0.001\t", "\tnan\t"),
        "q_value": ("\t0.01\t", "\t-inf\t"),
    }
    original, replacement = replacements[column]
    with pytest.raises(ContractError) as excinfo:
        parse_differential(GOOD.replace(original, replacement, 1))
    assert column in str(excinfo.value)


def test_integer_meta_returns_int_not_string():
    """任意の整数メタは int で返す（payload に生文字列を出さない）。"""
    result = parse_differential(GOOD)
    assert result.integer_meta("n_a") == 5
    assert isinstance(result.integer_meta("n_a"), int)
    assert result.integer_meta("n_b") == 5


def test_integer_meta_is_none_when_absent():
    """書かれていない値は 0 ではなく None（「0 人」と「未記載」は別物）。"""
    text = "\n".join(
        line for line in GOOD.splitlines()
        if not line.startswith("# group_a")
    ).replace("# group_b = gf_HFD\tn_b = 5",
              "# group_a = gf_AIN\n# group_b = gf_HFD")
    result = parse_differential(text)
    assert result.integer_meta("n_a") is None
    assert result.integer_meta("n_b") is None
