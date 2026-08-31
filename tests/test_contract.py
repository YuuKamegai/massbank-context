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
