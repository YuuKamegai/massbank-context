# 根拠バンドル v3

v2（[evidence-bundle-v2.md](evidence-bundle-v2.md)）に、試料ごとの値と検出を足したもの。ここには v2 との差分だけを書く。
massbank-context は 2 と 3 の両方を読む（v2 は samples / intensities を持たない）。

## manifest

- contract_version は整数 3。
- `comparison.polarity` は必須で `positive` / `negative`。v2 でも書いてあれば同じ検証をする。
- 最上位の `intensities_stage` は、役割 `intensities` があるとき必須で
  `pre_imputation`（補完の前の値）/ `post_imputation`（補完の後の値）。
- `population.detection_unavailable_reason` は任意（文字列または null）。`detected` を取り込めなかった理由
  （`sample_axis_mismatch`・`mz_mismatch` など）を書く。検出が取れていれば null。

## 役割

`samples` と `intensities` は v3 のときだけ各 1 個まで許し、**両方そろえるか両方省く**。v2 に書くと拒否する。

- samples: `sample_id, group, role, injection_order`
  - `group` は比較の group_a の試料が `a`、group_b の試料が `b`、比較に入らない試料（QC・ブランクなど）は空。
  - `sample_id` は空でなく重複しない。`role` は自由（`sample` / `qc` / `blank` など。不明は空）。
  - `injection_order` は整数または空。
  - `group` が `a` / `b` の試料数は `comparison.n_a` / `n_b` と一致させる。
- intensities: `feature_id, sample_id, value, detected`（縦持ち）
  - **全特徴 × group が `a` / `b` の試料**を過不足なく 1 行ずつ書く。欠け・重複は拒否する。
  - `group` が空の試料（QC など）の行は書かない。書くと拒否する（黙って捨てない）。
  - `feature_id` は features に、`sample_id` は samples にあること。
  - `value` は空か有限の数。負の値も可（対数変換後の値があるため）。段は `intensities_stage` で示す。
  - `detected` は `true`（実測の検出）/ `false`（未検出・gap-fill・補完）/ 空（不明）。

## features の任意列

`adduct`（MS-DIAL の付加体表記。例 `[M+H]1+`）。無ければ列ごと省くか空にする。空は読み込み時に `None`。

## 見本

見本のバンドルの原本は ms-data-parser の `tests/fixtures/evidence_bundle_v3/`（`tests/evidence_bundle_fixture.py` で
作り直す）。本リポジトリの写しは `tests/fixtures/evidence_v3/` で、`tests/test_evidence_bundle_v3.py` の
`test_shared_fixture_from_ms_data_parser_loads` が読む。両リポのテストが同じ見本を読むことで、契約の食い違いを検出する。
manifest がファイルの sha256 を持つので、写しは改行を変換しない（`.gitattributes` の `-text`）。

`intensities_stage` の位置は見本でも manifest の最上位で、この文書と一致する。

## 統計エンジンの結果のバンドル（ms-data-parser 2026-11 の構成見直し以降）

ms-data-parser main 4a4eee3 から、`dataset_export_evidence_bundle` が書くのは統計エンジン（`dataset_compare`）の結果だけになった
（旧 `arf_export_evidence_bundle` は無い）。契約の版は 3 のままで、massbank-context の読み取りは変えていない。
manifest の**値の意味**は次のように変わった（見本は 4a4eee3 の写しに更新済み）。

- `comparison.method` は実際の検定名（既定 `moderated_t`。以前は常に `Welch`）。massbank-context はこの値を見ない。
- `comparison` に `statistics`（検定・CI 水準・判定規則・BH の母集団・比較の単位）・`contrast_id`・`matrix_id`・`rule_hash`・
  `effect_size_definition`・`mean_definition` が増えた（未知の欄として素通しで読む）。完全交絡を越えた比較なら `unadjusted_confounded: true`。
- features の `significant` は q・|log2fc| の閾値に加え `verdict == "changed"` のときだけ真（`inconclusive` は有意にしない）。
- `mean_a` / `mean_b` は単位（被験者）ごとの前処理後の線形値の算術平均。`n_a` / `n_b` は注入の数のまま。
- `population.preprocessing_report` は **null**、`derivation_status` は `not_available`、`n_features_before` は null。
  そのため `get_comparison_evidence` の `population_derivation.caveats`（前処理の注意。例: プール QC の層別）は**常に空**になり、
  同じ所の `derivation_status` が `not_available` になる。「注意が無い」ではなく「上流が渡していない」。
  前処理の設定は `comparison.preprocess`（行列の設定）にある。
- `population.exclusions` は `{stage: excluded_samples | excluded_features, removed, names, n_more}` の 2 行。
- `intensities_stage` は常に `pre_imputation`（結果が読んだ行列の正規化後・log2 前の値。エンジンは補完しない）。

## dataset 経路の spectrum_links

ms-data-parser 8631538（2026-11-05）からは、同じアラインメントの `.arf2` / `.dcl` が見つかれば `dataset_export_evidence_bundle` が
`spectrum_links` を書き、見つからなければ空で `provenance.spectrum_links_reason` に理由を書く。以下はそれ以前の出力についての記録。

ms-data-parser の `dataset_export_evidence_bundle` は v3 でも `spectrum_links` を空で書く
（`provenance.spectrum_links_status` が `not_provided`）。v3 の契約としては正しいが、外来の分子の走査では MS/MS の筋全体が
`unavailable`（`msms_all: no_spectrum_links`）になり、水準 A・B が付かない。`verified_alignment` 以外で links が空のときは同じ扱いで、
`verified_alignment` で空のときだけ `no_msms` と読む（設計書 §6）。r2 までの較正では同じアラインメントの `.dcl` と ID で結んだ links を足した
（`scripts/build_xenobiotic_calibration_bundles.py`、設計書 §10.2）。
