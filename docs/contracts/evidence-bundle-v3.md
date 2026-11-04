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

見本のバンドルは ms-data-parser の `tests/fixtures/evidence_bundle_v3/` と共有する（写しは本リポジトリの
tests に置く。Task 9）。両リポのテストが同じ見本を読むことで、契約の食い違いを検出する。
