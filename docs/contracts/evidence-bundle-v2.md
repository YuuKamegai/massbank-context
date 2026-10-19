# 根拠バンドル v2

入口はmanifest.json。contract_nameはmassbank-context-evidence、contract_versionは整数2。
source_analysis_id/exporter（name,version）/exported_at/comparison/population/filesが必須。
comparisonはgroup_a/group_b/n_a/n_b/log2fc_sign/method/correction/q_threshold/log2fc_threshold/preprocess。
方向は `positive means group_b is higher` に固定する。
populationはn_features_total/population_complete/definition/exclusionsを持つ。
filesはrole/path/sha256。roleはfeatures/annotations/spectrum_linksが各1個、spectraは任意個。
pathはバンドル内の相対パス。TSVはUTF-8（BOM可）、csv規則でタブや改行を引用できる。

- features: feature_id,spot_id,mz,rt,log2fc,p_value,q_value,mean_a,mean_b,significant
- annotations: annotation_id,feature_id,source,name,inchikey,rank,provenance_json,curation_json,status
- 注釈の任意列: structure_ids_json,scores_json
- spectrum_links: feature_id,spectrum_file,spectrum_id,source_analysis_id,mapping_evidence

数値欠測は空欄、有意判定はtrue/false。構造IDは空欄可。状態はactive/rejected/superseded。
リンクのspectrum_idは既存ローダのIDと一致させる。名前やm/zで推測して接続しない。
候補なし・MS/MSなしでもfeatures行を残し、対応TSVはヘッダのみで表現する。
旧v1は欠落行を含むため、新しい経路では母集団不完全として扱う。
