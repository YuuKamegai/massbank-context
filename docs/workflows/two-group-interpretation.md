# 会話側LLMで2群比較を解釈する

1. 研究者との会話から生物種・組織・処置・群の意味・測定条件・問いを整理する。
   値のない項目はnull。研究者に内容を示して確認し、確認済み背景を登録する。
   speciesのvalueは、経路集計する場合は `{"label":"mouse","taxid":"10090"}` の形を使う。
2. ms-data-parserで前処理・2群比較を実行し、新しいevidence bundleを出力する。
   正のlog2fcはgroup_bが高い。比較方向と群の意味を確認する。
3. `load_evidence_bundle`で読み、analysis_idと全件数を確認する。
4. `get_comparison_evidence`で概要を確認する。非有意・未同定も母集団に含まれる。
5. `match_feature_spectra`で対応済みMS/MSを照合する。MS/MSなしの特徴量も残る。
6. `resolve_candidate_evidence`で候補別のRDF根拠を取得する。構造情報のない名前だけの候補は照会不能。
7. `get_comparison_evidence`で根拠版を固定して全必要ページを取得する。
   特徴量ページはnext_cursor、RDF等の詳細はnext_evidence_cursorを同じツールのcursorへ渡す。
   各詳細ページは最大64根拠。lookupsは状態と件数を返し、巨大なhitsを二重送信しない。
   特徴量ページは約128KiBを目安に件数を減らすため、要求limitより少ない行でnext_cursorが返る場合がある。
   cursorが失効したら初めから取得し直す。大量の根拠は特徴量IDを絞って段階的に読む。
8. `draft_comparison_findings`で事実と解釈の材料を取得する。経路検定は任意。

回答は「観測事実」「候補と根拠」「解釈仮説」「限界と追加検証」に分ける。
例: 「f1はBで増加。候補c1を仮定すると、RDF根拠e1の経路参加関係と整合する。
候補c2ならこの説明は成立しない。量変化だけで活性化とは断定できない。」

各主張にfeature_ids/candidate_ids/evidence_ids/context_revisionを対応させ、実在するIDを引用する。
一般知識はmodel_knowledgeと記す。RDF未ヒットを試料に存在しない証拠にしない。
候補順位を背景で変更せず、候補共通の解釈と候補依存の仮説を区別する。
根拠の参照整合性は`validate_claim_references`で確認できるが、意味の妥当性は研究者が評価する。

背景の確認はサーバが本人確認する仕組みではなく、会話側が研究者の返答を受けて記録する手続き。
背景変更後は旧背景に基づく回答を新背景での評価として提示しない。
