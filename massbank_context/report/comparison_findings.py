"""会話側LLMに事実・仮説・参照可能な根拠を渡す。"""
from collections import Counter
from massbank_context.evidence.selection import candidate_dicts
from massbank_context.report.comparison import evidence_page, evidence_rows, population_projection
from massbank_context.report.comparison_population import compute_pathways, taxid_key


def _excluded_features(bundle, exclude_feature_ids, include_pathways) -> frozenset:
    """経路集計の入力から除く特徴量。黙って無視しない（未知のIDと、効かない指定は誤り）。"""
    if not exclude_feature_ids:
        return frozenset()
    if not isinstance(exclude_feature_ids, list) or not all(isinstance(x, str) for x in exclude_feature_ids):
        raise ValueError('exclude_feature_idsは特徴量IDの配列で指定してください')
    if not include_pathways:
        raise ValueError('exclude_feature_idsは経路集計の入力から除く指定なので、include_pathways=Trueのときだけ使えます')
    unknown = sorted(set(exclude_feature_ids) - {f.feature_id for f in bundle.features})
    if unknown:
        raise ValueError(f'exclude_feature_idsに未知の特徴量があります: {unknown}')
    return frozenset(exclude_feature_ids)


def build_comparison_findings(store, analysis_id, organism='', include_pathways=False,
                              exclude_feature_ids=None):
    view = store.read_view(analysis_id)
    excluded = _excluded_features(view[0], exclude_feature_ids, include_pathways)
    page = evidence_page(view, analysis_id)
    context = page['context']
    confirmed = context.get('status') == 'confirmed'
    pathways = dict(tested=False, reason='not_requested')
    if include_pathways:
        species = context.get('context', {}).get('species', {})
        value = species.get('value')
        taxid = value.get('taxid') if isinstance(value, dict) else value
        if not confirmed or species.get('status') != 'confirmed' or not taxid:
            pathways = dict(tested=False, reason='organism_unconfirmed')
        elif taxid_key(taxid) != taxid_key(organism):
            pathways = dict(tested=False, reason='organism_context_mismatch')
        else:
            bundle, _, state = view
            matches = state.get('matches', {})
            projection = population_projection(bundle, candidate_dicts(bundle, matches), matches, excluded)
            pathways = compute_pathways(projection, state.get('lookups', {}).get('lookups', []), organism)
            pathways['exclusion_counts'] = dict(Counter(x['reason'] for x in projection.exclusions.values()))
            pathways.pop('exclusions', None)
            pathways.pop('representatives', None)
            all_rows = pathways.get('rows', [])
            pathways['n_rows_total'] = len(all_rows)
            pathways['rows_truncated'] = len(all_rows) > 10
            pathways['rows'] = all_rows[:10]
            for row in pathways['rows']:
                members = row['feature_ids']
                row.update(n_member_features=len(members), members_truncated=len(members) > 50,
                           feature_ids=members[:50])
    return dict(analysis_id=analysis_id, evidence_revision=page['evidence_revision'],
                context_revision=page['context_revision'], facts=page,
                hypotheses=[dict(slot='biological_interpretation', claim=None,
                                 feature_ids=[], candidate_ids=[], evidence_ids=[],
                                 candidate_dependency=None)], pathways=pathways,
                readiness=dict(context_confirmed=confirmed, population_complete=page['summary']['population_complete']),
                next_actions=([] if confirmed else ['研究者に背景を確認しset_analysis_contextで登録する'])
                             + (['get_comparison_evidenceで次ページを取得する'] if page['next_cursor'] else [])
                             + ['主張の引用IDをvalidate_comparison_claimsで検証する'],
                limits=['代謝物量の変化とパスウェイ参加関係は経路活性を直接証明しません。',
                        '候補順位はスペクトル由来のままです。既存注釈の評価値と合成しません。',
                        '未照会・ヒットなし・失敗・未提供は別の状態です。',
                        '非有意は同等性の証拠ではありません。',
                        'RDF情報はMassBankメタデータとPubChemパスウェイに限ります。',
                        '一般知識はmodel_knowledge、候補依存の説明は仮説として分けてください。'])


def reference_universe(view, analysis_id) -> dict:
    """ページに関係なく、現在の根拠版で引用できる全ID。"""
    bundle, revision, state = view
    candidates = candidate_dicts(bundle, state.get('matches', {}))
    return dict(analysis_id=analysis_id, evidence_revision=revision,
                context_revision=state.get('context', {}).get('context_revision'),
                features=[dict(feature_id=f.feature_id) for f in bundle.features],
                candidates=candidates, evidence=evidence_rows(candidates, state.get('lookups', {})))


def validate_claims(store, analysis_id, claims) -> dict:
    if not isinstance(claims, list) or not all(isinstance(c, dict) for c in claims):
        raise ValueError('claimsはオブジェクトの配列で指定してください')
    universe = reference_universe(store.read_view(analysis_id), analysis_id)
    issues = validate_claim_references(claims, universe)
    return dict(analysis_id=analysis_id, evidence_revision=universe['evidence_revision'],
                context_revision=universe['context_revision'], n_claims=len(claims),
                valid=not issues, issues=issues,
                note='IDの存在と対応だけを検証します。主張の科学的妥当性は研究者が評価してください。')


def validate_claim_references(claims, evidence):
    """ID存在と対応だけを検証する。主張の科学的正しさは保証しない。

    evidenceは全件（reference_universe）を渡す。1ページだけを渡すと、
    ページ外の正しいIDも未知として報告される。
    """
    page = evidence.get('facts', evidence)
    features = {f['feature_id'] for f in page['features']}
    candidates = {c['candidate_id']: c for c in page['candidates']}
    records = {e['evidence_id']: e for e in page['evidence']}
    issues = []
    for claim in claims:
        cid = claim.get('claim_id')
        fids = set(claim.get('feature_ids', []))
        cids = set(claim.get('candidate_ids', []))
        for eid in claim.get('evidence_ids', []):
            row = records.get(eid)
            if row is None:
                issues.append(dict(claim_id=cid, code='unknown_evidence', evidence_id=eid))
            elif row['candidate_id'] not in cids or row['feature_id'] not in fids:
                issues.append(dict(claim_id=cid, code='evidence_candidate_mismatch', evidence_id=eid))
        if fids - features:
            issues.append(dict(claim_id=cid, code='unknown_feature', feature_ids=sorted(fids - features)))
        for key in sorted(cids):
            if key not in candidates:
                issues.append(dict(claim_id=cid, code='unknown_candidate', candidate_id=key))
            elif candidates[key]['feature_id'] not in fids:
                issues.append(dict(claim_id=cid, code='candidate_feature_mismatch', candidate_id=key))
        if claim.get('context_revision') and claim['context_revision'] != page['context_revision']:
            issues.append(dict(claim_id=cid, code='context_revision_mismatch'))
        if claim.get('evidence_revision') and claim['evidence_revision'] != page['evidence_revision']:
            issues.append(dict(claim_id=cid, code='evidence_revision_mismatch'))
    return issues
