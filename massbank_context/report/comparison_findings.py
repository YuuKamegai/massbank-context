"""会話側LLMに事実・仮説・参照可能な根拠を渡す。"""
from massbank_context.report.comparison import get_evidence
from collections import Counter
from massbank_context.report.comparison_population import build_population, compute_pathways
from massbank_context.evidence.selection import candidate_dicts


def build_comparison_findings(store, analysis_id, organism='', include_pathways=False):
    page = get_evidence(store, analysis_id)
    context = page['context']
    confirmed = context.get('status') == 'confirmed'
    pathways = dict(tested=False, reason='not_requested')
    if include_pathways:
        species = context.get('context', {}).get('species', {})
        value = species.get('value')
        taxid = value.get('taxid') if isinstance(value, dict) else value
        if not confirmed or species.get('status') != 'confirmed' or not taxid:
            pathways = dict(tested=False, reason='organism_unconfirmed')
        elif str(taxid) != str(organism):
            pathways = dict(tested=False, reason='organism_context_mismatch')
        else:
            bundle, revision, state = store.read_view(analysis_id)
            if revision != page['evidence_revision']:
                raise ValueError('stale_cursor: 根拠版が変わりました')
            projection = build_population(bundle, candidate_dicts(bundle, state.get('matches', {})),
                                          state.get('matches', {}).get('features', {}))
            pathways = compute_pathways(projection, state.get('lookups', {}).get('lookups', []), organism)
            pathways['exclusion_counts'] = dict(Counter(x['reason'] for x in projection.exclusions.values()))
            pathways.pop('exclusions',None)
            pathways.pop('representatives',None)
            all_rows = pathways.get('rows', [])
            pathways['n_rows_total'] = len(all_rows)
            pathways['rows_truncated'] = len(all_rows)>10
            pathways['rows'] = all_rows[:10]
            for row in pathways['rows']:
                members = row['feature_ids']
                row.update(n_member_features=len(members),members_truncated=len(members)>50,
                           feature_ids=members[:50])
    return dict(analysis_id=analysis_id, evidence_revision=page['evidence_revision'],
                context_revision=page['context_revision'], facts=page,
                hypotheses=[dict(slot='biological_interpretation', claim=None,
                                 feature_ids=[], candidate_ids=[], evidence_ids=[],
                                 candidate_dependency=None)], pathways=pathways,
                readiness=dict(context_confirmed=confirmed, population_complete=page['summary']['population_complete']),
                next_actions=([] if confirmed else ['研究者に背景を確認しset_analysis_contextで登録する'])
                             + (['get_comparison_evidenceで次ページを取得する'] if page['next_cursor'] else []),
                limits=['代謝物量の変化とパスウェイ参加関係は経路活性を直接証明しません。',
                        '候補順位はスペクトル由来のままです。既存注釈の評価値と合成しません。',
                        '未照会・ヒットなし・失敗・未提供は別の状態です。',
                        '非有意は同等性の証拠ではありません。',
                        'RDF情報はMassBankメタデータとPubChemパスウェイに限ります。',
                        '一般知識はmodel_knowledge、候補依存の説明は仮説として分けてください。'])


def validate_claim_references(claims, evidence):
    """ID存在と対応だけを検証する。主張の科学的正しさは保証しない。"""
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
            issues.append(dict(claim_id=cid, code='unknown_feature'))
        for key in cids:
            if key not in candidates:
                issues.append(dict(claim_id=cid, code='unknown_candidate'))
            elif candidates[key]['feature_id'] not in fids:
                issues.append(dict(claim_id=cid, code='candidate_feature_mismatch'))
        if claim.get('context_revision') and claim['context_revision'] != page['context_revision']:
            issues.append(dict(claim_id=cid, code='context_revision_mismatch'))
    return issues
