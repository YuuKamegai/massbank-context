"""候補根拠と2群比較の公開層。解釈は会話側LLMが行う。"""
import sqlite3
from collections import Counter
from pathlib import Path
from massbank_context.core import config
from massbank_context.core.serialization import json_payload, round_floats
from massbank_context.evidence.store import EvidenceStore
from massbank_context.evidence.context import set_context
from massbank_context.evidence.matching import match_features
from massbank_context.experiment.bundle import load_bundle
from massbank_context.experiment.legacy_bundle import import_legacy
from massbank_context.kg.evidence import resolve_candidates
from massbank_context.report.comparison import get_evidence, tabulate
from massbank_context.report.comparison_findings import build_comparison_findings, validate_claims

__all__ = ['load_evidence_bundle', 'set_analysis_context', 'match_feature_spectra',
           'resolve_candidate_evidence', 'get_comparison_evidence', 'draft_comparison_findings',
           'validate_comparison_claims']


def _invoke(operation):
    try:
        # 小さなp値などを0へ丸めない。旧ツールの書式は変えない。
        return json_payload(round_floats(operation(EvidenceStore(config.ARTIFACTS_PATH)), keep_small=True))
    except (ValueError, OSError, sqlite3.Error) as exc:
        code = getattr(exc, 'code', 'invalid_evidence_request')
        for known in ('missing_index', 'incompatible_index', 'stale_cursor', 'invalid_cursor', 'stale_revision'):
            if str(exc).startswith(known):
                code = known
        return json_payload({'error': {'code': code, 'message': str(exc),
                                       'required_action': '入力・解析ID・根拠版を確認し、必要な取得を再実行してください。'}})


def load_evidence_bundle(path: str) -> str:
    """ms-data-parserのv2/v3 manifestまたは旧v1 TSVを読み込む。旧形式は母集団不完全として扱う。"""
    def operation(store):
        bundle = import_legacy(Path(path)) if Path(path).suffix.lower() == '.tsv' else load_bundle(Path(path))
        aid = store.register(bundle)
        return dict(analysis_id=aid, n_features=len(bundle.features),
                    population_complete=bundle.population_complete,
                    next_tools=['set_analysis_context', 'get_comparison_evidence'])
    return _invoke(operation)


def set_analysis_context(analysis_id: str, context: dict, confirmed: bool = False) -> str:
    """会話から抽出した背景を登録する。研究者の確認後だけconfirmed=Trueにする。不明はnull。"""
    return _invoke(lambda store: set_context(store, analysis_id, context, confirmed))


def match_feature_spectra(analysis_id: str, feature_ids: list[str] | None = None,
                          top_k: int = 5, acquisition_class: str = 'LC-ESI-MS2',
                          precursor_tolerance: float = .01) -> str:
    """対応済みMS/MSをローカル照合。順位を背景で変更しない。MS/MSなし・候補打切りを報告する。"""
    def operation(store):
        result = match_features(store, analysis_id, feature_ids,top_k=top_k,
                                acquisition_class=acquisition_class, precursor_tolerance=precursor_tolerance)
        features = result.pop('features')
        result['status_counts'] = dict(Counter(v['status'] for v in features.values()))
        result['n_truncated_features'] = sum(v.get('truncated',False) for v in features.values())
        result['details_tool'] = 'get_comparison_evidence'
        return result
    return _invoke(operation)


def resolve_candidate_evidence(analysis_id: str, candidate_ids: list[str] | None = None,
                               refresh: bool = False) -> str:
    """候補のRDF関連を取得。構造IDだけをPubChemへ送信し、empty/failed/未照会を区別する。"""
    return _invoke(lambda store: resolve_candidates(store, analysis_id, candidate_ids, refresh=refresh))


def get_comparison_evidence(analysis_id: str, feature_ids: list[str] | None = None,
                            cursor: str | None = None, limit: int = 50) -> str:
    """全件概要と特徴量別根拠を取得。最大200特徴量・64根拠。features/candidates/evidence/lookupsは列名1回のTSV（入れ子はJSON）。next_cursorとnext_evidence_cursorをcursorへ渡す。版変更時は取り直す。"""
    return _invoke(lambda store: tabulate(get_evidence(store, analysis_id, feature_ids=feature_ids,
                                                       cursor=cursor, limit=limit)))


def draft_comparison_findings(analysis_id: str, organism: str = '',
                              include_pathways: bool = False) -> str:
    """解釈用の事実と空の仮説を返す。経路集計は任意。活性断定をせず根拠IDを引用して会話側で解釈する。"""
    def operation(store):
        result = build_comparison_findings(store, analysis_id, organism, include_pathways)
        return dict(result, facts=tabulate(result['facts']))
    return _invoke(operation)


def validate_comparison_claims(analysis_id: str, claims: list[dict]) -> str:
    """回答の各主張(claim_id, feature_ids, candidate_ids, evidence_ids, context_revision, evidence_revision)の引用IDが現行根拠版に実在し対応するか全件で検証する。意味の妥当性は保証しない。"""
    return _invoke(lambda store: validate_claims(store, analysis_id, claims))
