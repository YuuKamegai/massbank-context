"""全特徴量を保持したv2根拠バンドルの検証。"""
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from massbank_context.core.evidence_ids import stable_id, canonical_json
from massbank_context.core.evidence_models import Candidate, EvidenceBundle, Feature, SpectrumLink
from massbank_context.query.loader import load_spectra_file

INCHIKEY_PATTERN = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')
FEATURE_COLUMNS = ['feature_id', 'spot_id', 'mz', 'rt', 'log2fc', 'p_value',
                   'q_value', 'mean_a', 'mean_b', 'significant']
ANNOTATION_COLUMNS = ['annotation_id', 'feature_id', 'source', 'name', 'inchikey',
                      'rank', 'provenance_json', 'curation_json', 'status']
LINK_COLUMNS = ['feature_id', 'spectrum_file', 'spectrum_id',
                'source_analysis_id', 'mapping_evidence']


class BundleError(ValueError):
    def __init__(self, reason, field='', code='contract_mismatch'):
        super().__init__(reason)
        self.reason, self.field, self.code = reason, field, code


def read_json(text):
    def unique_pairs(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise BundleError(f'JSONのキーが重複しています: {key}')
            result[key] = value
        return result
    return json.loads(text, object_pairs_hook=unique_pairs,
                      parse_constant=lambda x: (_ for _ in ()).throw(BundleError(f'不正な数値: {x}')))


def _number(value, column):
    if value is None or value == '':
        return None
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise BundleError(f'数値ではありません: {column}', column) from exc
    if not math.isfinite(number) or (column in ('p_value', 'q_value') and not 0 <= number <= 1):
        raise BundleError(f'数値の範囲が不正です: {column}', column)
    return number


def _table(path, columns):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        names = reader.fieldnames or []
        if len(names) != len(set(names)) or not set(columns) <= set(names):
            raise BundleError(f'必須列または列名が不正です: {path.name}')
        rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise BundleError('TSV列数が一致しません')
    return rows


def _required(data, keys):
    if not isinstance(data, dict) or any(key not in data for key in keys):
        raise BundleError(f'必須項目がありません: {keys}')


def load_bundle(path: Path) -> EvidenceBundle:
    try:
        return _load_bundle(Path(path).resolve())
    except BundleError:
        raise
    except (OSError, ValueError, TypeError, KeyError, AttributeError) as exc:
        raise BundleError(f'バンドルを読み込めません: {exc}') from exc


def _load_bundle(path):
    manifest = read_json(path.read_text(encoding='utf-8-sig'))
    _required(manifest, ['contract_name', 'contract_version', 'source_analysis_id',
                         'exporter', 'exported_at', 'comparison', 'population', 'files'])
    if manifest['contract_name'] != 'massbank-context-evidence' or manifest['contract_version'] != 2:
        raise BundleError('対応する版はmassbank-context-evidence v2です')
    if not isinstance(manifest['source_analysis_id'], str) or not manifest['source_analysis_id'].strip():
        raise BundleError('出力元解析IDがありません')
    c = manifest['comparison']
    _required(c, ['group_a', 'group_b', 'n_a', 'n_b', 'log2fc_sign', 'method',
                  'correction', 'q_threshold', 'log2fc_threshold', 'preprocess'])
    if not c['group_a'] or not c['group_b'] or c['group_a'] == c['group_b']:
        raise BundleError('異なる2群が必要です')
    if c['log2fc_sign'] != 'positive means group_b is higher':
        raise BundleError('比較方向が不一致です')
    if any(type(c[k]) is not int or c[k] < 1 for k in ('n_a', 'n_b')):
        raise BundleError('群別nは正の整数で指定してください')
    if not 0 <= _number(c['q_threshold'], 'q_value') <= 1 or _number(c['log2fc_threshold'], 'threshold') < 0:
        raise BundleError('閾値が不正です')
    _required(manifest['exporter'], ['name', 'version'])
    _required(manifest['population'], ['n_features_total', 'population_complete', 'definition', 'exclusions'])
    root, roles, spectra, seen_files = path.parent, {}, {}, set()
    for item in manifest['files']:
        _required(item, ['role', 'path', 'sha256'])
        raw = Path(item['path'])
        file = (root / raw).resolve()
        if raw.is_absolute() or not file.is_relative_to(root) or file == path or file in seen_files:
            raise BundleError('参照ファイルが範囲外または重複しています')
        seen_files.add(file)
        if hashlib.sha256(file.read_bytes()).hexdigest() != item['sha256']:
            raise BundleError('ファイルハッシュが一致しません', item['path'])
        role = item['role']
        if role == 'spectra':
            for spectrum in load_spectra_file(file):
                key = item['path'] + '::' + spectrum.spectrum_id
                if key in spectra:
                    raise BundleError('スペクトルIDが重複しています')
                spectra[key] = spectrum
        elif role in ('features', 'annotations', 'spectrum_links') and role not in roles:
            roles[role] = file
        else:
            raise BundleError('ファイル種別が不正または重複しています')
    if set(roles) != {'features', 'annotations', 'spectrum_links'}:
        raise BundleError('必須TSVファイルがありません')
    aid = stable_id('analysis', [canonical_json(manifest)])
    features, ids = [], set()
    for row in _table(roles['features'], FEATURE_COLUMNS):
        fid = row['feature_id']
        if not fid or fid in ids or row['significant'] not in ('true', 'false'):
            raise BundleError('特徴量IDまたは有意判定が不正です')
        ids.add(fid)
        numbers = {key: _number(row[key], key) for key in FEATURE_COLUMNS[2:-1]}
        features.append(Feature(fid, row['spot_id'], **numbers, significant=row['significant'] == 'true'))
    population = manifest['population']
    if type(population['n_features_total']) is not int or population['n_features_total'] != len(features):
        raise BundleError('特徴量の総数と実際の行数が不一致です')
    if type(population['population_complete']) is not bool:
        raise BundleError('母集団の完全性は真偽値で指定してください')
    annotations, annotation_ids = [], set()
    for row in _table(roles['annotations'], ANNOTATION_COLUMNS):
        if row['feature_id'] not in ids or not row['annotation_id'] or row['annotation_id'] in annotation_ids:
            raise BundleError('注釈の参照またはIDが不正です')
        annotation_ids.add(row['annotation_id'])
        ik = row['inchikey'] or None
        if ik and not INCHIKEY_PATTERN.fullmatch(ik):
            raise BundleError('InChIKeyの形式が不正です')
        if row['status'] not in ('active', 'rejected', 'superseded') or not row['source']:
            raise BundleError('候補の状態または出典が不正です')
        rank = int(row['rank']) if row['rank'] else None
        if rank is not None and rank < 1:
            raise BundleError('順位は1以上です')
        provenance, curation = read_json(row['provenance_json']), read_json(row['curation_json'])
        structure_ids = read_json(row.get('structure_ids_json') or '{}')
        scores = read_json(row.get('scores_json') or '{}')
        if not all(isinstance(x, dict) for x in (provenance, curation, structure_ids, scores)):
            raise BundleError('候補メタデータはJSONオブジェクトです')
        annotations.append(Candidate(stable_id('candidate', [aid, row['annotation_id']]),
                                     row['feature_id'], row['source'], row['annotation_id'],
                                     name=row['name'], inchikey=ik, rank=rank, status=row['status'],
                                     provenance=provenance, curation=curation,
                                     structure_ids=structure_ids, scores=scores))
    links, seen_links = [], set()
    for row in _table(roles['spectrum_links'], LINK_COLUMNS):
        key = row['spectrum_file'] + '::' + row['spectrum_id']
        identity = (row['feature_id'], key)
        if (row['feature_id'] not in ids or key not in spectra or identity in seen_links
                or row['source_analysis_id'] != manifest['source_analysis_id'] or not row['mapping_evidence']):
            raise BundleError('スペクトルとの対応が不正です')
        seen_links.add(identity)
        links.append(SpectrumLink(**{k: row[k] for k in LINK_COLUMNS}))
    return EvidenceBundle(aid, manifest, tuple(features), tuple(annotations), tuple(links),
                          spectra, population['population_complete'], path)
