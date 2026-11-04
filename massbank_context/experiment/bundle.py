"""全特徴量を保持したv2/v3根拠バンドルの検証。"""
import csv
import hashlib
import io
import json
import math
import re
from collections.abc import Sequence
from pathlib import Path
from massbank_context.core.evidence_ids import stable_id, canonical_json
from massbank_context.core.evidence_models import Candidate, EvidenceBundle, Feature, Sample, SpectrumLink
from massbank_context.query.loader import load_spectra_file

INCHIKEY_PATTERN = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')
FEATURE_COLUMNS = ['feature_id', 'spot_id', 'mz', 'rt', 'log2fc', 'p_value',
                   'q_value', 'mean_a', 'mean_b', 'significant']
ANNOTATION_COLUMNS = ['annotation_id', 'feature_id', 'source', 'name', 'inchikey',
                      'rank', 'provenance_json', 'curation_json', 'status']
LINK_COLUMNS = ['feature_id', 'spectrum_file', 'spectrum_id',
                'source_analysis_id', 'mapping_evidence']
SAMPLE_COLUMNS = ['sample_id', 'group', 'role', 'injection_order']
INTENSITY_COLUMNS = ['feature_id', 'sample_id', 'value', 'detected']
SUPPORTED_CONTRACT_VERSIONS = (2, 3)
POLARITIES = ('positive', 'negative')
INTENSITIES_STAGES = ('pre_imputation', 'post_imputation')


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


def _table_text(text, name, columns):
    reader = csv.DictReader(io.StringIO(text, newline=''), delimiter='\t')
    names = reader.fieldnames or []
    if len(names) != len(set(names)) or not set(columns) <= set(names):
        raise BundleError(f'必須列または列名が不正です: {name}')
    rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise BundleError('TSV列数が一致しません')
    return rows


def _table(path, columns):
    return _table_text(path.read_text(encoding='utf-8-sig'), path.name, columns)


def parse_samples(text: str) -> tuple[Sample, ...]:
    """samples.tsv を検証して読む。group は a / b / 空のみ。"""
    samples, seen = [], set()
    for row in _table_text(text, 'samples.tsv', SAMPLE_COLUMNS):
        sid = row['sample_id']
        if not sid or sid in seen:
            raise BundleError('試料IDが空または重複しています', 'sample_id')
        seen.add(sid)
        if row['group'] not in ('a', 'b', ''):
            raise BundleError('試料のgroupはa/b/空のいずれかです', 'group')
        try:
            order = int(row['injection_order']) if row['injection_order'] != '' else None
        except ValueError as exc:
            raise BundleError('注入順は整数または空です', 'injection_order') from exc
        samples.append(Sample(sid, row['group'], row['role'], order))
    return tuple(samples)


def parse_intensities(text: str, feature_ids: set[str], samples: Sequence[Sample]):
    """intensities.tsv を検証して読む。全特徴 x 群a/bの試料が過不足なく揃う契約。

    戻り値は feature_id -> sample_id -> (value, detected)。detected の空は None。
    groupが空の試料（QCなど）の行は契約違反として拒否する（黙って捨てない）。
    """
    by_id = {s.sample_id: s for s in samples}
    expected = {s.sample_id for s in samples if s.group}
    matrix: dict = {}
    for row in _table_text(text, 'intensities.tsv', INTENSITY_COLUMNS):
        fid, sid = row['feature_id'], row['sample_id']
        if fid not in feature_ids:
            raise BundleError('値の特徴量IDがfeaturesにありません', 'feature_id')
        if sid not in by_id:
            raise BundleError('値の試料IDがsamplesにありません', 'sample_id')
        if not by_id[sid].group:
            raise BundleError('比較に入らない試料の値は書けません', 'sample_id')
        if sid in matrix.setdefault(fid, {}):
            raise BundleError('特徴量と試料の組が重複しています', 'sample_id')
        if row['detected'] not in ('true', 'false', ''):
            raise BundleError('detectedはtrue/false/空のいずれかです', 'detected')
        detected = None if row['detected'] == '' else row['detected'] == 'true'
        matrix[fid][sid] = (_number(row['value'], 'value'), detected)
    if any(set(matrix.get(fid, {})) != expected for fid in feature_ids):
        raise BundleError('全特徴量と全比較試料の値が揃っていません')
    return matrix


def role_path(manifest: dict, role: str) -> str | None:
    """manifest の files から役割のパスを返す。無ければ None。"""
    return next((f['path'] for f in manifest.get('files', []) if f.get('role') == role), None)


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
    if manifest['contract_name'] != 'massbank-context-evidence' or type(manifest['contract_version']) is not int \
            or manifest['contract_version'] not in SUPPORTED_CONTRACT_VERSIONS:
        raise BundleError('対応する版はmassbank-context-evidence v2 / v3です')
    v3 = manifest['contract_version'] == 3
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
    if (v3 or 'polarity' in c) and c.get('polarity') not in POLARITIES:
        raise BundleError('極性はpositive/negativeで指定してください', 'polarity')
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
        elif (role in ('features', 'annotations', 'spectrum_links')
              or (v3 and role in ('samples', 'intensities'))) and role not in roles:
            roles[role] = file
        else:
            raise BundleError('ファイル種別が不正または重複しています')
    if not {'features', 'annotations', 'spectrum_links'} <= set(roles):
        raise BundleError('必須TSVファイルがありません')
    if ('samples' in roles) != ('intensities' in roles):
        raise BundleError('samplesとintensitiesは両方そろえるか両方省いてください')
    if 'intensities' in roles and manifest.get('intensities_stage') not in INTENSITIES_STAGES:
        raise BundleError('intensities_stageはpre_imputation/post_imputationで指定してください')
    aid = stable_id('analysis', [canonical_json(manifest)])
    features, ids = [], set()
    for row in _table(roles['features'], FEATURE_COLUMNS):
        fid = row['feature_id']
        if not fid or fid in ids or row['significant'] not in ('true', 'false'):
            raise BundleError('特徴量IDまたは有意判定が不正です')
        ids.add(fid)
        numbers = {key: _number(row[key], key) for key in FEATURE_COLUMNS[2:-1]}
        features.append(Feature(fid, row['spot_id'], **numbers, significant=row['significant'] == 'true',
                                adduct=row.get('adduct') or None))
    population = manifest['population']
    if type(population['n_features_total']) is not int or population['n_features_total'] != len(features):
        raise BundleError('特徴量の総数と実際の行数が不一致です')
    if type(population['population_complete']) is not bool:
        raise BundleError('母集団の完全性は真偽値で指定してください')
    if 'samples' in roles:
        # 行列は保持しない。検証のためだけに読む（原本は store.read_source で読める）。
        samples = parse_samples(roles['samples'].read_text(encoding='utf-8-sig'))
        for key, group in (('n_a', 'a'), ('n_b', 'b')):
            if sum(x.group == group for x in samples) != c[key]:
                raise BundleError(f'試料数と{key}が一致しません', key)
        parse_intensities(roles['intensities'].read_text(encoding='utf-8-sig'), ids, samples)
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
