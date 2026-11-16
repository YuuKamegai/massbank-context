"""較正の期待（tests/fixtures/xenobiotic/calibration_expectations.json）を機械的に作る（spec §8）。

走査を回す前に、report.md の表の行とバンドルの値だけから特徴 ID を決める。規則は
``SELECTION_RULE`` に書き、JSON にも写す。結果を見て変えない。

期待は結果を見て変えないので、`OUTPUT` が既にあれば `--force` なしでは止める。
r3 の前の改訂（設計書 §10.9）は JSON を直接直した。`--force` で作り直すと改訂が消える。

使い方: `.venv/Scripts/python.exe scripts/build_xenobiotic_calibration_expectations.py [--force]`
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from massbank_context.core.formula import ELECTRON_MASS, monoisotopic_mass, parse_formula  # noqa: E402

OUTPUT = REPO / 'tests' / 'fixtures' / 'xenobiotic' / 'calibration_expectations.json'
DATASETS = Path(r'C:\Users\yuu18\datasets\reanalysis')
STUDY = Path(r'C:\Users\yuu18\reanalysis-study\runs')

SELECTION_RULE = (
    '(1) report.md の表に行がある対象は、その行の元の results/novel_candidates.tsv'
    '（メトホルミンは results/claims_check.json の C4 の chosen）の spot_id を特徴 ID とし、'
    'バンドルの rt・log2fc が行の値と一致することを確かめる。'
    '(2) 行の無い本体・代謝物は、run の判断ログ 2026-10-11 の規則と同じく、理論イオンから ±5 ppm 以内で '
    'MS/MS のある特徴のうち、登録の InChIKey の先頭ブロックの注釈があればそれに絞り、'
    'mean_a + mean_b が最大のものを取る。無ければ「検出されていない」として skipped に書く。'
    '(3) spec §8 により、ST002901 の陰イオンの 2 件は expect に入れず skipped に書いて水準の記録だけにする。'
)
DECOY_RULE = (
    'おとりは A・B が 0。C の上限は極性ごとに 20 件（走査を回す前に決めた。spec §8 の「少数」を、'
    '手順書がレポートで印を付ける上位の件数 20 に読んだ）。'
)

# 判断ログ 2026-10-11 の理論イオン（add: 足す質量の組成、loss: 引く組成、n: 分子数、z: 価数）
IONS = {
    'positive': (('[M+H]+', 1, 1, 'H', ''), ('[M+Na]+', 1, 1, 'Na', ''), ('[M+NH4]+', 1, 1, 'NH4', ''),
                 ('[M+H-H2O]+', 1, 1, 'H', 'H2O'), ('[M+H-NH3]+', 1, 1, 'H', 'NH3'),
                 ('[2M+H]+', 2, 1, 'H', ''), ('[M+2H]2+', 1, 2, 'H2', '')),
    'negative': (('[M-H]-', 1, -1, '', 'H'), ('[M-H2O-H]-', 1, -1, '', 'H3O'), ('[M+Cl]-', 1, -1, 'Cl', ''),
                 ('[M+FA-H]-', 1, -1, 'CHO2', ''), ('[2M-H]-', 2, -1, '', 'H'),
                 ('[M+Na-2H]-', 1, -1, 'Na', 'H2')),
}

METFORMIN = dict(entry_id='X1', role='drug', name='metformin', formula='C4H11N5',
                 inchikey='XZWYZXLIPXDOLR-UHFFFAOYSA-N', received_by='test')
FLUOXETINE = dict(entry_id='X1', role='drug', name='fluoxetine', formula='C17H18F3NO',
                  inchikey='RTHCYVBBDHJXIQ-UHFFFAOYSA-N', received_by='test')
ST002901_ENTRIES = [
    FLUOXETINE,
    dict(entry_id='X2', role='drug_metabolite', name='norfluoxetine', formula='C16H16F3NO',
         inchikey='WIQRCHMSJFFONW-UHFFFAOYSA-N', parent_entry_id='X1', received_by='test'),
    dict(entry_id='X3', role='anesthetic', name='ketamine', formula='C13H16ClNO',
         inchikey='YQEZLKZALYSWHR-UHFFFAOYSA-N', received_by='both'),
    dict(entry_id='X4', role='drug_metabolite', name='norketamine', formula='C12H14ClNO',
         inchikey='BEQZHFIKTBVCAU-UHFFFAOYSA-N', parent_entry_id='X3', received_by='both'),
    dict(entry_id='X5', role='drug_metabolite', name='hydroxynorketamine', formula='C12H14ClNO2',
         parent_entry_id='X3', received_by='both'),
    dict(entry_id='X6', role='drug_metabolite', name='dehydronorketamine', formula='C12H12ClNO',
         parent_entry_id='X3', received_by='both'),
    dict(entry_id='X7', role='anesthetic', name='xylazine', formula='C12H16N2S', received_by='both'),
]
MTBLS9365_ENTRIES = [
    dict(entry_id='X1', role='drug', name='amlodipine', formula='C20H25ClN2O5',
         inchikey='HTIQEAQVCYTUBX-UHFFFAOYSA-N', received_by='test'),
    dict(FLUOXETINE, entry_id='X2'),
]
# run の判断ログ「質量精度と許容幅」の帯（陽陰とも <200 は ±10、≥200 は ±5）
MTBLS9365_BANDS = [{'min_mz': None, 'max_mz': 200.0, 'ppm': 10.0},
                   {'min_mz': 200.0, 'max_mz': None, 'ppm': 5.0}]

# (表示名, 組成式, 期待の水準)。report.md の新しい知見の候補の表の行。
ST003841_REPORT = (('1-Methylbiguanide', 'C'), ('N-Isopropylimidodicarbonimidic diamide', 'C'),
                   ("N-Cyano-N',N'-dimethylguanidine", 'C'), ('6-Propyl-1,3,5-triazine-2,4-diamine', 'C'))
ST002901_REPORT_POS = (('(.+/-.)-Norketamine', 'B', 'norketamine'),
                       ('(2R,6R)-Hydroxynorketamine', 'B', 'hydroxynorketamine'),
                       ('Dehydronorketamine', 'B', 'dehydronorketamine'),
                       ('2-Ethyl-2-(p-tolyl)malonamide', 'C', None))
ST002901_REPORT_NEG = ('Autophinib', '4-Bromo-2-{[(1-propyl-1H-benzimidazol-2-yl)amino]methyl}phenol')


def _table(path: Path) -> list[dict]:
    with path.open(encoding='utf-8', newline='') as stream:
        lines = [line for line in stream if not line.startswith('#')]
    return list(csv.DictReader(lines, delimiter='\t'))


class Bundle:
    def __init__(self, accession: str, pol: str):
        self.dir = DATASETS / accession / 'xenobiotic_calibration' / pol / 'bundle'
        self.features = {row['feature_id']: row for row in _table(self.dir / 'features.tsv')}
        self.links = {row['feature_id'] for row in _table(self.dir / 'spectrum_links.tsv')}
        self.keys: dict[str, set[str]] = {}
        for row in _table(self.dir / 'annotations.tsv'):
            if row['inchikey']:
                self.keys.setdefault(row['feature_id'], set()).add(row['inchikey'].split('-')[0])


def _ion_mz(formula: str, n: int, z: int, add: str, loss: str) -> float:
    mass = n * monoisotopic_mass(parse_formula(formula))
    mass += monoisotopic_mass(parse_formula(add)) if add else 0.0
    mass -= monoisotopic_mass(parse_formula(loss)) if loss else 0.0
    return (mass - z * ELECTRON_MASS) / abs(z)


def _from_report_row(bundle: Bundle, row: dict, source: str) -> tuple[str, str]:
    """(特徴 ID, 根拠)。バンドルの値が表の行と食い違えば止める。"""
    feature = bundle.features[row['spot_id']]
    if abs(float(feature['rt']) - float(row['rt'])) > 1e-3 \
            or abs(float(feature['log2fc']) - float(row['log2fc'])) > 1e-4:
        raise SystemExit(f'バンドルの値が行と違う: {row}')
    return row['spot_id'], (f"{source}: {row.get('name')} spot {row['spot_id']}、"
                            f"{row.get('nearest_ion') or row.get('adduct')} {row.get('ppm_error')} ppm、"
                            f"RT {float(row['rt']):.3f} 分、log2FC {float(row['log2fc']):.2f}、"
                            f"バンドルの m/z {float(feature['mz']):.4f}")


def _by_rule(bundle: Bundle, entry: dict, polarity: str) -> tuple[str | None, str]:
    """規則 (2)。(特徴 ID か None, 根拠)。"""
    block = (entry.get('inchikey') or '').split('-')[0]
    found = []
    for fid, row in bundle.features.items():
        if not row['mz'] or fid not in bundle.links:
            continue
        for name, n, z, add, loss in IONS[polarity]:
            theory = _ion_mz(entry['formula'], n, z, add, loss)
            ppm = (float(row['mz']) - theory) / theory * 1e6
            if abs(ppm) <= 5:
                total = float(row['mean_a'] or 0) + float(row['mean_b'] or 0)
                found.append((block in bundle.keys.get(fid, ()), total, fid, name, ppm))
                break
    if not found:
        return None, f"理論イオン（{polarity}）±5 ppm に MS/MS のある特徴が無い"
    if any(item[0] for item in found):
        found = [item for item in found if item[0]]
    keyed, total, fid, name, ppm = max(found, key=lambda item: item[1])
    row = bundle.features[fid]
    return fid, (f"規則 (2): 候補 {len(found)} 件のうち mean_a+mean_b 最大。spot {fid}、{name} {ppm:+.2f} ppm、"
                 f"RT {float(row['rt']):.3f} 分、m/z {float(row['mz']):.4f}、"
                 f"{'注釈の InChIKey あり' if keyed else '注釈の InChIKey なし（質量と MS/MS の有無だけで選んだ）'}")


def _novel(accession: str) -> dict[tuple[str, str], dict]:
    rows = _table(STUDY / accession / 'results' / 'novel_candidates.tsv')
    return {(row['polarity'], row['name']): row for row in rows}


def st003841() -> dict:
    bundle = Bundle('ST003841', 'pos')
    claims = json.loads((STUDY / 'ST003841' / 'results' / 'claims_check.json').read_text(encoding='utf-8'))
    chosen = next(c for c in claims['claims'] if c['id'] == 'C4')['chosen']
    fid, basis = _from_report_row(bundle, chosen, 'report.md 1 節 C4 / claims_check.json')
    expect = [dict(feature_id=fid, label='metformin', min_level='A', basis=basis)]
    novel = _novel('ST003841')
    for name, level in ST003841_REPORT:
        fid, basis = _from_report_row(bundle, novel[('pos', name)], 'report.md 2 節 / novel_candidates.tsv')
        expect.append(dict(feature_id=fid, label=name, min_level=level, basis=basis))
    return dict(accession='ST003841', polarity='positive', bundle_dir=str(bundle.dir),
                test_group='VilCre_Met', ppm=10.0, entries=[METFORMIN], expect=expect, skipped=[])


def st002901() -> list[dict]:
    novel = _novel('ST002901')
    bundle = Bundle('ST002901', 'pos')
    expect, skipped = [], []
    for entry, level in ((ST002901_ENTRIES[0], 'A'), (ST002901_ENTRIES[1], 'B'), (ST002901_ENTRIES[2], 'A')):
        fid, basis = _by_rule(bundle, entry, 'positive')
        if fid is None:
            skipped.append(dict(label=entry['name'], reason=basis))
        else:
            expect.append(dict(feature_id=fid, label=entry['name'], min_level=level, basis=basis))
    for name, level, label in ST002901_REPORT_POS:
        fid, basis = _from_report_row(bundle, novel[('pos', name)], 'report.md 2 節 / novel_candidates.tsv')
        expect.append(dict(feature_id=fid, label=label or name, min_level=level, basis=basis))
    positive = dict(accession='ST002901', polarity='positive', bundle_dir=str(bundle.dir),
                    test_group='CSIS_Flx', ppm=10.0, entries=ST002901_ENTRIES, expect=expect,
                    skipped=skipped)

    bundle = Bundle('ST002901', 'neg')
    skipped = []
    for name in ST002901_REPORT_NEG:
        fid, basis = _from_report_row(bundle, novel[('neg', name)], 'report.md 2 節 / novel_candidates.tsv')
        skipped.append(dict(feature_id=fid, label=name, basis=basis,
                            reason='spec §8: 索引のフルオキセチンは陽イオンのレコードだけで、陰イオンでは MS/MS の共有の筋が '
                                   '使えないので合格の条件から外す。水準の記録だけにする'))
    # キシラジン（X7）は spec §8 の本体にも既知の代謝物にも入らないので、どの極性でも期待にしない。
    for entry in ST002901_ENTRIES[:6]:
        fid, basis = _by_rule(bundle, entry, 'negative')
        if fid is None:
            skipped.append(dict(label=entry['name'], reason=f'検出されていない: {basis}'))
        else:
            raise SystemExit(f'陰イオンの扱いを決めていない: {entry["name"]} {basis}')
    negative = dict(accession='ST002901', polarity='negative', bundle_dir=str(bundle.dir),
                    test_group='CSIS_Flx', ppm=10.0, entries=ST002901_ENTRIES, expect=[],
                    skipped=skipped)
    return [positive, negative]


def mtbls9365() -> list[dict]:
    return [dict(accession='MTBLS9365', polarity=polarity, bundle_dir=str(Bundle('MTBLS9365', pol).dir),
                 test_group='methionine_choline_deficient_MCD', ppm=10.0, ppm_bands=MTBLS9365_BANDS,
                 entries=MTBLS9365_ENTRIES, expect=[], skipped=[], max_A=0, max_B=0, max_C=20)
            for pol, polarity in (('pos', 'positive'), ('neg', 'negative'))]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--force', action='store_true', help='既にある期待を書き直す')
    args = parser.parse_args()
    if OUTPUT.exists() and not args.force:
        raise SystemExit(f'{OUTPUT} が既にあります。書き直すなら --force')
    data = {'selection_rule': SELECTION_RULE, 'decoy_rule': DECOY_RULE,
            'runs': [st003841(), *st002901(), *mtbls9365()]}
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(data, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    for run in data['runs']:
        print(run['accession'], run['polarity'])
        for item in run['expect']:
            print('  expect', item['feature_id'], item['label'], item['min_level'], '|', item['basis'])
        for item in run['skipped']:
            print('  skip', item.get('feature_id'), item['label'], '|', item['reason'][:70])


if __name__ == '__main__':
    main()
