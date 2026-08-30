"""ClassyFire の化合物クラスを読む（spec §6.6 の層別軸）。

出所は MassBank-RDF リポジトリの `script/massbank_inchikey_classyfire.tsv`。
ヘッダ無しの 6 列で、レコード URI / InChIKey / kingdom / superclass / class /
subclass が並ぶ。行はレコード単位なので、同一 InChIKey が装置・CE 違いの
レコード数だけ重複して現れる。

`mbo:ch_compound_class` は使えるレコードと "N/A; Environmental Standard" のような
使えないレコードが混在するため（spec §6.6）、層別の主軸にはこちらを使う。
両者は `compound_class.source` で区別され、共存できる。

階層の 4 レベルは source に持たせる（`classyfire:superclass` など）。
`compound_class(inchikey, source, label)` の DDL を変えずにレベルを区別でき、
層別背景がどの粒度のクラスを使ったかを後から再現できる。

依存グラフの leaf（stdlib のみ）。
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterator

#: TSV の 3 列目以降に並ぶ階層。順序は列順そのもの。
CLASSYFIRE_LEVELS = ("kingdom", "superclass", "class", "subclass")

SOURCE_PREFIX = "classyfire"

_INCHIKEY_COLUMN = 1
_FIRST_LEVEL_COLUMN = 2

#: ラベルとして扱わない値。ClassyFire が分類できなかったことを表すもので、
#: これを層別の軸に載せると「N/A クラス」という偽のクラスができる。
#: 比較は writer._USELESS_CLASS_LABELS と同じく大小文字非依存で行う。
_NOT_A_LABEL = frozenset({"", "N/A", "NA", "NULL", "-"})


def source_for(level: str) -> str:
    """階層名を compound_class.source の値にする。"""
    return f"{SOURCE_PREFIX}:{level}"


def read_classyfire_tsv(path: str | Path) -> Iterator[tuple[str, str, str]]:
    """(inchikey, source, label) を重複無しで返す。

    重複除去はここで行う。TSV は 54,183 レコード分の行を持つ一方、
    InChIKey の異なり数はその数分の一なので、畳まずに INSERT に流すと
    無意味な書き込みが数十万回走る。
    """
    seen: set[tuple[str, str, str]] = set()
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            columns = line.rstrip("\n").split("\t")
            if len(columns) <= _FIRST_LEVEL_COLUMN:
                continue
            inchikey = columns[_INCHIKEY_COLUMN].strip()
            if not inchikey:
                continue
            for offset, level in enumerate(CLASSYFIRE_LEVELS):
                index = _FIRST_LEVEL_COLUMN + offset
                if index >= len(columns):
                    break
                label = columns[index].strip()
                if label.upper() in _NOT_A_LABEL:
                    continue
                entry = (inchikey, source_for(level), label)
                if entry not in seen:
                    seen.add(entry)
                    yield entry
