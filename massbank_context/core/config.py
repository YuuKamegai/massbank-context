"""パス解決と環境変数（すべて任意・上書き用）。

可変状態の正準はこのモジュール。参照側は `config.INDEX_PATH` のように
module 修飾で読むこと。`from ... import INDEX_PATH` はスナップショット束縛に
なり、テストでの差し替えが効かなくなる。
"""
from __future__ import annotations

import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent.parent

INDEX_PATH = Path(os.environ.get(
    "MASSBANK_CONTEXT_INDEX", BASE_DIR / "index" / "massbank.sqlite"))

BUILD_COMMAND = ("CLI: C:/Python314/python.exe -m massbank_context.index build"
                 " --ttl-dir <dir> --output <sqlite>")
