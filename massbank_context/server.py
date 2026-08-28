"""MCP サーバのファサード。

import 副作用で @mcp.tool を登録し、各モジュールの __all__ を再エクスポートする。
新しいツールを足すときは、実体を該当モジュールに書き __all__ に載せること。
"""
from __future__ import annotations

from mcp.server.fastmcp import FastMCP

from massbank_context.tools.index_tools import index_status as _index_status
from massbank_context.tools.match_tools import match_spectra as _match_spectra
from massbank_context.tools.spectra_tools import (
    inspect_spectrum as _inspect_spectrum,
    load_spectra as _load_spectra,
)

mcp = FastMCP("massbank-context")

# structured_output=False は必須。付けないと FastMCP が outputSchema を導出し、
# MCP が同じ内容を content と structuredContent の両方で送る（＝2 倍）。
mcp.tool(structured_output=False)(_index_status)
mcp.tool(structured_output=False)(_load_spectra)
mcp.tool(structured_output=False)(_inspect_spectrum)
mcp.tool(structured_output=False)(_match_spectra)


if __name__ == "__main__":
    mcp.run()
