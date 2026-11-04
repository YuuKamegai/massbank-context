"""実stdio初期化と再起動を通して解析を復元する。"""
import asyncio
import json
import os
import sys
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from evidence_helpers import make_bundle


def test_stdio_restart_keeps_analysis(tmp_path):
    path = make_bundle(tmp_path / 'input')

    async def call(aid=None):
        params = StdioServerParameters(command=sys.executable,
                                       args=['-m', 'massbank_context.server'],
                                       env={**os.environ, 'MASSBANK_CONTEXT_ARTIFACTS': str(tmp_path / 'store')})
        async with stdio_client(params) as (read, write):
            async with ClientSession(read, write) as session:
                await session.initialize()
                tools = await session.list_tools()
                assert len(tools.tools) == 18
                if aid is None:
                    result = await session.call_tool('load_evidence_bundle', {'path': str(path)})
                else:
                    result = await session.call_tool('get_comparison_evidence', {'analysis_id': aid})
                assert not result.isError
                return json.loads(result.content[0].text)

    loaded = asyncio.run(call())
    page = asyncio.run(call(loaded['analysis_id']))
    assert page['summary']['n_features_total'] == 8
