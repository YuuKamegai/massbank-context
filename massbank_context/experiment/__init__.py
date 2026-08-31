"""実験側の入力（ms-data-parser の契約）を読む層。"""

from massbank_context.experiment.contract import (
    ContractError,
    DifferentialSet,
    parse_differential,
)

__all__ = ["ContractError", "DifferentialSet", "parse_differential"]
