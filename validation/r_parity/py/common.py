"""Shared helpers of the R parity harness: a result row and the elementary comparisons."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

TOLERANCE = 1e-12
"""Largest absolute difference still reported as ``equal``. Both sides compute in double
precision with different libraries and a different order of operations."""


@dataclass
class Row:
    """One compared quantity of one case."""

    case: str
    variant: str  # which R run the port is compared with: verbatim, patched, ...
    quantity: str
    n: int  # number of compared items
    n_mismatch: int  # items that differ beyond the tolerance, or names on one side only
    max_abs_diff: float
    status: str  # equal | differs | r_error | info
    detail: str = ""


def status_of(n_mismatch: int) -> str:
    return "equal" if n_mismatch == 0 else "differs"


def compare_vectors(
    case: str,
    variant: str,
    quantity: str,
    r_values: pd.Series,
    py_values: pd.Series,
    tolerance: float = TOLERANCE,
) -> Row:
    """Two named numeric vectors: same names, and values within the tolerance.

    A value missing on both sides counts as equal; missing on one side only is a mismatch.
    """
    r_names, py_names = set(r_values.index), set(py_values.index)
    only_r, only_py = r_names - py_names, py_names - r_names
    common = sorted(r_names & py_names)
    a = r_values.reindex(common).to_numpy(dtype=float)
    b = py_values.reindex(common).to_numpy(dtype=float)
    both_nan = np.isnan(a) & np.isnan(b)
    one_nan = np.isnan(a) ^ np.isnan(b)
    diff = np.where(both_nan | one_nan, 0.0, np.abs(a - b))
    over = (diff > tolerance) | one_nan
    n_mismatch = int(over.sum()) + len(only_r) + len(only_py)
    parts = []
    if only_r:
        parts.append(f"{len(only_r)} names only in R, e.g. {sorted(only_r)[:3]}")
    if only_py:
        parts.append(f"{len(only_py)} names only in Python, e.g. {sorted(only_py)[:3]}")
    if one_nan.any():
        parts.append(f"{int(one_nan.sum())} values missing on one side only")
    if over.any() and not one_nan.all():
        worst = int(np.argmax(diff))
        parts.append(
            f"largest difference at {common[worst]}: R {a[worst]:.6g}, Python {b[worst]:.6g}"
        )
    return Row(
        case,
        variant,
        quantity,
        len(common) + len(only_r) + len(only_py),
        n_mismatch,
        float(diff.max()) if len(diff) else 0.0,
        status_of(n_mismatch),
        "; ".join(parts),
    )


def compare_sets(
    case: str, variant: str, quantity: str, r_items: Iterable[str], py_items: Iterable[str]
) -> Row:
    r_set, py_set = {str(i) for i in r_items}, {str(i) for i in py_items}
    only_r, only_py = r_set - py_set, py_set - r_set
    parts = [f"R {len(r_set)}, Python {len(py_set)}, shared {len(r_set & py_set)}"]
    if only_r:
        parts.append(f"only in R e.g. {sorted(only_r)[:3]}")
    if only_py:
        parts.append(f"only in Python e.g. {sorted(only_py)[:3]}")
    n_mismatch = len(only_r) + len(only_py)
    return Row(
        case,
        variant,
        quantity,
        len(r_set | py_set),
        n_mismatch,
        0.0,
        status_of(n_mismatch),
        "; ".join(parts),
    )


def compare_scalars(
    case: str, variant: str, quantity: str, r_value: object, py_value: object
) -> Row:
    same = r_value == py_value
    return Row(
        case,
        variant,
        quantity,
        1,
        0 if same else 1,
        0.0,
        status_of(0 if same else 1),
        f"R {r_value}, Python {py_value}",
    )


def read_named(path: Path, column: str) -> pd.Series:
    """A two-column table written by the R side: names (the literal NA allowed) and values."""
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    values = pd.to_numeric(table[column].replace("NA", "nan"), errors="coerce")
    return pd.Series(values.to_numpy(dtype=float), index=table["gene"].astype(str))


def read_r_matrix(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    table = table.set_index("row")
    return table.replace("NA", "nan").apply(pd.to_numeric, errors="coerce")


def read_lines(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [line for line in path.read_text().splitlines() if line != ""]


def write_report(rows: list[Row], path: Path) -> pd.DataFrame:
    frame = pd.DataFrame([asdict(r) for r in rows])
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    return frame


def print_report(frame: pd.DataFrame) -> None:
    with pd.option_context(
        "display.width", 250, "display.max_rows", 2000, "display.max_colwidth", 90
    ):
        shown = frame.copy()
        shown["max_abs_diff"] = shown["max_abs_diff"].map(lambda v: f"{v:.1e}")
        print(shown.to_string(index=False))
