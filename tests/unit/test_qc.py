import numpy as np
import pandas as pd

from search.data import ExpressionSet
from search.qc import format_warnings, qc_expression


def codes(expr: ExpressionSet) -> set[str]:
    return {w.code for w in qc_expression(expr, "T")}


def test_warnings_from_names_and_values() -> None:
    samples = ["GB19.1", "GB19.2", "GbExp47_ctrl", "X1", "X1"]
    values = np.random.default_rng(0).normal(size=(4, 5))
    values[1] = 2.0
    values[2, 0] = np.nan
    expr = ExpressionSet(values, pd.Index(list("abcd")), pd.Index(samples))
    found = codes(expr)
    assert {
        "sample_table_missing",
        "missing_values",
        "constant_genes",
        "duplicate_samples",
        "control_like_names",
        "replicate_like_names",
    } <= found
    assert "QC warnings" in format_warnings(qc_expression(expr))


def test_warnings_from_sample_table() -> None:
    values = np.random.default_rng(1).normal(size=(3, 6))
    samples = [f"s{i}" for i in range(6)]
    table = pd.DataFrame(
        {
            "patient": ["p1", "p1", "p2", "p3", "p4", "p5"],
            "batch": ["a"] * 5 + ["b"],
            "tissue": ["tumour"] * 5 + ["normal"],
            "age": [50, np.nan, 60, 70, 80, 90],
        },
        index=samples,
    )
    expr = ExpressionSet(values, pd.Index(list("xyz")), pd.Index(samples), sample_table=table)
    found = codes(expr)
    assert {
        "repeated_patients",
        "unbalanced_batches",
        "non_tumour_tissue",
        "missing_covariate_values",
    } <= found
    assert "sample_table_missing" not in found


def test_clean_data_has_no_warnings() -> None:
    values = np.random.default_rng(2).normal(size=(3, 6))
    samples = ["a1", "b2", "c3", "d4", "e5", "f6"]
    table = pd.DataFrame({"batch": ["x", "x", "x", "y", "y", "y"]}, index=samples)
    expr = ExpressionSet(values, pd.Index(list("xyz")), pd.Index(samples), sample_table=table)
    assert qc_expression(expr) == []
    assert format_warnings([]) == "QC: no warnings"
