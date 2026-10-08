import numpy as np
import pandas as pd
import pytest

from search.combine import combine_cohorts, legacy_two_cohort_names
from search.config import CombineConfig, ScoreConfig
from search.data import ExpressionSet
from search.errors import DataError
from search.signature import median_abs_rank01, score_samples


def cohort_table(name: str, rows: list[dict]) -> pd.DataFrame:
    frame = pd.DataFrame(rows)
    frame.insert(0, "dataset", name)
    return frame


def test_combination_logic_and_ordering() -> None:
    a = cohort_table(
        "A",
        [
            {
                "gene": "seed",
                "C_raw": 1,
                "C_rank": 1.0,
                "C_rank_boot_mean": 1.0,
                "C_rank_boot_var": 0.0,
                "C_rank_boot_sd": 0.0,
                "p_MC_C": 1e-4,
                "positive_consistent": True,
                "weighted_seed_rho": 0.9,
                "is_seed": True,
            },
            {
                "gene": "good",
                "C_raw": 0.5,
                "C_rank": 0.9,
                "C_rank_boot_mean": 0.9,
                "C_rank_boot_var": 0.01,
                "C_rank_boot_sd": 0.1,
                "p_MC_C": 1e-3,
                "positive_consistent": True,
                "weighted_seed_rho": 0.6,
                "is_seed": False,
            },
            {
                "gene": "neg",
                "C_raw": 0.5,
                "C_rank": 0.95,
                "C_rank_boot_mean": 0.95,
                "C_rank_boot_var": 0.01,
                "C_rank_boot_sd": 0.1,
                "p_MC_C": 1e-3,
                "positive_consistent": False,
                "weighted_seed_rho": -0.6,
                "is_seed": False,
            },
            {
                "gene": "weak",
                "C_raw": 0.1,
                "C_rank": 0.6,
                "C_rank_boot_mean": 0.6,
                "C_rank_boot_var": 0.02,
                "C_rank_boot_sd": 0.14,
                "p_MC_C": 0.5,
                "positive_consistent": True,
                "weighted_seed_rho": 0.3,
                "is_seed": False,
            },
            {
                "gene": "onlyA",
                "C_raw": 0.4,
                "C_rank": 0.8,
                "C_rank_boot_mean": 0.8,
                "C_rank_boot_var": 0.01,
                "C_rank_boot_sd": 0.1,
                "p_MC_C": 1e-3,
                "positive_consistent": True,
                "weighted_seed_rho": 0.5,
                "is_seed": False,
            },
        ],
    )
    b = cohort_table(
        "B",
        [
            {
                "gene": "seed",
                "C_raw": 1,
                "C_rank": 1.0,
                "C_rank_boot_mean": 1.0,
                "C_rank_boot_var": 0.0,
                "C_rank_boot_sd": 0.0,
                "p_MC_C": 1e-4,
                "positive_consistent": True,
                "weighted_seed_rho": 0.9,
                "is_seed": True,
            },
            {
                "gene": "good",
                "C_raw": 0.5,
                "C_rank": 0.8,
                "C_rank_boot_mean": 0.8,
                "C_rank_boot_var": 0.04,
                "C_rank_boot_sd": 0.2,
                "p_MC_C": 2e-3,
                "positive_consistent": True,
                "weighted_seed_rho": 0.5,
                "is_seed": False,
            },
            {
                "gene": "neg",
                "C_raw": 0.5,
                "C_rank": 0.9,
                "C_rank_boot_mean": 0.9,
                "C_rank_boot_var": 0.01,
                "C_rank_boot_sd": 0.1,
                "p_MC_C": 1e-3,
                "positive_consistent": True,
                "weighted_seed_rho": 0.6,
                "is_seed": False,
            },
            {
                "gene": "weak",
                "C_raw": 0.1,
                "C_rank": 0.7,
                "C_rank_boot_mean": 0.7,
                "C_rank_boot_var": 0.02,
                "C_rank_boot_sd": 0.14,
                "p_MC_C": 1e-3,
                "positive_consistent": True,
                "weighted_seed_rho": 0.3,
                "is_seed": False,
            },
        ],
    )
    result = combine_cohorts(
        {"A": a, "B": b}, CombineConfig(enabled=True), symbols={"good": "GOOD"}
    )
    table = result.table.set_index("gene")
    assert list(result.signature["gene"]) == ["good"]
    assert table.loc["good", "pC_boot"] == pytest.approx(0.72)
    assert table.loc["good", "symbol"] == "GOOD"
    assert not table.loc["neg", "candidate"]  # negative edge in A
    assert not table.loc["weak", "candidate"]  # p too large in A
    assert not table.loc["seed", "candidate"]  # a seed
    assert not table.loc["onlyA", "candidate"] and np.isnan(table.loc["onlyA", "C_rank_B"])
    ivw = (0.9 / 0.01 + 0.8 / 0.04) / (1 / 0.01 + 1 / 0.04)
    assert table.loc["good", "C_IVW"] == pytest.approx(ivw)
    assert result.table["gene"].iloc[0] == "good"  # candidates first
    relaxed = combine_cohorts(
        {"A": a, "B": b},
        CombineConfig(
            enabled=True,
            require_positive_consistent=False,
            mc_p_cutoff=0.6,
            product_cutoff=0.3,
            exclude_seeds=False,
        ),
    )
    assert set(relaxed.signature["gene"]) == {"seed", "good", "neg", "weak"}
    legacy = legacy_two_cohort_names(result.table, result.cohorts)
    assert "Buffa_exact_candidate" in legacy.columns and "pC_hat_AxB_boot" in legacy.columns


def test_combination_without_bootstrap_or_null_is_tolerant() -> None:
    a = cohort_table(
        "A",
        [
            {
                "gene": "x",
                "C_raw": 0.5,
                "C_rank": 0.9,
                "positive_consistent": True,
                "weighted_seed_rho": 0.5,
                "is_seed": False,
                "n_edges": 2,
            }
        ],
    )
    result = combine_cohorts({"A": a}, CombineConfig(enabled=True))
    assert result.table.loc[0, "pC_boot"] == 0.9 and result.table.loc[0, "candidate"]
    assert any("no bootstrap" in n for n in result.notes) and any(
        "no null" in n for n in result.notes
    )


def test_median_abs_rank01_by_hand() -> None:
    values = np.array([[1.0, -4.0, 2.0], [3.0, 2.0, 2.0], [9.0, 9.0, 9.0], [0.0, 0.0, 0.0]])
    expr = ExpressionSet(values, pd.Index(["a", "b", "c", "d"]), pd.Index(["s1", "s2", "s3"]))
    table = median_abs_rank01(expr, ["a", "b", "c", "zz"], min_genes=3, dataset="T")
    np.testing.assert_allclose(table["E_median_abs_expression"], [3.0, 4.0, 2.0])
    np.testing.assert_allclose(table["HS_rank01"], [0.5, 1.0, 0.0])
    assert (
        table["n_signature_genes_present"].eq(3).all()
        and table["n_signature_genes_total"].eq(4).all()
    )
    with pytest.raises(DataError):
        median_abs_rank01(expr, ["a", "b"], min_genes=3)
    same = score_samples(expr, ["a", "b", "c"], ScoreConfig(enabled=True), "T")
    np.testing.assert_allclose(same["HS_rank01"], table["HS_rank01"])
