"""Combination of several training cohorts and selection of the signature.

Works on the per-cohort gene tables (observed connectivity, bootstrap summary, null p-values),
so it can be re-run from saved tables with different cutoffs.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from search.config import CombineConfig

PER_COHORT = [
    "C_raw",
    "C_rank",
    "C_rank_boot_mean",
    "C_rank_boot_var",
    "C_rank_boot_sd",
    "p_MC_C",
    "positive_consistent",
    "weighted_seed_rho",
]


@dataclass
class CombinedResult:
    table: pd.DataFrame
    signature: pd.DataFrame
    top: pd.DataFrame
    cohorts: list[str]
    notes: list[str] = field(default_factory=list)


def combine_cohorts(
    tables: dict[str, pd.DataFrame],
    cfg: CombineConfig,
    symbols: dict[str, str] | None = None,
) -> CombinedResult:
    """Outer join of the cohort tables on gene, product and inverse-variance scores, selection."""
    if not tables:
        raise ValueError("no cohort tables to combine")
    cohorts = list(tables)
    notes: list[str] = []
    merged: pd.DataFrame | None = None
    seed_flags: pd.Series | None = None
    for name, table in tables.items():
        part = table.set_index("gene")
        missing = [c for c in PER_COHORT if c not in part.columns]
        for column in missing:
            part[column] = np.nan
        if "C_rank_boot_mean" in missing:
            part["C_rank_boot_mean"] = part["C_rank"]
            notes.append(f"{name}: no bootstrap, observed ranks used as bootstrap means")
        if "p_MC_C" in missing:
            notes.append(f"{name}: no null model, the Monte Carlo criterion is skipped")
        if "positive_consistent" in missing:
            part["positive_consistent"] = True
        is_seed = (
            part["is_seed"].astype(bool)
            if "is_seed" in part.columns
            else pd.Series(False, index=part.index)
        )
        renamed = part[PER_COHORT].rename(columns={c: f"{c}_{name}" for c in PER_COHORT})
        merged = renamed if merged is None else merged.join(renamed, how="outer")
        seed_flags = (
            is_seed
            if seed_flags is None
            else seed_flags.reindex(merged.index, fill_value=False)
            | is_seed.reindex(merged.index, fill_value=False)
        )
    assert merged is not None and seed_flags is not None
    merged["is_seed"] = seed_flags.reindex(merged.index, fill_value=False).astype(bool)

    product_obs = np.ones(len(merged))
    product_boot = np.ones(len(merged))
    ivw_num = np.zeros(len(merged))
    ivw_den = np.zeros(len(merged))
    positive_all = np.ones(len(merged), dtype=bool)
    mc_all = np.ones(len(merged), dtype=bool)
    for name in cohorts:
        var = merged[f"C_rank_boot_var_{name}"].to_numpy(dtype=float)
        var_safe = np.where(np.isfinite(var), np.maximum(var, cfg.var_floor), cfg.var_floor)
        merged[f"var_{name}_safe"] = var_safe
        mean = merged[f"C_rank_boot_mean_{name}"].to_numpy(dtype=float)
        product_obs = product_obs * merged[f"C_rank_{name}"].to_numpy(dtype=float)
        product_boot = product_boot * mean
        ivw_num = ivw_num + mean / var_safe
        ivw_den = ivw_den + 1.0 / var_safe
        positive_all &= merged[f"positive_consistent_{name}"].fillna(False).astype(bool).to_numpy()
        p = merged[f"p_MC_C_{name}"].to_numpy(dtype=float)
        if np.isfinite(p).any():
            with np.errstate(invalid="ignore"):
                mc_all &= p < cfg.mc_p_cutoff
    merged["pC_obs"] = product_obs
    merged["pC_boot"] = product_boot
    with np.errstate(invalid="ignore", divide="ignore"):
        merged["C_IVW"] = ivw_num / ivw_den
    merged["positive_consistent_all"] = positive_all
    merged["MC_significant_all"] = mc_all
    with np.errstate(invalid="ignore"):
        merged["product_pass"] = merged["pC_boot"].to_numpy(dtype=float) > cfg.product_cutoff
    candidate = merged["MC_significant_all"] & merged["product_pass"]
    if cfg.require_positive_consistent:
        candidate &= merged["positive_consistent_all"]
    if cfg.exclude_seeds:
        candidate &= ~merged["is_seed"]
    merged["candidate"] = candidate.astype(bool)

    merged = merged.reset_index().rename(columns={"index": "gene"})
    if symbols:
        merged.insert(1, "symbol", merged["gene"].map(symbols))
    merged = merged.sort_values(
        ["candidate", "pC_boot", "pC_obs", "C_IVW"],
        ascending=[False, False, False, False],
        kind="stable",
        na_position="last",
    ).reset_index(drop=True)
    signature = merged[merged["candidate"]].reset_index(drop=True)
    top = signature.head(cfg.top_n).reset_index(drop=True)
    return CombinedResult(merged, signature, top, cohorts, notes)


def legacy_two_cohort_names(table: pd.DataFrame, cohorts: list[str]) -> pd.DataFrame:
    """Column aliases with the names of the GBM script for the two-cohort case (IT x FR)."""
    if len(cohorts) != 2:
        return table
    a, b = cohorts
    out = table.copy()
    aliases = {
        "pC_obs": f"pC_hat_{a}x{b}_obs",
        "pC_boot": f"pC_hat_{a}x{b}_boot",
        "C_IVW": f"C_hat_IVW_{a}{b}",
        "positive_consistent_all": f"positive_consistent_{a}{b}",
        "MC_significant_all": f"MC_significant_{a}{b}",
        "product_pass": "Buffa_product_pass",
        "candidate": "Buffa_exact_candidate",
    }
    for new, old in aliases.items():
        out[old] = out[new]
    return out
