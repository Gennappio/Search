"""Comparison of two runs: how much their rankings, p-values and signatures agree."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from search.export.tables import read_cohort_tables
from search.run import resolve_run_dir


@dataclass
class RunComparison:
    run_a: Path
    run_b: Path
    cohort_metrics: pd.DataFrame
    combined_metrics: pd.DataFrame
    signature_jaccard: float
    only_a: list[str]
    only_b: list[str]

    def summary(self) -> str:
        lines = [f"A: {self.run_a}", f"B: {self.run_b}", ""]
        if len(self.cohort_metrics):
            lines.append(self.cohort_metrics.to_string(index=False))
            lines.append("")
        if len(self.combined_metrics):
            lines.append(self.combined_metrics.to_string(index=False))
            lines.append("")
        lines.append(
            f"signature Jaccard {self.signature_jaccard:.3f}; only in A: {len(self.only_a)}; "
            f"only in B: {len(self.only_b)}"
        )
        return "\n".join(lines)


def _correlations(a: pd.Series, b: pd.Series) -> tuple[float, float, int]:
    joined = pd.concat([a.astype(float), b.astype(float)], axis=1, join="inner").dropna()
    if len(joined) < 3:
        return float("nan"), float("nan"), int(len(joined))
    pearson = float(joined.iloc[:, 0].corr(joined.iloc[:, 1]))
    spearman = float(joined.iloc[:, 0].corr(joined.iloc[:, 1], method="spearman"))
    return pearson, spearman, int(len(joined))


def _top_overlap(a: pd.Series, b: pd.Series, n: int) -> float:
    top_a = set(a.dropna().sort_values(ascending=False).head(n).index)
    top_b = set(b.dropna().sort_values(ascending=False).head(n).index)
    return len(top_a & top_b) / n if n else float("nan")


def compare_runs(run_a: Path | str, run_b: Path | str, top_n: int = 50) -> RunComparison:
    a_dir, b_dir = resolve_run_dir(run_a), resolve_run_dir(run_b)
    tables_a, tables_b = (
        read_cohort_tables(a_dir / "results"),
        read_cohort_tables(b_dir / "results"),
    )
    rows = []
    for name in sorted(set(tables_a) & set(tables_b)):
        ta, tb = tables_a[name].set_index("gene"), tables_b[name].set_index("gene")
        for column in (
            "C_raw",
            "C_rank",
            "C_rank_boot_mean",
            "C_raw_boot_mean",
            "p_MC_C",
            "D_final",
        ):
            if column in ta.columns and column in tb.columns:
                pearson, spearman, n = _correlations(ta[column], tb[column])
                rows.append(
                    {
                        "cohort": name,
                        "column": column,
                        "n_genes": n,
                        "pearson": pearson,
                        "spearman": spearman,
                        f"top{top_n}_overlap": _top_overlap(ta[column], tb[column], top_n),
                    }
                )
    cohort_metrics = pd.DataFrame(rows)

    combined_rows = []
    sig_a: set[str] = set()
    sig_b: set[str] = set()
    full_a, full_b = (
        a_dir / "results" / "combined_full_ranking.csv",
        b_dir / "results" / "combined_full_ranking.csv",
    )
    if full_a.exists() and full_b.exists():
        fa, fb = pd.read_csv(full_a).set_index("gene"), pd.read_csv(full_b).set_index("gene")
        for column in ("pC_boot", "pC_obs", "C_IVW"):
            if column in fa.columns and column in fb.columns:
                pearson, spearman, n = _correlations(fa[column], fb[column])
                combined_rows.append(
                    {
                        "column": column,
                        "n_genes": n,
                        "pearson": pearson,
                        "spearman": spearman,
                        f"top{top_n}_overlap": _top_overlap(fa[column], fb[column], top_n),
                    }
                )
        if "candidate" in fa.columns and "candidate" in fb.columns:
            sig_a = set(fa.index[fa["candidate"].astype(bool)].astype(str))
            sig_b = set(fb.index[fb["candidate"].astype(bool)].astype(str))
    union = sig_a | sig_b
    jaccard = len(sig_a & sig_b) / len(union) if union else float("nan")
    return RunComparison(
        a_dir,
        b_dir,
        cohort_metrics,
        pd.DataFrame(combined_rows),
        jaccard,
        sorted(sig_a - sig_b),
        sorted(sig_b - sig_a),
    )


def variants_grid(grid: dict[str, list[object]]) -> list[dict[str, object]]:
    """Every combination of the grid, as override dictionaries."""
    combos: list[dict[str, object]] = [{}]
    for key, values in grid.items():
        combos = [{**combo, key: value} for combo in combos for value in values]
    return combos


def variant_label(overrides: dict[str, object]) -> str:
    return (
        "__".join(f"{key.split('.')[-1]}-{value}" for key, value in overrides.items()) or "baseline"
    )


def summarise_variants(runs: dict[str, Path], baseline: str, top_n: int = 50) -> pd.DataFrame:
    """One row per variant: candidates, agreement with the baseline, top genes."""
    rows = []
    base_dir = resolve_run_dir(runs[baseline])
    base_full = pd.read_csv(base_dir / "results" / "combined_full_ranking.csv").set_index("gene")
    for label, path in runs.items():
        run_dir = resolve_run_dir(path)
        full = pd.read_csv(run_dir / "results" / "combined_full_ranking.csv").set_index("gene")
        comparison = compare_runs(base_dir, run_dir, top_n)
        score = "pC_boot" if "pC_boot" in full.columns else "pC_obs"
        pearson, spearman, _ = _correlations(base_full[score], full[score])
        top = full[score].dropna().sort_values(ascending=False).head(10)
        names = (
            full["symbol"].reindex(top.index)
            if "symbol" in full.columns
            else pd.Series(top.index, index=top.index)
        )
        rows.append(
            {
                "variant": label,
                "n_candidates": int(full["candidate"].astype(bool).sum())
                if "candidate" in full.columns
                else np.nan,
                "spearman_vs_baseline": spearman,
                f"top{top_n}_overlap_vs_baseline": _top_overlap(
                    base_full[score], full[score], top_n
                ),
                "signature_jaccard_vs_baseline": comparison.signature_jaccard,
                "top10": ", ".join(str(n) if pd.notna(n) else str(g) for g, n in names.items()),
                "run_dir": str(run_dir),
            }
        )
    return pd.DataFrame(rows)
