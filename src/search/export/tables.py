"""Result tables on disk. Column names follow the GBM script so downstream analyses keep working."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from search.bootstrap import BootstrapSummary
from search.combine import CombinedResult, legacy_two_cohort_names
from search.core.corepass import CoreResult
from search.growth import GrowthResult
from search.null import NullResult


def cohort_gene_table(
    name: str, core: CoreResult, bootstrap: BootstrapSummary | None, null: NullResult | None
) -> pd.DataFrame:
    table = core.gene_table()
    if bootstrap is not None:
        table = table.merge(bootstrap.table(), on="gene", how="left")
    if null is not None:
        table = table.merge(null.table(), on="gene", how="left")
    table.insert(0, "dataset", name)
    return table


def seed_cluster_table(core: CoreResult) -> pd.DataFrame:
    rows = []
    for j, seed in enumerate(core.seeds):
        mask = core.delta[:, j]
        rows.append(
            pd.DataFrame(
                {
                    "seed": seed,
                    "gene": np.asarray(core.genes[mask]),
                    "stat": core.stat[mask, j],
                    "p": core.p[mask, j] if core.p is not None else np.nan,
                    "gamma": core.gamma[mask, j],
                }
            )
        )
    return (
        pd.concat(rows, ignore_index=True)
        if rows
        else pd.DataFrame(columns=["seed", "gene", "stat", "p", "gamma"])
    )


def write_cohort_tables(
    results_dir: Path,
    name: str,
    core: CoreResult,
    bootstrap: BootstrapSummary | None,
    null: NullResult | None,
    gene_table: pd.DataFrame | None = None,
) -> list[Path]:
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    def save(frame: pd.DataFrame, filename: str) -> None:
        path = results_dir / filename
        frame.to_csv(path, index=False)
        written.append(path)

    table = gene_table if gene_table is not None else cohort_gene_table(name, core, bootstrap, null)
    save(table, f"{name}_gene_connectivity.csv")
    pairs = core.seed_pairs()
    pairs.insert(0, "dataset", name)
    save(pairs[["dataset", "seed_x", "seed_y", "S"]], f"{name}_shared_neighbourhood_S.csv")
    save(pairs[["dataset", "seed_x", "seed_y", "F"]], f"{name}_sign_concordance_F.csv")
    save(pairs[["dataset", "seed_x", "seed_y", "O"]], f"{name}_seed_overlap_O.csv")
    clusters = seed_cluster_table(core)
    clusters.insert(0, "dataset", name)
    save(clusters, f"{name}_seed_clusters.csv")
    if null is not None and len(null.random_seed_summary):
        save(null.random_seed_summary, f"{name}_MC_random_seed_S_summary.csv")
    seeds = pd.DataFrame(
        {
            "dataset": name,
            "seed": core.seeds,
            "cluster_size": core.cluster_size,
            "nodes_degree": core.nodes_degree,
        }
    )
    save(seeds, f"{name}_seeds.csv")
    if core.coef_table is not None:
        coef = core.coef_table.copy()
        coef.insert(0, "dataset", name)
        save(coef, f"{name}_covariate_coefficients.csv")
    return written


def write_combined_tables(
    results_dir: Path, combined: CombinedResult, prefix: str = "combined"
) -> list[Path]:
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    full = legacy_two_cohort_names(combined.table, combined.cohorts)
    paths = {
        f"{prefix}_full_ranking.csv": full,
        f"{prefix}_signature.csv": legacy_two_cohort_names(combined.signature, combined.cohorts),
        f"{prefix}_signature_top.csv": legacy_two_cohort_names(combined.top, combined.cohorts),
    }
    written = []
    for filename, frame in paths.items():
        path = results_dir / filename
        frame.to_csv(path, index=False)
        written.append(path)
    return written


def write_scores(results_dir: Path, scores: pd.DataFrame, filename: str) -> Path:
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    path = results_dir / filename
    scores.to_csv(path, index=False)
    return path


def read_cohort_tables(results_dir: Path) -> dict[str, pd.DataFrame]:
    """The per-cohort gene tables of a run directory, keyed by cohort name."""
    tables: dict[str, pd.DataFrame] = {}
    for path in sorted(Path(results_dir).glob("*_gene_connectivity.csv")):
        frame = pd.read_csv(path)
        name = (
            str(frame["dataset"].iloc[0])
            if "dataset" in frame.columns and len(frame)
            else path.name.split("_gene_connectivity")[0]
        )
        tables[name] = frame
    return tables


def write_growth_tables(results_dir: Path, name: str, growth: GrowthResult) -> list[Path]:
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    frames = {
        f"{name}_growth_iterations.csv": growth.iteration_table(),
        f"{name}_growth_events.csv": growth.event_table(),
        f"{name}_growth_connectivity.csv": growth.connectivity_table(),
        f"{name}_growth_final_seeds.csv": pd.DataFrame(
            {
                "gene": growth.final_seeds,
                "initial": [g in set(growth.initial_seeds) for g in growth.final_seeds],
            }
        ),
    }
    written = []
    for filename, frame in frames.items():
        frame.insert(0, "dataset", name)
        path = results_dir / filename
        frame.to_csv(path, index=False)
        written.append(path)
    return written


def write_external_null(results_dir: Path, name: str, table: pd.DataFrame) -> Path:
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    path = results_dir / f"{name}_external_null.csv"
    table.to_csv(path, index=False)
    return path
