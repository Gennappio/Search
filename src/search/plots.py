"""Figures from a run directory. Optional: needs the ``plots`` extra (matplotlib, networkx).
They illustrate results; nothing here is validated."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from search.export.tables import read_cohort_tables


def _require() -> tuple[Any, Any]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import networkx as nx
    except ImportError as exc:  # pragma: no cover
        raise ImportError("plots need matplotlib and networkx: pip install search[plots]") from exc
    return plt, nx


def seed_network_plot(
    pairs: pd.DataFrame, column: str, path: Path, title: str, threshold: float = 0.0
) -> Path:
    """Seeds as nodes, pairwise overlap (S or O) or correlation as weighted edges."""
    plt, nx = _require()
    graph = nx.Graph()
    seeds = sorted(set(pairs["seed_x"]) | set(pairs["seed_y"]))
    graph.add_nodes_from(seeds)
    for _, row in pairs.iterrows():
        a, b, w = row["seed_x"], row["seed_y"], row[column]
        if a < b and pd.notna(w) and abs(float(w)) > threshold:
            graph.add_edge(a, b, weight=float(w))
    pos = nx.circular_layout(graph) if len(seeds) <= 30 else nx.spring_layout(graph, seed=0)
    weights = (
        np.array([abs(graph.edges[e]["weight"]) for e in graph.edges])
        if graph.edges
        else np.array([])
    )
    colors = ["#d62728" if graph.edges[e]["weight"] > 0 else "#1f77b4" for e in graph.edges]
    fig, ax = plt.subplots(figsize=(7, 7))
    nx.draw_networkx_nodes(graph, pos, node_color="#9ecae1", node_size=600, ax=ax)
    if len(weights):
        nx.draw_networkx_edges(
            graph, pos, width=1 + 5 * weights / weights.max(), edge_color=colors, ax=ax
        )
    nx.draw_networkx_labels(graph, pos, font_size=7, ax=ax)
    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def rank_scatter_plot(
    full: pd.DataFrame, cohorts: list[str], path: Path, product_cutoff: float
) -> Path:
    plt, _ = _require()
    a, b = cohorts[:2]
    frame = full[~full["is_seed"].astype(bool)]
    fig, ax = plt.subplots(figsize=(6, 5))
    colors = np.where(frame["candidate"].astype(bool), "#d62728", "#7f7f7f")
    ax.scatter(frame[f"C_rank_{a}"], frame[f"C_rank_{b}"], s=6, c=colors, alpha=0.6)
    cut = np.sqrt(product_cutoff)
    ax.axhline(cut, ls="--", c="k", lw=0.8)
    ax.axvline(cut, ls="--", c="k", lw=0.8)
    ax.set_xlabel(f"{a} fractional rank of C")
    ax.set_ylabel(f"{b} fractional rank of C")
    ax.set_title("Connectivity ranks, candidates in red")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def product_histogram(full: pd.DataFrame, path: Path, product_cutoff: float) -> Path:
    plt, _ = _require()
    values = full.loc[~full["is_seed"].astype(bool), "pC_boot"].dropna()
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(values, bins=60, color="#9ecae1", edgecolor="white")
    ax.axvline(product_cutoff, ls="--", c="k")
    ax.set_xlabel("product of bootstrap-mean ranks")
    ax.set_ylabel("genes")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def score_violin(scores: pd.DataFrame, path: Path) -> Path:
    plt, _ = _require()
    groups = [g["HS_rank01"].dropna().to_numpy() for _, g in scores.groupby("dataset", sort=True)]
    labels = [str(name) for name, _ in scores.groupby("dataset", sort=True)]
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.violinplot(groups, showmedians=True)
    ax.set_xticks(range(1, len(labels) + 1), labels)
    ax.set_ylabel("signature score, rank 0-1")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_run(run_dir: Path | str, product_cutoff: float = 0.5) -> list[Path]:
    """Every figure the results of a run allow; returns the written paths."""
    run_dir = Path(run_dir)
    results = run_dir / "results"
    out = run_dir / "plots"
    out.mkdir(exist_ok=True)
    written: list[Path] = []
    for name in read_cohort_tables(results):
        for kind, column in (("shared_neighbourhood_S", "S"), ("seed_overlap_O", "O")):
            path = results / f"{name}_{kind}.csv"
            if path.exists():
                pairs = pd.read_csv(path)
                if pairs[column].notna().any():
                    written.append(
                        seed_network_plot(
                            pairs,
                            column,
                            out / f"{name}_seed_network_{column}.png",
                            f"{name}: seed {column}",
                        )
                    )
    full_path = results / "combined_full_ranking.csv"
    if full_path.exists():
        full = pd.read_csv(full_path)
        cohorts = [
            c.removeprefix("C_rank_")
            for c in full.columns
            if c.startswith("C_rank_") and not c.startswith("C_rank_boot")
        ]
        if len(cohorts) >= 2 and "candidate" in full.columns:
            written.append(
                rank_scatter_plot(full, cohorts, out / "rank_scatter.png", product_cutoff)
            )
        if "pC_boot" in full.columns:
            written.append(product_histogram(full, out / "product_histogram.png", product_cutoff))
    scores_path = results / "scores_signature_all.csv"
    if scores_path.exists():
        written.append(score_violin(pd.read_csv(scores_path), out / "score_violin.png"))
    return written
