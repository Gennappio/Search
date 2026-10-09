"""The engine: one pass from an expression matrix and a seed set to memberships,
connectivity and seed overlaps. Growth, bootstrap and the null model call this repeatedly."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from search.config import RunConfig
from search.core.association import associate
from search.core.connectivity import compute_connectivity, edge_statistics, nodes_degree
from search.core.edges import compute_edges
from search.core.membership import compute_membership
from search.core.overlap import jaccard_binary, shared_neighbourhood
from search.data import ExpressionSet
from search.errors import SeedError
from search.stats import rank01

FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]


@dataclass
class CoreResult:
    genes: pd.Index
    seeds: list[str]
    seed_idx: NDArray[np.intp]
    n_samples: int
    stat: FloatArray
    p: FloatArray | None
    delta: BoolArray
    gamma: FloatArray
    cluster_size: NDArray[np.int64]
    C_raw: FloatArray
    C_rank: FloatArray | None
    n_edges: NDArray[np.int64]
    n_positive_edges: NDArray[np.int64]
    n_negative_edges: NDArray[np.int64]
    positive_consistent: NDArray[np.bool_]
    weighted_seed_rho: FloatArray
    O: FloatArray  # noqa: E741 - the symbol of the paper
    S: FloatArray | None
    F: FloatArray | None
    nodes_degree: FloatArray
    config_hash: str
    coef_table: pd.DataFrame | None = None

    @property
    def n_genes(self) -> int:
        return len(self.genes)

    @property
    def n_seeds(self) -> int:
        return len(self.seeds)

    @property
    def is_seed(self) -> NDArray[np.bool_]:
        flags = np.zeros(self.n_genes, dtype=bool)
        flags[self.seed_idx] = True
        return flags

    def gene_table(self) -> pd.DataFrame:
        """One row per gene, same columns as the GBM script's gene table."""
        data: dict[str, Any] = {
            "gene": np.asarray(self.genes),
            "C_raw": self.C_raw,
            "C_rank": self.C_rank if self.C_rank is not None else np.full(self.n_genes, np.nan),
            "n_edges": self.n_edges,
            "n_positive_edges": self.n_positive_edges,
            "n_negative_edges": self.n_negative_edges,
            "positive_consistent": self.positive_consistent,
            "weighted_seed_rho": self.weighted_seed_rho,
            "is_seed": self.is_seed,
        }
        return pd.DataFrame(data)

    def neighbourhoods(self) -> dict[str, pd.Series]:
        """Per seed, the association values of the genes in its cluster (SEARCH's CL)."""
        out: dict[str, pd.Series] = {}
        for j, seed in enumerate(self.seeds):
            mask = self.delta[:, j]
            out[seed] = pd.Series(self.stat[mask, j], index=self.genes[mask], name=seed)
        return out

    def seed_pairs(self) -> pd.DataFrame:
        """Long table of seed pairs with O, S and F."""
        rows = []
        for a, sa in enumerate(self.seeds):
            for b, sb in enumerate(self.seeds):
                rows.append(
                    {
                        "seed_x": sa,
                        "seed_y": sb,
                        "O": self.O[a, b],
                        "S": self.S[a, b] if self.S is not None else np.nan,
                        "F": self.F[a, b] if self.F is not None else np.nan,
                    }
                )
        return pd.DataFrame(rows)


def resolve_seeds(
    expr: ExpressionSet, seeds: list[str], universe: list[str] | None
) -> tuple[ExpressionSet, list[str]]:
    """Restrict the matrix to the universe and check that the seeds are inside it."""
    unique_seeds = list(dict.fromkeys(str(s) for s in seeds))
    if len(unique_seeds) < 2:
        raise SeedError("at least two seeds are needed: with one seed connectivity is undefined")
    work = expr
    if universe is not None:
        wanted = {str(g) for g in universe}
        keep = expr.genes[expr.genes.isin(wanted)]
        work = expr.subset_genes(keep)
    missing = [s for s in unique_seeds if s not in work.genes]
    if missing:
        where = "the universe" if universe is not None else "the matrix"
        raise SeedError(f"{len(missing)} seeds are not in {where}: {missing[:10]}")
    return work, unique_seeds


def run_core(
    expr: ExpressionSet,
    seeds: list[str],
    cfg: RunConfig,
    *,
    universe: list[str] | None = None,
    sample_idx: NDArray[np.intp] | None = None,
) -> CoreResult:
    """One pass of the engine.

    ``sample_idx`` selects columns by position (bootstrap or jackknife); passing it explicitly
    makes any two implementations comparable on identical resamples.
    """
    work, seed_list = resolve_seeds(expr, seeds, universe)
    if sample_idx is not None:
        work = work.take_samples(sample_idx)
    seed_idx = work.gene_positions(seed_list)

    assoc = associate(work.values, seed_idx, cfg.association, work.sample_table, seed_list)
    delta = compute_edges(assoc, cfg.edge, seed_idx)
    gamma, cluster_size = compute_membership(delta, assoc.stat, cfg.membership)
    c_raw = compute_connectivity(gamma, seed_idx, cfg.connectivity)
    c_rank = rank01(c_raw) if cfg.connectivity.rank == "fractional" else None
    stats = edge_statistics(delta, assoc.stat, gamma)
    o = jaccard_binary(delta)
    want_sf = cfg.overlap.shared_neighbourhood == "always" or (
        cfg.overlap.shared_neighbourhood == "auto" and len(seed_list) <= cfg.overlap.auto_max_seeds
    )
    s, f = shared_neighbourhood(gamma, assoc.stat, seed_idx) if want_sf else (None, None)

    return CoreResult(
        genes=work.genes,
        seeds=seed_list,
        seed_idx=seed_idx,
        n_samples=work.n_samples,
        stat=assoc.stat,
        p=assoc.p,
        delta=delta,
        gamma=gamma,
        cluster_size=cluster_size,
        C_raw=c_raw,
        C_rank=c_rank,
        n_edges=np.asarray(stats["n_edges"], dtype=np.int64),
        n_positive_edges=np.asarray(stats["n_positive_edges"], dtype=np.int64),
        n_negative_edges=np.asarray(stats["n_negative_edges"], dtype=np.int64),
        positive_consistent=np.asarray(stats["positive_consistent"], dtype=bool),
        weighted_seed_rho=np.asarray(stats["weighted_seed_rho"], dtype=np.float64),
        O=o,
        S=s,
        F=f,
        nodes_degree=nodes_degree(delta, seed_idx),
        config_hash=cfg.core_hash(),
        coef_table=assoc.coef_table,
    )
