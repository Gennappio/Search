"""Recursive growth of the seed set (SEARCH 2014).

Each iteration: the neighbourhoods of the new seeds are computed with the engine; new seeds whose
neighbourhood overlaps too little with the others are pruned (Omax < O_lower); the connectivity D
of every universe gene to the current seeds is computed; seeds with D < D_lower are pruned; genes
with D > D_higher become seeds. Pruned genes leave the universe for good. The loop stops when
nothing changes, when the universe is empty, or at max_iter.

Differences from the R original are deliberate and listed in docs/decisions.md: pruned seeds are
removed from every neighbourhood, Omax is the maximum over every other seed, matrices are kept
for every iteration, and the standard-deviation prefilter is keyed by gene name.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from search.config import RunConfig
from search.core.connectivity import compute_connectivity
from search.core.corepass import resolve_seeds, run_core
from search.core.membership import compute_membership
from search.core.overlap import jaccard_binary, omax
from search.data import ExpressionSet
from search.stats import rank01, row_sds

FloatArray = NDArray[np.float64]
log = logging.getLogger("search")


@dataclass
class GrowthIteration:
    iteration: int
    seeds_in: list[str]
    universe_size: int
    D: pd.Series  # over the universe after the D pruning of this iteration
    Omax: pd.Series  # over the seeds checked for overlap (empty at iteration 1)
    removed_O: list[str]
    removed_D: list[str]
    added: list[str]


@dataclass
class GrowthResult:
    initial_seeds: list[str]
    final_seeds: list[str]
    universe: pd.Index  # the initial universe, in matrix order
    final_universe: pd.Index
    iterations: list[GrowthIteration]
    removed: list[dict[str, object]] = field(default_factory=list)
    stop_reason: str = ""
    neighbourhoods: dict[str, pd.Series] = field(default_factory=dict)

    @property
    def n_iterations(self) -> int:
        return len(self.iterations)

    def final_D(self) -> FloatArray:
        """Connectivity of the last iteration, aligned to the initial universe (NaN if pruned)."""
        if not self.iterations:
            return np.full(len(self.universe), np.nan)
        return self.iterations[-1].D.reindex(self.universe).to_numpy(dtype=np.float64)

    def iteration_table(self) -> pd.DataFrame:
        rows = [
            {
                "iteration": it.iteration,
                "n_seeds_in": len(it.seeds_in),
                "n_universe": it.universe_size,
                "n_removed_O": len(it.removed_O),
                "n_removed_D": len(it.removed_D),
                "n_added": len(it.added),
                "n_seeds_out": len(it.seeds_in)
                - len(it.removed_O)
                - len(it.removed_D)
                + len(it.added),
            }
            for it in self.iterations
        ]
        table = pd.DataFrame(rows)
        table["stop_reason"] = ""
        if len(table):
            table.loc[table.index[-1], "stop_reason"] = self.stop_reason
        return table

    def event_table(self) -> pd.DataFrame:
        rows: list[dict[str, object]] = []
        for it in self.iterations:
            rows.extend(
                {
                    "iteration": it.iteration,
                    "gene": g,
                    "event": "removed_O",
                    "value": it.Omax.get(g, np.nan),
                }
                for g in it.removed_O
            )
            rows.extend(
                {"iteration": it.iteration, "gene": g, "event": "removed_D", "value": np.nan}
                for g in it.removed_D
            )
            rows.extend(
                {
                    "iteration": it.iteration,
                    "gene": g,
                    "event": "added",
                    "value": it.D.get(g, np.nan),
                }
                for g in it.added
            )
        return pd.DataFrame(rows, columns=["iteration", "gene", "event", "value"])

    def connectivity_table(self) -> pd.DataFrame:
        parts = [
            pd.DataFrame({"iteration": it.iteration, "gene": it.D.index, "D": it.D.to_numpy()})
            for it in self.iterations
        ]
        return (
            pd.concat(parts, ignore_index=True)
            if parts
            else pd.DataFrame(columns=["iteration", "gene", "D"])
        )


def _membership_from_neighbourhoods(
    universe: pd.Index, seeds: list[str], neighbourhoods: dict[str, pd.Series]
) -> tuple[NDArray[np.bool_], FloatArray]:
    """delta and stat matrices (universe x seeds) rebuilt from the stored neighbourhoods."""
    delta = np.zeros((len(universe), len(seeds)), dtype=bool)
    stat = np.full((len(universe), len(seeds)), np.nan)
    for j, seed in enumerate(seeds):
        hood = neighbourhoods[seed]
        hood = hood[hood.index.isin(universe)]
        rows = universe.get_indexer(hood.index)
        delta[rows, j] = True
        stat[rows, j] = hood.to_numpy(dtype=np.float64)
    return delta, stat


def run_growth(
    expr: ExpressionSet,
    seeds: list[str],
    cfg: RunConfig,
    *,
    universe: list[str] | None = None,
    sample_idx: NDArray[np.intp] | None = None,
) -> GrowthResult:
    gc = cfg.growth
    work, seed_list = resolve_seeds(expr, seeds, universe)
    if sample_idx is not None:
        work = work.take_samples(sample_idx)
    initial_universe = work.genes
    removed: list[dict[str, object]] = []

    # prefilter on standard deviation, keyed by gene name
    if gc.sd_threshold > 0:
        sds = pd.Series(row_sds(work.values), index=work.genes)
        low = [str(g) for g in sds.index[sds < gc.sd_threshold]]
        if low:
            removed.extend({"gene": g, "reason": "low_sd", "iteration": 0} for g in low)
            work = work.subset_genes([g for g in work.genes if g not in set(low)])
            seed_list = [s for s in seed_list if s not in set(low)]
            log.info("growth: %d genes removed by the sd threshold", len(low))

    current_seeds = list(seed_list)
    neighbourhoods: dict[str, pd.Series] = {}
    added: list[str] = []
    iterations: list[GrowthIteration] = []
    stop_reason = "max_iter"

    for iteration in range(1, gc.max_iter + 1):
        if len(current_seeds) < 2:
            stop_reason = "fewer_than_two_seeds"
            break
        seeds_in = list(current_seeds)
        to_compute = seeds_in if iteration == 1 else added
        removed_o: list[str] = []
        omax_series = pd.Series(dtype=float)
        if to_compute:
            core = run_core(work, to_compute, cfg) if len(to_compute) >= 2 else None
            if core is None:
                # a single new seed: compute it together with an existing seed, keep its column
                partner = next(s for s in seeds_in if s != to_compute[0])
                core = run_core(work, [to_compute[0], partner], cfg)
            hoods = core.neighbourhoods()
            for seed in to_compute:
                neighbourhoods[seed] = hoods[seed]
        if iteration > 1:
            delta, _ = _membership_from_neighbourhoods(work.genes, current_seeds, neighbourhoods)
            o = jaccard_binary(delta)
            all_omax = pd.Series(omax(o), index=current_seeds)
            candidates = current_seeds if not gc.keep_seeds else added
            omax_series = all_omax.reindex(candidates)
            removed_o = [str(g) for g in omax_series.index[omax_series < gc.O_lower]]
            if removed_o:
                current_seeds = [s for s in current_seeds if s not in set(removed_o)]
                work = work.subset_genes([g for g in work.genes if g not in set(removed_o)])
                for seed in removed_o:
                    neighbourhoods.pop(seed, None)
                removed.extend(
                    {"gene": g, "reason": "omax", "iteration": iteration} for g in removed_o
                )
        if len(current_seeds) < 2:
            iterations.append(
                GrowthIteration(
                    iteration,
                    seeds_in,
                    work.n_genes,
                    pd.Series(dtype=float),
                    omax_series,
                    removed_o,
                    [],
                    [],
                )
            )
            stop_reason = "fewer_than_two_seeds"
            break

        delta, stat = _membership_from_neighbourhoods(work.genes, current_seeds, neighbourhoods)
        seed_idx = work.gene_positions(current_seeds)
        gamma, _ = compute_membership(delta, stat, cfg.membership)
        d_values = compute_connectivity(gamma, seed_idx, cfg.connectivity)
        d = pd.Series(d_values, index=work.genes)

        seed_d = d.reindex(current_seeds)
        removed_d = [str(g) for g in seed_d.index[seed_d < gc.D_lower]]
        if removed_d:
            current_seeds = [s for s in current_seeds if s not in set(removed_d)]
            work = work.subset_genes([g for g in work.genes if g not in set(removed_d)])
            for seed in removed_d:
                neighbourhoods.pop(seed, None)
            d = d.drop(removed_d)
            removed.extend(
                {"gene": g, "reason": "D_lower", "iteration": iteration} for g in removed_d
            )

        not_seed = d[~d.index.isin(current_seeds)]
        added = [str(g) for g in not_seed.index[not_seed > gc.D_higher]]
        current_seeds = current_seeds + added
        iterations.append(
            GrowthIteration(
                iteration, seeds_in, work.n_genes, d, omax_series, removed_o, removed_d, added
            )
        )
        log.info(
            "growth iteration %d: %d seeds in, %d removed by O, %d by D, %d added, universe %d",
            iteration,
            len(seeds_in),
            len(removed_o),
            len(removed_d),
            len(added),
            work.n_genes,
        )
        if len(added) + len(removed_d) + len(removed_o) == 0:
            stop_reason = "converged"
            break
        if work.n_genes == 0:
            stop_reason = "universe_empty"
            break
    else:
        stop_reason = "max_iter"

    return GrowthResult(
        initial_seeds=list(seed_list),
        final_seeds=list(current_seeds),
        universe=initial_universe,
        final_universe=work.genes,
        iterations=iterations,
        removed=removed,
        stop_reason=stop_reason,
        neighbourhoods=neighbourhoods,
    )


def growth_statistic(result: GrowthResult, statistic: str) -> FloatArray:
    """Final connectivity aligned to the initial universe; fractional rank if asked."""
    d = result.final_D()
    return rank01(d) if statistic == "C_rank" else d
