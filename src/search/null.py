"""Null model: the same pass with random seed sets, repeated, to estimate how much connectivity
a gene reaches by chance. Pooled or per-gene empirical p-values."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from search.bootstrap import stage_hash, statistic_of
from search.checkpoint import ChunkMeta, ChunkStore, chunk_ranges
from search.config import RunConfig
from search.core.corepass import resolve_seeds, run_core
from search.data import ExpressionSet
from search.parallel import parallel_map
from search.resampling import derive_rng, eligible_null_genes, random_seed_set

FloatArray = NDArray[np.float64]
log = logging.getLogger("search")
TOLERANCE = 1e-12


@dataclass
class NullResult:
    statistic: str
    pooling: str
    n_rep: int
    genes: pd.Index
    p: FloatArray
    count_ge: NDArray[np.int64]
    n_null: NDArray[np.int64]
    random_seed_summary: pd.DataFrame

    def table(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "gene": np.asarray(self.genes),
                "p_MC_C": self.p,
                "MC_count_ge": self.count_ge,
                "MC_n_null": self.n_null,
            }
        )


@dataclass
class _Job:
    work: ExpressionSet
    eligible: NDArray[np.str_]
    k: int
    cfg: RunConfig
    rep: int

    def run(self) -> tuple[int, FloatArray, float, float, list[str]]:
        rng = derive_rng(self.cfg.null.seed_base, self.rep)
        random_seeds = random_seed_set(self.eligible, self.k, rng)
        res = run_core(self.work, random_seeds, self.cfg)
        values = statistic_of(res, self.cfg.null.statistic)
        mean_s = median_s = float("nan")
        if res.S is not None and self.k > 1:
            upper = res.S[np.triu_indices(self.k, 1)]
            upper = upper[np.isfinite(upper)]
            if upper.size:
                mean_s, median_s = float(upper.mean()), float(np.median(upper))
        return self.rep, values, mean_s, median_s, random_seeds


def _run_job(job: _Job) -> tuple[int, FloatArray, float, float, list[str]]:
    return job.run()


def empirical_p(count_ge: NDArray[np.int64], n_null: NDArray[np.int64]) -> FloatArray:
    """(1 + count_ge) / (1 + n_null), the pseudo-count form that never returns 0."""
    with np.errstate(invalid="ignore", divide="ignore"):
        p = (1.0 + count_ge) / (1.0 + n_null)
    return np.where(n_null > 0, p, np.nan)


def pooled_counts(observed: FloatArray, chunk: FloatArray) -> tuple[NDArray[np.int64], int]:
    """How many values of the chunk, over all genes, are >= every observed value."""
    values = np.sort(chunk[np.isfinite(chunk)])
    n = len(values)
    ok = np.isfinite(observed)
    count = np.zeros(len(observed), dtype=np.int64)
    count[ok] = n - np.searchsorted(values, observed[ok] - TOLERANCE, side="right")
    return count, n


def per_gene_counts(
    observed: FloatArray, chunk: FloatArray
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    finite = np.isfinite(chunk)
    with np.errstate(invalid="ignore"):
        ge = finite & (chunk >= observed[:, None] - TOLERANCE)
    return ge.sum(axis=1).astype(np.int64), finite.sum(axis=1).astype(np.int64)


def run_internal_null(
    expr: ExpressionSet,
    seeds: list[str],
    cfg: RunConfig,
    *,
    universe: list[str] | None,
    observed: FloatArray,
    ckpt_dir: Path | str,
    input_hash: str,
    resume: bool = True,
    n_jobs: int = 1,
) -> NullResult:
    nc = cfg.null
    work, seed_list = resolve_seeds(expr, seeds, universe)
    if len(observed) != work.n_genes:
        raise ValueError("observed statistic must be aligned with the universe genes")
    eligible = eligible_null_genes(work, seed_list, nc.eligible_min_sd)
    k = len(seed_list)
    store = ChunkStore(ckpt_dir, "null", cfg.checkpoint_dtype)
    config_hash = stage_hash(cfg, "null", nc.statistic, str(nc.eligible_min_sd))
    ranges = chunk_ranges(nc.n_rep, nc.chunk_size)
    for chunk_no, (start, end) in enumerate(ranges, start=1):
        meta = ChunkMeta(
            "null",
            start,
            end,
            nc.seed_base,
            config_hash,
            input_hash,
            {"k": k, "n_eligible": int(len(eligible))},
        )
        if resume and store.is_valid(chunk_no, meta):
            log.info("null chunk %d/%d already done, skipping", chunk_no, len(ranges))
            continue
        log.info("null chunk %d/%d, replicates %d-%d", chunk_no, len(ranges), start, end)
        jobs = [_Job(work, eligible, k, cfg, rep) for rep in range(start, end + 1)]
        outputs = parallel_map(_run_job, jobs, n_jobs)
        columns = {rep: values for rep, values, _, _, _ in outputs}
        summary = pd.DataFrame(
            {
                "rep": [rep for rep, *_ in outputs],
                "mean_S": [m for _, _, m, _, _ in outputs],
                "median_S": [m for _, _, _, m, _ in outputs],
                "seeds": [",".join(s) for *_, s in outputs],
            }
        )
        store.write(chunk_no, work.genes, columns, meta, {"seedsummary": summary})

    count_ge = np.zeros(work.n_genes, dtype=np.int64)
    n_null = np.zeros(work.n_genes, dtype=np.int64)
    summaries = []
    for chunk_no in range(1, len(ranges) + 1):
        chunk = store.read(chunk_no).reindex(work.genes).to_numpy(dtype=np.float64)
        if nc.pooling == "pooled":
            counts, n = pooled_counts(observed, chunk)
            count_ge += counts
            n_null += n
        else:
            counts, n_g = per_gene_counts(observed, chunk)
            count_ge += counts
            n_null += n_g
        summaries.append(store.read_aux(chunk_no, "seedsummary"))
    if nc.pooling == "pooled":
        n_null = np.where(np.isfinite(observed), n_null, 0)
    p = empirical_p(count_ge, n_null)
    seed_summary = (
        pd.concat(summaries, ignore_index=True)
        if summaries
        else pd.DataFrame(columns=["rep", "mean_S", "median_S", "seeds"])
    )
    return NullResult(
        nc.statistic, nc.pooling, nc.n_rep, work.genes, p, count_ge, n_null, seed_summary
    )
