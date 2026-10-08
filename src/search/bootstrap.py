"""Bootstrap or jackknife over samples of a per-gene statistic, chunked and resumable."""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from search.checkpoint import ChunkMeta, ChunkStore, chunk_ranges
from search.config import RunConfig
from search.core.corepass import CoreResult, resolve_seeds, run_core
from search.data import ExpressionSet
from search.parallel import parallel_map
from search.resampling import resample_indices

FloatArray = NDArray[np.float64]
log = logging.getLogger("search")


@dataclass
class BootstrapSummary:
    statistic: str
    scheme: str
    n_rep: int
    genes: pd.Index
    mean: FloatArray
    var: FloatArray
    sd: FloatArray
    n_ok: NDArray[np.int64]

    def table(self) -> pd.DataFrame:
        s = self.statistic
        return pd.DataFrame(
            {
                "gene": np.asarray(self.genes),
                f"{s}_boot_mean": self.mean,
                f"{s}_boot_var": self.var,
                f"{s}_boot_sd": self.sd,
                "boot_n_ok": self.n_ok,
            }
        )


def statistic_of(result: CoreResult, statistic: str) -> FloatArray:
    if statistic == "C_rank":
        if result.C_rank is None:
            raise ValueError("statistic C_rank needs connectivity.rank = fractional")
        return np.asarray(result.C_rank, dtype=np.float64)
    return np.asarray(result.C_raw, dtype=np.float64)


def stage_hash(cfg: RunConfig, *parts: str) -> str:
    text = "|".join([cfg.core_hash(), *parts])
    return hashlib.sha256(text.encode()).hexdigest()[:16]


@dataclass
class _Job:
    work: ExpressionSet
    seeds: list[str]
    cfg: RunConfig
    rep: int

    def run(self) -> tuple[int, NDArray[np.intp], FloatArray]:
        idx = resample_indices(
            self.cfg.resampling.scheme, self.work.n_samples, self.cfg.resampling.seed_base, self.rep
        )
        res = run_core(self.work, self.seeds, self.cfg, sample_idx=idx)
        return self.rep, idx, statistic_of(res, self.cfg.resampling.statistic)


def _run_job(job: _Job) -> tuple[int, NDArray[np.intp], FloatArray]:
    return job.run()


def run_bootstrap(
    expr: ExpressionSet,
    seeds: list[str],
    cfg: RunConfig,
    *,
    universe: list[str] | None,
    ckpt_dir: Path | str,
    input_hash: str,
    resume: bool = True,
    n_jobs: int = 1,
    store_indices: bool = True,
) -> BootstrapSummary:
    rc = cfg.resampling
    work, seed_list = resolve_seeds(expr, seeds, universe)
    store = ChunkStore(ckpt_dir, "bootstrap", cfg.checkpoint_dtype)
    config_hash = stage_hash(cfg, rc.scheme, rc.statistic)
    ranges = chunk_ranges(rc.n_rep, rc.chunk_size)
    for chunk_no, (start, end) in enumerate(ranges, start=1):
        meta = ChunkMeta(
            "bootstrap",
            start,
            end,
            rc.seed_base,
            config_hash,
            input_hash,
            {"n_samples": work.n_samples},
        )
        if resume and store.is_valid(chunk_no, meta):
            log.info("bootstrap chunk %d/%d already done, skipping", chunk_no, len(ranges))
            continue
        log.info("bootstrap chunk %d/%d, replicates %d-%d", chunk_no, len(ranges), start, end)
        jobs = [_Job(work, seed_list, cfg, rep) for rep in range(start, end + 1)]
        outputs = parallel_map(_run_job, jobs, n_jobs)
        columns = {rep: values for rep, _, values in outputs}
        aux = None
        if store_indices:
            idx_frame = pd.DataFrame({f"rep_{rep:05d}": idx for rep, idx, _ in outputs})
            aux = {"indices": idx_frame}
        store.write(chunk_no, work.genes, columns, meta, aux)
    if not ranges:
        empty = np.full(work.n_genes, np.nan)
        return BootstrapSummary(
            rc.statistic,
            rc.scheme,
            0,
            work.genes,
            empty,
            empty,
            empty,
            np.zeros(work.n_genes, dtype=np.int64),
        )
    matrix = store.read_all(len(ranges)).reindex(work.genes).to_numpy(dtype=np.float64)
    with np.errstate(invalid="ignore"):
        mean = np.nanmean(matrix, axis=1)
        var = (
            np.nanvar(matrix, axis=1, ddof=1)
            if matrix.shape[1] > 1
            else np.full(work.n_genes, np.nan)
        )
    n_ok = np.isfinite(matrix).sum(axis=1).astype(np.int64)
    return BootstrapSummary(
        rc.statistic, rc.scheme, rc.n_rep, work.genes, mean, var, np.sqrt(var), n_ok
    )
