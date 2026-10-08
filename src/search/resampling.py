"""Reproducible resampling: every replicate derives its own generator from a base seed, so
replicates are independent of chunking, parallelism and execution order."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from search.data import ExpressionSet
from search.stats import row_sds


def derive_rng(seed_base: int, rep: int) -> np.random.Generator:
    """Generator of replicate ``rep`` (1-based)."""
    return np.random.default_rng(seed_base + rep)


def bootstrap_indices(n_samples: int, rng: np.random.Generator) -> NDArray[np.intp]:
    """n draws with replacement."""
    return np.asarray(rng.integers(0, n_samples, size=n_samples), dtype=np.intp)


def jackknife_indices(n_samples: int, rng: np.random.Generator) -> NDArray[np.intp]:
    """n - 1 distinct samples in random order."""
    return np.asarray(rng.permutation(n_samples)[: n_samples - 1], dtype=np.intp)


def resample_indices(scheme: str, n_samples: int, seed_base: int, rep: int) -> NDArray[np.intp]:
    rng = derive_rng(seed_base, rep)
    if scheme == "bootstrap":
        return bootstrap_indices(n_samples, rng)
    if scheme == "jackknife":
        return jackknife_indices(n_samples, rng)
    raise ValueError(f"unknown resampling scheme {scheme}")


def eligible_null_genes(
    expr: ExpressionSet, seeds: list[str], min_sd: float = 0.0
) -> NDArray[np.str_]:
    """Genes that may serve as random seeds: variable enough and not real seeds."""
    sds = row_sds(expr.values)
    keep = sds > min_sd
    seed_set = set(seeds)
    names = np.asarray([str(g) for g in expr.genes])
    mask = keep & np.array([g not in seed_set for g in names])
    return names[mask]


def random_seed_set(eligible: NDArray[np.str_], k: int, rng: np.random.Generator) -> list[str]:
    if len(eligible) < k:
        raise ValueError(f"not enough eligible genes ({len(eligible)}) for {k} random seeds")
    return [str(g) for g in rng.choice(eligible, size=k, replace=False)]
