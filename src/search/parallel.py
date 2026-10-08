"""Replicate-level parallelism. Sequential when n_jobs is 1; joblib processes otherwise, with
BLAS pinned to one thread per worker so that workers do not oversubscribe the cores."""

from __future__ import annotations

import os
from collections.abc import Callable, Iterable


def default_n_jobs() -> int:
    value = os.environ.get("SLURM_CPUS_PER_TASK")
    if value and value.isdigit():
        return max(1, int(value))
    return 1


def _single_thread_call[T, R](fn: Callable[[T], R], item: T) -> R:
    from threadpoolctl import threadpool_limits

    with threadpool_limits(limits=1):
        return fn(item)


def parallel_map[T, R](fn: Callable[[T], R], items: Iterable[T], n_jobs: int = 1) -> list[R]:
    """Apply ``fn`` to every item, in order; processes when n_jobs > 1."""
    items = list(items)
    if n_jobs <= 1 or len(items) <= 1:
        return [fn(item) for item in items]
    from joblib import Parallel, delayed

    results = Parallel(n_jobs=n_jobs, prefer="processes")(
        delayed(_single_thread_call)(fn, item) for item in items
    )
    return list(results)
