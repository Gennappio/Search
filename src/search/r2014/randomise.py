"""``SEARCH.randomise``, ``.SEARCH.randomise.inner`` and ``write.data.files`` of SEARCH 0.0.2.

``SEARCH.randomise`` is the entry point of the package: it resamples the columns of the data
``n.boot`` times, runs ``SEARCH`` on every resample and writes the results of each to
``<output.directory>/output/``. It never runs ``SEARCH`` on the data as they are.
"""

from __future__ import annotations

import gzip
import os
import random
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np

from search.r2014 import platform
from search.r2014.core import (
    SEARCH,
    AssociationFunction,
    SearchResult,
    estimate_nodes_degree,
)
from search.r2014.rio import table_of_matrix, table_of_strings, table_of_vector, write_and_gzip
from search.r2014.rng import RRandom, SampleKind
from search.r2014.rsem import (
    NA_REAL,
    FloatArray,
    Name,
    NamedVector,
    RError,
    RMatrix,
    intersect,
    is_na_real,
)
from search.r2014.rstats import r_round

Summary = dict[str, Any]


def _three_decimals(values: FloatArray) -> FloatArray:
    """``as.numeric(sprintf("%.3f", v))``: through text and back."""
    out = np.empty(values.shape, dtype=np.float64)
    na = is_na_real(values.ravel()).reshape(values.shape)
    for index, v in np.ndenumerate(values):
        out[index] = NA_REAL if na[index] else float(f"{float(v):.3f}")
    return out


def _through_text(m: RMatrix) -> RMatrix:
    """The round trip ``write.data.files`` gives a matrix before writing it::

        m <- matrix(sprintf("%.3f", m), nrow = nrow(m), dimnames = list(rownames(m), colnames(m)));
        m <- t(apply(m, 1, as.numeric)); rownames(m) <- tmp.rownames; colnames(m) <- tmp.colnames;

    ``apply`` returns a matrix to transpose only when there are at least two columns. With one
    column it returns a plain vector and ``t()`` makes a single row of it; with none it returns
    an empty vector. Putting the old row names on what comes out then fails unless there was
    one row.
    """
    values = _three_decimals(m.values)
    if m.ncol >= 2:
        return RMatrix(values, m.rownames, m.colnames)
    if m.ncol == 1:
        if m.rownames is not None and m.nrow != 1:
            raise RError("length of 'dimnames' [1] not equal to array extent", "dimnames(x) <- dn")
        if m.colnames is not None and m.nrow != 1:
            raise RError("length of 'dimnames' [2] not equal to array extent", "dimnames(x) <- dn")
        return RMatrix(values.reshape(1, m.nrow), m.rownames, m.colnames)
    if m.rownames is not None and m.nrow != 1:
        raise RError("length of 'dimnames' [1] not equal to array extent", "dimnames(x) <- dn")
    return RMatrix(np.empty((1, 0)), m.rownames, m.colnames)


def _sort_decreasing(vec: NamedVector) -> NamedVector:
    """``sort(vec, decreasing = TRUE)``: missing values are dropped, ties keep their order."""
    keep = ~np.isnan(vec.values)
    values = vec.values[keep]
    names = None if vec.names is None else [n for n, k in zip(vec.names, keep, strict=True) if k]
    order = np.argsort(-values, kind="stable")
    return NamedVector(values[order], None if names is None else [names[i] for i in order])


def write_data_files(results: SearchResult, data_directory: str, iteration: int) -> SearchResult:
    """``write.data.files``: writes the files of one replicate and returns the results as
    changed on the way (matrices through three decimals, D rounded to five and sorted)."""

    def path(kind: str, iter_i: str) -> str:
        return f"{data_directory}{kind}__iteration_{iteration}_{iter_i}.txt"

    for iter_i in list(results.correlationMatrix):
        results.correlationMatrix[iter_i] = _through_text(results.correlationMatrix[iter_i])
        write_and_gzip(
            table_of_matrix(results.correlationMatrix[iter_i]),
            path("subnetworks_correlation", iter_i),
        )

    for iter_i in list(results.O):
        results.O[iter_i] = _through_text(results.O[iter_i])
        write_and_gzip(table_of_matrix(results.O[iter_i]), path("shared_neighbourhood", iter_i))

    for iter_i in list(results.D):
        d = results.D[iter_i]
        results.D[iter_i] = _sort_decreasing(NamedVector(r_round(d.values, 5), d.names))
        write_and_gzip(table_of_vector(results.D[iter_i]), path("connectivity", iter_i))

    for iter_i, added in results.network_growth.items():
        write_and_gzip(table_of_strings(added), path("network_growth", iter_i))

    for iter_i in list(results.nodes_degree):
        nd = results.nodes_degree[iter_i]
        results.nodes_degree[iter_i] = NamedVector(r_round(nd.values, 5), nd.names)
        write_and_gzip(
            table_of_vector(_sort_decreasing(results.nodes_degree[iter_i])),
            path("nodes_degree", iter_i),
        )

    for iter_i, table in results.other_data.items():
        write_and_gzip(table_of_matrix(table), path("other_data", iter_i))

    return results


def _summary(results: SearchResult) -> Summary:
    """What ``.SEARCH.randomise.inner`` returns: sizes only, everything else is in the files."""
    return {
        "removed": [len(results.removed)],
        "correlationLists": [len(results.correlationLists)],
        "correlationMatrix": {k: [m.nrow, m.ncol] for k, m in results.correlationMatrix.items()},
        "D": {k: [len(v)] for k, v in results.D.items()},
        "O": {k: [m.nrow, m.ncol] for k, m in results.O.items()},
        "network.growth": {k: [len(v)] for k, v in results.network_growth.items()},
    }


def resample_columns(
    n: int,
    iteration: int,
    sampling: str,
    sampling_seed: float | str | None,
    sample_kind: SampleKind,
    rng: RRandom | None,
) -> np.ndarray:
    """The columns ``.SEARCH.randomise.inner`` draws for one replicate, 1-based as in R."""
    if sampling_seed is not None:
        try:
            seed = float(sampling_seed)
        except ValueError:
            raise RError("supplied seed is not a valid integer", "set.seed") from None
        seed = iteration + seed - 1
        if seed != seed or abs(seed) > 2147483647:
            raise RError("supplied seed is not a valid integer", "set.seed")
        rng = RRandom(int(seed), sample_kind)
    elif rng is None:
        # no seed: R continues the stream of the session, which nothing here can know
        rng = RRandom(random.SystemRandom().randrange(2**31), sample_kind)
    if sampling == "bootstrap":
        return rng.sample(n, replace=True)
    if sampling == "jacknife":
        return rng.sample(n, size=n - 1, replace=False)
    raise RError("\nQuitting. unknown sampling method")


def SEARCH_randomise_inner(
    iteration: int,
    x: RMatrix,
    Alist: Sequence[Name],
    Blist: Sequence[Name],
    covariates: Any,
    sampling: str,
    sampling_seed: float | str | None,
    output_directory: str,
    sample_kind: SampleKind = "Rejection",
    rng: RRandom | None = None,
    **search_args: Any,
) -> Summary:
    """``.SEARCH.randomise.inner``: one replicate."""
    out_dir = f"{output_directory}/output/"
    graphs_dir = f"{output_directory}/graphs/"
    logs_dir = f"{output_directory}/logs/"
    if not os.path.exists(out_dir):
        for directory in (out_dir, graphs_dir, logs_dir):
            os.makedirs(directory, exist_ok=True)

    random_samples = resample_columns(x.ncol, iteration, sampling, sampling_seed, sample_kind, rng)
    columns = random_samples - 1
    resampled = RMatrix(
        np.ascontiguousarray(x.values[:, columns]),
        x.rownames,
        None if x.colnames is None else [x.colnames[j] for j in columns],
    )
    covariates_sample = None if covariates is None else covariates.rows(columns)

    results = SEARCH(
        resampled, Alist=Alist, Blist=Blist, covariates=covariates_sample, **search_args
    )
    results.nodes_degree = estimate_nodes_degree(results.correlationMatrix)
    results = write_data_files(results, out_dir, iteration)
    # make.default.plots() would run here: the figures are not ported
    return _summary(results)


def SEARCH_randomise(
    data: RMatrix,
    seed_nodes: Sequence[Name],
    universe_nodes: Sequence[Name],
    covariates: Any = None,
    scale: bool = True,
    D_lower: float = 0.1,
    D_higher: float = 0.24,
    O_lower: float = 0.01,
    keep_seeds: bool = True,
    max_iter: int = 1,
    verbose: bool = True,
    sd_threshold: float = 0.1,
    n_boot: int = 10,
    n_boot_start: int = 1,
    sampling: str = "bootstrap",
    sampling_seed: float | str | None = None,
    output_directory: str = ".",
    nodes_association_function: str | AssociationFunction = "nodes.correlation",
    sample_kind: SampleKind | None = None,
    **association_args: Any,
) -> list[Summary]:
    """``SEARCH.randomise``.

    ``sample_kind`` is not an argument of the package: it says which ``sample()`` of R to
    reproduce, the one of R 3.6.0 and later (``"Rejection"``) or the one of the R the package
    was written for (``"Rounding"``). By default it is the one of the current platform
    (:mod:`search.r2014.platform`). Without ``sampling_seed`` the resamples are random and not
    reproducible, as in R.
    """
    if sample_kind is None:
        sample_kind = platform.current().sample_kind
    if data.nrow < 2 or data.ncol < 2:
        raise RError("cant proceed, data is either not a matrix OR does not contain data")
    # restrict seed.nodes to features present in expression data
    seeds = intersect(list(seed_nodes), data.rownames or [])
    if len(seeds) == 0:
        raise RError("cant proceed, seed.nodes and rownames in data have no overlap")

    # i <- n.boot.start:(n.boot.start+n.boot-1): with n.boot = 0 the sequence runs backwards
    last = n_boot_start + n_boot - 1
    step = 1 if last >= n_boot_start else -1
    rng = (
        None
        if sampling_seed is not None
        else RRandom(random.SystemRandom().randrange(2**31), sample_kind)
    )
    summaries: list[Summary] = []
    for iteration in range(n_boot_start, last + step, step):
        summaries.append(
            SEARCH_randomise_inner(
                iteration,
                data,
                Alist=seeds,
                Blist=universe_nodes,
                covariates=covariates,
                sampling=sampling,
                sampling_seed=sampling_seed,
                output_directory=output_directory,
                sample_kind=sample_kind,
                rng=rng,
                scale=scale,
                D_lower=D_lower,
                D_higher=D_higher,
                O_lower=O_lower,
                keep_Alist=keep_seeds,
                max_iter=max_iter,
                verbose=verbose,
                sd_threshold=sd_threshold,
                nodes_association_function=nodes_association_function,
                **association_args,
            )
        )
    return summaries


def read_output_files(root: str | Path) -> dict[str, str]:
    """Every file under ``root`` as text, compressed files uncompressed: relative path to
    content. For comparing what two runs wrote."""
    base = Path(root)
    out: dict[str, str] = {}
    for path in sorted(p for p in base.rglob("*") if p.is_file()):
        raw = path.read_bytes()
        out[path.relative_to(base).as_posix()] = (
            gzip.decompress(raw) if path.suffix == ".gz" else raw
        ).decode("utf-8")
    return out
