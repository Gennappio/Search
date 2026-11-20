"""``summarise.iterations`` of SEARCH 0.0.2.

It reads the connectivity files written by ``SEARCH.randomise``: those of a run with the real
seeds, and those of runs with random seeds that the user made by hand into folders named
``batch_<i>_0_FALSE_0_0``. For every growth iteration it averages D over the replicates and
gives every gene the share of the null values of that gene that reach the average.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from search.r2014 import platform, x87
from search.r2014.rio import read_connectivity
from search.r2014.rsem import NA_REAL, FloatArray, Names, RError, RMatrix, union


@dataclass
class IterationSummary:
    """The list ``summarise.iterations`` returns: ``real`` and ``NULL``, each keyed by the
    iteration number as a string."""

    real: dict[str, RMatrix] = field(default_factory=dict)
    null: dict[str, RMatrix] = field(default_factory=dict)


def _colon(start: int, end: int | None) -> list[int]:
    """``start:end`` of R: it counts down when ``end`` is smaller, so ``1:0`` is ``c(1, 0)``."""
    if end is None:
        raise RError("argument of length 0", f"{start}:NULL")
    step = 1 if end >= start else -1
    return list(range(start, end + step, step))


def _row_means(values: FloatArray) -> FloatArray:
    """``rowMeans(values, na.rm = TRUE)``: the columns are added from the first to the last in
    the accumulator of the platform, a double on arm64 and the 80-bit type on x86, and the sum
    is divided by the number of values. A row without values gives 0 / 0.

    The accumulator matters: the mean of ten copies of 0.1 is 0.1 in 80 bits and
    0.09999999999999999 in doubles, and the package then asks whether a null value of 0.1 is
    at least as large as that mean.
    """
    if platform.current().extended:
        out = np.empty(values.shape[0])
        for i, row in enumerate(values):
            kept = [float(v) for v in row if not math.isnan(v)]
            if not kept:
                out[i] = float("nan")
            elif not all(math.isfinite(v) for v in kept):
                out[i] = sum(kept)  # an infinite D: the mean is it, or NaN for both signs
            else:
                total_ext = x87.ZERO
                for v in kept:
                    total_ext = x87.add(total_ext, x87.from_double(v))
                out[i] = x87.to_double(x87.div(total_ext, x87.from_int(len(kept))))
        return out
    total = np.zeros(values.shape[0])
    count = np.zeros(values.shape[0])
    for j in range(values.shape[1]):
        column = values[:, j]
        present = ~np.isnan(column)
        total = np.where(present, total + np.where(present, column, 0.0), total)
        count = count + present
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.asarray(total / count, dtype=np.float64)


def _collect(files: dict[str, RMatrix], features: Names) -> RMatrix:
    """The matrix of the replicates read so far::

        m <- matrix(NA, nrow = length(features.all), ncol = length(l),
                    dimnames = list(features.all, 1:length(l)));
        for (i in colnames(m)) m[rownames(l[[i]]), i] <- l[[i]][, 1];

    R defect R9: the columns are numbered 1, 2, ... up to the number of files found, and each
    is filled from the file of the replicate with that number. When a replicate in the middle
    has no file (it stopped at an earlier iteration), the replicates after it are shifted out:
    their files are read and never used, and the column of the missing one stays empty.
    """
    n_files = len(files)
    if n_files == 0:
        # 1:0 is c(1, 0): two names for no column
        raise RError("length of 'dimnames' [2] not equal to array extent")
    row_of = {name: i for i, name in enumerate(features)}
    values = np.full((len(features), n_files), NA_REAL)
    for j in range(n_files):
        table = files.get(str(j + 1))
        if table is None or table.rownames is None:
            continue
        for name, value in zip(table.rownames, table.values[:, 0], strict=True):
            if name not in row_of:
                # R defect R8: a gene of a null run that no real run has
                raise RError("subscript out of bounds")
            values[row_of[name], j] = value
    # no gene at all: the row names are those of an empty vector, that is none
    return RMatrix(
        values, list(features) if features else None, [str(j + 1) for j in range(n_files)]
    )


def summarise_iterations(
    output_directory: str = ".",
    max_iter: int = 1,
    n_boot: int = 1,
    output_directory_NULL: str | None = None,
    n_boot_NULL: int | None = None,
    batches_NULL: int | None = None,
) -> IterationSummary:
    """``summarise.iterations``. The directories are pasted in front of the file names as they
    are, so they must end with a slash."""
    features_all: Names = []
    out = IterationSummary()
    null_root = "" if output_directory_NULL is None else output_directory_NULL

    for iter_i in (str(k) for k in _colon(1, max_iter)):
        # real runs
        conn_all: dict[str, RMatrix] = {}
        for i in (str(k) for k in _colon(1, n_boot)):
            table = read_connectivity(
                f"{output_directory}output/connectivity__iteration_{i}_{iter_i}.txt.gz"
            )
            if table is not None:
                conn_all[i] = table
                features_all = union(features_all, table.rownames or [])
        real = _collect(conn_all, features_all)

        # conn.score.mean[[iter]] <- cbind(..., "mean" = rowMeans(..., na.rm = TRUE));
        mean = _row_means(real.values)
        # ordered by the mean, largest first; rows without a mean last; ties keep their order
        order = np.argsort(-mean, kind="stable")
        rows = [features_all[k] for k in order]
        mean_sorted = mean[order]
        real_sorted = real.values[order]

        # null runs
        null_all: dict[str, RMatrix] = {}
        for batch in _colon(1, batches_NULL):
            for replicate in _colon(1, n_boot_NULL):
                table = read_connectivity(
                    f"{null_root}batch_{batch}_0_FALSE_0_0/output/"
                    f"connectivity__iteration_{replicate}_{iter_i}.txt.gz"
                )
                if table is not None:
                    assert n_boot_NULL is not None
                    null_all[str(replicate + n_boot_NULL * (batch - 1))] = table

        if real.nrow == 1:
            # one gene: m[order(...), ] is no longer a matrix, and nrow() of it is NULL
            raise RError("invalid 'times' argument")
        p_value = np.full(len(rows), NA_REAL)
        if len(null_all) > 0:
            null = _collect_null(null_all, rows)
            # (number of null values >= the mean, + 1) / (number of null columns + 1): a gene
            # absent from a null run is not counted above and is counted below
            with np.errstate(invalid="ignore"):
                reached = (null.values >= mean_sorted[:, None]).sum(axis=1)
            p_value = (reached + 1) / (null.ncol + 1)
            out.null[iter_i] = null

        out.real[iter_i] = RMatrix(
            np.column_stack([real_sorted, mean_sorted, p_value])
            if len(rows)
            else np.empty((0, real.ncol + 2)),
            list(rows) if rows else None,
            [*(real.colnames or []), "mean", "P"],
        )
    return out


def _collect_null(files: dict[str, RMatrix], rows: Names) -> RMatrix:
    """The null matrix, with its rows in the order of the sorted real one. With a single null
    file the reordering ``m[rownames(real), ]`` drops the matrix to a vector, and the function
    stops at the first gene."""
    if len(files) == 1:
        raise RError("incorrect number of dimensions")
    return _collect(files, rows)
