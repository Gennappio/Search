"""``nodes.correlation`` of SEARCH 0.0.2, line by line.

The association functions of the package take the seeds, the data matrix and the universe and
return, for every seed, its neighbourhood: a named vector with the genes associated with the
seed and the value of the association. ``nodes.glm`` and ``nodes.fastLm`` are in
:mod:`search.r2014.glm`.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from search.r2014.rsem import (
    Name,
    NamedVector,
    RError,
    RList,
    RMatrix,
    UnportedState,
    compare,
    subset_logical,
)
from search.r2014.rstats import r_cor


@dataclass
class AssociationData:
    """What an association function returns: ``CL``, the neighbourhoods, and ``other.data``
    (the table of all coefficients for the regression functions, ``NULL`` for correlation)."""

    CL: RList
    other_data: RMatrix | None = None


def universe_rows(x: RMatrix, Blist: Sequence[Name]) -> list[int]:
    """``rownames(x) %in% Blist``: the rows of the matrix inside the universe, in matrix
    order."""
    if x.rownames is None:
        raise UnportedState("a data matrix without row names")
    members = set(Blist)
    return [i for i, name in enumerate(x.rownames) if name in members]


def nodes_correlation(
    seeds: Sequence[Name],
    x: RMatrix,
    Blist: Sequence[Name],
    covariates: Any = None,
    correlation_threshold: float = 0.5,
    correlation_method: str = "spearman",
    correlation_type: str = "standard",
    correlation_use: str = "everything",
) -> AssociationData:
    """``nodes.correlation``. ``covariates`` is accepted and ignored, as in R."""
    if correlation_type != "standard":
        raise RError("correlation.type not supported!")

    # x.slimmed <- t(x[rownames(x) %in% Blist, ]);
    rows = universe_rows(x, Blist)
    assert x.rownames is not None
    if len(rows) < 2:
        # One row: R drops the matrix to a vector, t() makes a single row of it, and the
        # seeds are then looked up among the sample names. No row: nothing to look up in.
        sample_names = set(x.colnames or [])
        if len(rows) == 1 and len(seeds) > 0 and all(s in sample_names for s in seeds):
            raise UnportedState("a universe of one gene whose seeds are named like samples")
        raise RError("subscript out of bounds", "x.slimmed[, seeds]")
    genes = [x.rownames[i] for i in rows]
    position: dict[Name, int] = {}
    for i, name in enumerate(genes):
        position.setdefault(name, i)
    if len(position) != len(genes):
        raise UnportedState("duplicated row names in the data matrix")

    # x.cor <- cor(x.slimmed[, seeds], x.slimmed, method = correlation.method,
    #              use = correlation.use);
    seed_pos: list[int] = []
    for seed in seeds:
        if seed is None or seed not in position:
            raise RError("subscript out of bounds", "x.slimmed[, seeds]")
        seed_pos.append(position[seed])
    values = x.values[rows]
    x_cor = r_cor(values[seed_pos], values, correlation_method, correlation_use)

    # x.cor.filtered <- lapply(x.cor.list, function(q) { q[abs(q) > correlation.threshold] } )
    # R defect R7: where the correlation is missing the condition is NA, and q[NA] is not an
    # element left out but a missing value under a missing name.
    items: list[NamedVector] = []
    for k in range(len(seed_pos)):
        row = x_cor[k]
        if np.isnan(row).any() or correlation_threshold != correlation_threshold:
            q = NamedVector(row, list(genes))
            items.append(subset_logical(q, compare(np.abs(row), ">", correlation_threshold)))
        else:  # nothing missing: the same selection, without the element-by-element walk
            chosen = np.flatnonzero(np.abs(row) > correlation_threshold)
            items.append(NamedVector(row[chosen], [genes[i] for i in chosen]))
    return AssociationData(RList(items, list(seeds)), None)
