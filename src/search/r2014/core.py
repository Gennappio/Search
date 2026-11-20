"""``SEARCH()`` of SEARCH 0.0.2 and its helpers, line by line, defects included.

This is not the engine of Search (:mod:`search.growth`), which computes the method as it was
meant. It is a transliteration of the R code as it was written, kept to answer one question:
can a result of the 2014 package be reproduced exactly? The R statement each block stands for
is quoted in a comment, and the places where the package misbehaves are marked ``R defect``
with the number they have in ``docs/porting.md``.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from search.r2014.association import AssociationData, nodes_correlation, universe_rows
from search.r2014.rsem import (
    NA_REAL,
    FloatArray,
    Name,
    NamedVector,
    Names,
    RError,
    RList,
    RMatrix,
    UnportedState,
    compare,
    concat_lists,
    get_by_names,
    intersect,
    is_in,
    remove_by_name_list,
    remove_by_name_vector,
    removeNode,
    set_by_names,
    set_names,
    setdiff,
    subset_logical,
    union,
)
from search.r2014.rstats import row_sd

AssociationFunction = Callable[..., AssociationData]


def _association_functions() -> dict[str, AssociationFunction]:
    from search.r2014.glm import nodes_fastLm, nodes_glm

    return {
        "nodes.correlation": nodes_correlation,
        "nodes.glm": nodes_glm,
        "nodes.fastLm": nodes_fastLm,
    }


@dataclass
class SearchResult:
    """The list ``SEARCH()`` returns. The dictionaries are keyed by the iteration number as a
    string, as the R lists are."""

    removed: Names
    correlationLists: RList
    correlationMatrix: dict[str, RMatrix]
    D: dict[str, NamedVector]
    O: dict[str, RMatrix]  # noqa: E741  (the name the package gives it)
    network_growth: dict[str, Names]
    other_data: dict[str, RMatrix] = field(default_factory=dict)
    nodes_degree: dict[str, NamedVector] = field(default_factory=dict)


# ---- overlap ----------------------------------------------------------------------------------
def getO(x: NamedVector, y: NamedVector) -> float:
    """``getO``: the members two neighbourhoods share over the members of either."""
    xn = set(x.names or [])
    yn = set(y.names or [])
    n_union = len(xn | yn)
    if n_union == 0:
        return float("nan")  # 0 / 0
    return len(xn & yn) / n_union


def _incidence(lists: Sequence[NamedVector]) -> tuple[FloatArray, FloatArray] | None:
    """Neighbourhoods as rows of a 0/1 matrix over their distinct member names, with the
    number of distinct members of each. ``None`` when a list has no names at all."""
    column: dict[Name, int] = {}
    rows: list[list[int]] = []
    for item in lists:
        if item.names is None:
            return None
        rows.append(sorted({column.setdefault(name, len(column)) for name in item.names}))
    m = np.zeros((len(lists), max(len(column), 1)), dtype=np.float64)
    for i, cols in enumerate(rows):
        m[i, cols] = 1.0
    return m, m.sum(axis=1)


def _overlap_block(xs: Sequence[NamedVector], ys: Sequence[NamedVector]) -> FloatArray:
    """``getO`` of every list of ``xs`` with every list of ``ys``. Counts are whole numbers,
    so the matrix product gives the same ratios as the pairwise function."""
    both = _incidence([*xs, *ys])
    if both is None:
        return np.array([[getO(a, b) for b in ys] for a in xs], dtype=np.float64).reshape(
            len(xs), len(ys)
        )
    m, size = both
    mx, my = m[: len(xs)], m[len(xs) :]
    inter = mx @ my.T
    n_union = size[: len(xs), None] + size[None, len(xs) :] - inter
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.asarray(inter / n_union, dtype=np.float64)


def getOmax(x: RList, y: RList) -> NamedVector | None:
    """``getOmax``: for every list of ``x``, its largest overlap with the lists of ``y``.

    R defect R2: list ``i`` of ``x`` is compared only with the lists of ``y`` whose position
    is greater than ``i``, so a list is never compared with those before it, and the last one
    of a symmetric call is compared with nothing and gets 0.
    """
    if len(x) == 0:
        return None  # c(), that is NULL
    overlap = _overlap_block(x.items, y.items)
    r = np.zeros(len(x))
    for i in range(len(x)):
        # r[x.index] <<- max(getO(x.i, y.j), r[x.index]) over the j after i, starting from 0;
        # max() of R is NaN as soon as an argument is
        later = overlap[i, i + 1 :]
        if later.size:
            r[i] = float("nan") if np.isnan(later).any() else max(float(later.max()), 0.0)
    return NamedVector(r, None)


def get_matrix_format(x: RList, all_features: Sequence[Name]) -> tuple[RMatrix, RMatrix]:
    """``get.matrix.format``: the O matrix of the lists and the matrix of their members'
    values. Returns ``(O, subnetworks.correlation)``."""
    a = x.names
    if a is None or len(a) == 0:
        raise UnportedState("get.matrix.format on a list without names")
    if any(name is None for name in a) or len(set(a)) != len(a):
        raise UnportedState("neighbourhood lists with missing or duplicated names")
    n = len(a)
    feature_pos: dict[Name, int] = {}
    for i, name in enumerate(all_features):
        feature_pos.setdefault(name, i)
    if len(feature_pos) != len(all_features):
        raise UnportedState("duplicated row names in the data matrix")

    # correlation.matrix[ a[i], names(x[[a[i]]]) ] <- x[[a[i]]];
    used: set[int] = set()
    entries: list[tuple[list[int], FloatArray]] = []
    for item in x.items:
        if item.names is None:
            entries.append(([], np.empty(0)))  # a zero-length index: nothing is assigned
            continue
        cols: list[int] = []
        for name in item.names:
            if name is None or name not in feature_pos:
                raise RError(
                    "subscript out of bounds",
                    "`[<-`(`*tmp*`, a[i], names(x[[a[i]]]), value = x[[a[i]]])",
                )
            cols.append(feature_pos[name])
        entries.append((cols, item.values))
        used.update(cols)

    # for(j in i:length(a)) { O[i,j] <- getO(x[[i]], x[[j]]); }: the upper triangle only
    overlap = np.triu(_overlap_block(x.items, x.items))
    o_matrix = RMatrix(
        np.where(np.triu(np.ones((n, n), dtype=bool)), overlap, 0.0), list(a), list(a)
    )

    # the columns that are not all NA; a value of NA does not count as a member
    kept = sorted(used)
    column_of = {col: k for k, col in enumerate(kept)}
    values = np.full((n, len(kept)), NA_REAL)
    for i, (cols, vals) in enumerate(entries):
        if len(set(cols)) == len(cols):
            values[i, [column_of[col] for col in cols]] = vals
        else:
            for col, value in zip(cols, vals, strict=True):
                values[i, column_of[col]] = value  # a name given twice: the last value stays
    present = ~np.isnan(values).all(axis=0) if len(kept) else np.zeros(0, dtype=bool)
    values = values[:, present]
    names: Names = [all_features[col] for col, keep in zip(kept, present, strict=True) if keep]
    if values.shape[1] == 1:
        # one column left: the data frame gives it back as a plain vector and as.matrix()
        # makes a column without any name of it
        return o_matrix, RMatrix(values, None, None)
    if values.shape[1] == 0:
        return o_matrix, RMatrix(values, list(a), None)
    return o_matrix, RMatrix(values, list(a), names)


def estimate_nodes_degree(correlation_matrix_list: dict[str, RMatrix]) -> dict[str, NamedVector]:
    """``estimate.nodes.degree``: for every seed, the share of the seeds found among the
    members of its neighbourhood."""
    out: dict[str, NamedVector] = {}
    for iter_i, cm in correlation_matrix_list.items():
        common = intersect(cm.rownames or [], cm.colnames or [])
        if len(common) == 0:
            out[iter_i] = NamedVector(np.empty(0), None)
            continue
        assert cm.rownames is not None and cm.colnames is not None
        rows = [cm.rownames.index(name) for name in common]
        cols = [cm.colnames.index(name) for name in common]
        sub = cm.values[np.ix_(rows, cols)]
        # length(which(x != "NA")): NA is not counted, NaN is
        from search.r2014.rsem import is_na_real

        degree = (~is_na_real(sub.ravel()).reshape(sub.shape)).sum(axis=1) / sub.shape[1]
        # a single common node: the 1 x 1 selection loses its names
        out[iter_i] = NamedVector(degree, list(common) if len(common) > 1 else None)
    return out


# ---- SEARCH -----------------------------------------------------------------------------------
def SEARCH(
    x: RMatrix,
    Alist: Sequence[Name],
    Blist: Sequence[Name],
    covariates: Any = None,
    scale: bool = True,
    D_lower: float = 0.0,
    D_higher: float = 0.0,
    O_lower: float = 0.0,
    keep_Alist: bool = True,
    max_iter: int = 1,
    sd_threshold: float = 0.0,
    verbose: bool = True,
    nodes_association_function: str | AssociationFunction = "nodes.correlation",
    **association_args: Any,
) -> SearchResult:
    """``SEARCH(x, Alist, Blist, ...)``.

    ``association_args`` are the arguments R passes through ``...`` to the association
    function, with dots written as underscores (``correlation_threshold = 0.5``).
    ``verbose`` only controls what R prints and is ignored.
    """
    alist: Names = list(Alist)
    if len(set(alist)) != len(alist):
        raise UnportedState("duplicated seeds")
    correlation_lists = RList([], None)
    removed: Names = []
    added: Names = []
    network_growth: dict[str, Names] = {}
    network_connectivity: dict[str, NamedVector] = {}
    network_correlation: dict[str, RMatrix] = {}
    network_overlap: dict[str, RMatrix] = {}
    other_data: dict[str, RMatrix] = {}
    association_data: AssociationData | None = None

    # Blist <- union(Alist, Blist);
    blist: Names = union(alist, list(Blist))

    if sd_threshold > 0:
        # sds <- apply(x[rownames(x) %in% Blist,], 1, sd);
        rows = universe_rows(x, blist)
        if len(rows) == 1:
            raise RError(
                "dim(X) must have a positive length", "apply(x[rownames(x) %in% Blist, ], 1, sd)"
            )
        sds = NamedVector(row_sd(x.values[rows]), None)
        # names(sds) <- Blist;
        # R defect R1: the standard deviations are in the row order of the matrix and get the
        # names of the universe in the order of the universe, seeds first.
        named = set_names(sds, blist)
        assert named is not None
        # to.remove <- names(sds[sds < sd.threshold]);
        to_remove = subset_logical(named, compare(named.values, "<", sd_threshold)).names or []
        alist = removeNode(to_remove, alist)
        blist = removeNode(to_remove, blist)

    iteration = 1
    while iteration <= max_iter:
        n_rem_o = 0

        if len(alist) == 0:
            correlation_lists = RList([], None)
            break

        if isinstance(nodes_association_function, str):
            functions = _association_functions()
            if nodes_association_function not in functions:
                raise RError(
                    f"object '{nodes_association_function}' of mode 'function' was not found",
                    'get(nodes.association.function, mode = "function")',
                )
            dynamic_function = functions[nodes_association_function]
        else:
            dynamic_function = nodes_association_function

        if iteration == 1:
            association_data = dynamic_function(
                seeds=alist, x=x, Blist=blist, covariates=covariates, **association_args
            )
            correlation_lists = association_data.CL
        elif len(added) > 0:
            # just do the new stuff
            association_data = dynamic_function(
                seeds=added, x=x, Blist=blist, covariates=covariates, **association_args
            )
            new_cl = association_data.CL
            if not keep_Alist:
                # Omax <- getOmax(x = c(newCL, correlationLists), y = c(newCL, correlationLists));
                # names(Omax) <- c(added, names(correlationLists));
                both = concat_lists(new_cl, correlation_lists)
                omax = set_names(getOmax(both, both), [*added, *(correlation_lists.names or [])])
            else:
                # Omax <- getOmax(newCL, c(newCL, correlationLists)); names(Omax) <- added;
                omax = set_names(getOmax(new_cl, concat_lists(new_cl, correlation_lists)), added)

            # to.remove <- names(to.remove[to.remove < O.lower]);
            to_remove = (
                []
                if omax is None
                else subset_logical(omax, compare(omax.values, "<", O_lower)).names or []
            )
            alist = removeNode(to_remove, alist)
            blist = removeNode(to_remove, blist)
            # newCL <- removeByName(to.remove, newCL);
            # R defect R3: only the new lists are filtered. With keep.Alist = FALSE an old seed
            # removed here leaves the seeds and the universe and keeps its list.
            new_cl = remove_by_name_list(to_remove, new_cl)
            n_rem_o = len(to_remove)
            correlation_lists = concat_lists(new_cl, correlation_lists)
        assert association_data is not None

        # counts <- rep(0, length(Blist)); names(counts) <- Blist;
        position = {name: i for i, name in enumerate(blist)}
        counts = np.zeros(len(blist))
        if len(correlation_lists) == 0:
            # for(i in 1:length(toCount)): 1:0 is c(1, 0), and toCount[[1]] does not exist
            raise RError("subscript out of bounds", "toCount[[i]]")
        for item in correlation_lists.items:
            # all.names %in% names(tmp.x): a list without names contains nothing
            members = {position[name] for name in set(item.names or []) if name in position}
            if members:
                counts[sorted(members)] += 1.0
        d = NamedVector(counts, list(blist))

        # counts[names(toCount)] <- counts[names(toCount)] - 1;
        # A list name that is no longer in the universe is not an error here: reading it gives
        # NA, and assigning it appends an element with that name (R defect R3).
        list_names: Names = list(correlation_lists.names or [])
        d = set_by_names(d, list_names, get_by_names(d, list_names).values - 1.0)
        if scale:
            with np.errstate(divide="ignore", invalid="ignore"):
                # counts[names(toCount)] <- counts[names(toCount)]/(length(Alist)-1);
                d = set_by_names(
                    d, list_names, get_by_names(d, list_names).values / float(len(alist) - 1)
                )
                # counts[setdiff(Blist, names(toCount))] <- .../length(Alist);
                others = setdiff(blist, list_names)
                d = set_by_names(d, others, get_by_names(d, others).values / float(len(alist)))

        # to.remove <- D[names(D) %in% Alist]; to.remove <- names(to.remove[to.remove < D.lower]);
        assert d.names is not None
        of_seeds = np.array(is_in(d.names, alist), dtype=bool)
        seeds_d = NamedVector(
            d.values[of_seeds], [n for n, k in zip(d.names, of_seeds, strict=True) if k]
        )
        to_remove = subset_logical(seeds_d, compare(seeds_d.values, "<", D_lower)).names or []
        alist = removeNode(to_remove, alist)
        blist = removeNode(to_remove, blist)
        # correlationLists <- removeByName(to.remove, correlationLists);
        # R defect R4: the lists of the removed seeds go, but the removed genes stay as members
        # of the other lists.
        correlation_lists = remove_by_name_list(to_remove, correlation_lists)
        d = remove_by_name_vector(to_remove, d)
        network_connectivity[str(iteration)] = d
        removed = [*to_remove, *removed]
        n_rem_d = len(to_remove)

        # not.already <- !(names(D) %in% Alist); to.add <- D[not.already];
        # to.add <- names(to.add[to.add > D.higher]); Alist <- c(Alist, to.add);
        assert d.names is not None
        not_already = ~np.array(is_in(d.names, alist), dtype=bool)
        candidates = NamedVector(
            d.values[not_already], [n for n, k in zip(d.names, not_already, strict=True) if k]
        )
        to_add = subset_logical(candidates, compare(candidates.values, ">", D_higher)).names or []
        alist = [*alist, *to_add]
        n_add_d = len(to_add)
        added = to_add
        network_growth[str(iteration)] = list(to_add)

        # other.data[[as.character(iter)]] <- association.data[["other.data"]];
        if association_data.other_data is not None:
            other_data[str(iteration)] = association_data.other_data

        if (n_add_d + n_rem_d + n_rem_o) == 0 or len(blist) == 0 or len(correlation_lists) == 0:
            break

        # the matrices are stored only when another iteration may follow
        o_matrix, correlation_matrix = get_matrix_format(correlation_lists, x.rownames or [])
        network_correlation[str(iteration)] = correlation_matrix
        network_overlap[str(iteration)] = o_matrix
        iteration += 1

    return SearchResult(
        removed=removed,
        correlationLists=correlation_lists,
        correlationMatrix=network_correlation,
        D=network_connectivity,
        O=network_overlap,
        network_growth=network_growth,
        other_data=other_data,
    )
