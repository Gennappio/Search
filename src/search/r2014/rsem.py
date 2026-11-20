"""The parts of R's semantics that SEARCH 0.0.2 leans on.

Most defects of the package are not arithmetic: they come from what R does with a named vector
when a name is missing, with a logical index that holds ``NA``, or with a list element set to
``NULL``. A port that reproduces the defects has to reproduce those rules, so they are written
here once, each next to the R behaviour it stands for, and the functions of the package are
transliterated on top of them.

Conventions:

- a missing string (``NA_character_``) is ``None``;
- a missing number (``NA_real_``) is the NaN R uses for it, with payload 1954, so that it
  survives arithmetic as it does in R and can be told from an ordinary ``NaN`` when printed;
- a vector without the ``names`` attribute has ``names = None``, which is not the same as a
  vector whose names are all missing.
"""

from __future__ import annotations

import struct
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from search.errors import SearchError

FloatArray = NDArray[np.float64]
Name = str | None
Names = list[Name]

#: ``NA_real_``: a quiet NaN whose low word is 1954 (``R_NaReal`` in ``arithmetic.c``).
NA_REAL: float = struct.unpack("<d", struct.pack("<Q", 0x7FF80000000007A2))[0]


class RError(SearchError):
    """The point where the R code stops with an error.

    ``message`` is what ``conditionMessage`` returns in R (the text after ``Error in ... :``),
    ``call`` the expression R reports, when the port knows it.
    """

    def __init__(self, message: str, call: str | None = None) -> None:
        super().__init__(message if call is None else f"Error in {call} : {message}")
        self.message = message
        self.call = call


class UnportedState(SearchError):
    """A state the package can reach only from inputs this port does not model (duplicated
    row names, duplicated seeds given to ``SEARCH`` directly, and the like). The port stops
    instead of guessing what R would do."""


class NotReproducible(UserWarning):
    """The R package gives, for these inputs, a result decided by rounding error: another
    order of the samples, another machine or another numerical library changes it, in R as
    well. The port returns a result of the same kind and warns that it is not R's."""


def is_na(values: FloatArray) -> NDArray[np.bool_]:
    """``is.na`` of R: true for ``NA`` and for ``NaN``."""
    return np.isnan(values)


def is_na_real(values: FloatArray) -> NDArray[np.bool_]:
    """True where the value is ``NA_real_`` and not an ordinary ``NaN`` (``R_IsNA``)."""
    v = np.ascontiguousarray(values, dtype=np.float64)
    low = v.view(np.uint64) & np.uint64(0xFFFFFFFF)
    return np.asarray(np.isnan(v) & (low == np.uint64(1954)), dtype=bool)


@dataclass
class NamedVector:
    """A numeric vector with optional names."""

    values: FloatArray
    names: Names | None = None
    _index: dict[str, int] | None = field(default=None, repr=False, compare=False)

    def __post_init__(self) -> None:
        self.values = np.asarray(self.values, dtype=np.float64)
        if self.names is not None and len(self.names) != len(self.values):
            raise ValueError("names and values of different length")

    def __len__(self) -> int:
        return len(self.values)

    def copy(self) -> NamedVector:
        return NamedVector(self.values.copy(), None if self.names is None else list(self.names))

    def first_index(self) -> dict[str, int]:
        """Position of the first element carrying each name (missing names never match)."""
        if self._index is None:
            index: dict[str, int] = {}
            for i, name in enumerate(self.names or []):
                if name is not None and name not in index:
                    index[name] = i
            self._index = index
        return self._index


@dataclass
class RList:
    """A list of named vectors, itself with names: the ``correlationLists`` of the package."""

    items: list[NamedVector] = field(default_factory=list)
    names: Names | None = None

    def __len__(self) -> int:
        return len(self.items)

    def get(self, name: Name) -> NamedVector | None:
        """``x[[name]]``: the first element with that name, ``NULL`` when there is none."""
        if name is None or self.names is None:
            return None
        for item_name, item in zip(self.names, self.items, strict=True):
            if item_name == name:
                return item
        return None


@dataclass
class RMatrix:
    """A numeric matrix with optional row and column names."""

    values: FloatArray
    rownames: Names | None = None
    colnames: Names | None = None

    def __post_init__(self) -> None:
        self.values = np.asarray(self.values, dtype=np.float64)
        if self.values.ndim != 2:
            raise ValueError("a matrix needs two dimensions")

    @property
    def nrow(self) -> int:
        return int(self.values.shape[0])

    @property
    def ncol(self) -> int:
        return int(self.values.shape[1])


# ---- sets -------------------------------------------------------------------------------------
def is_in(x: Sequence[Name], table: Iterable[Name]) -> list[bool]:
    """``x %in% table``. A missing string matches a missing string, as in ``match``."""
    members = set(table)
    return [value in members for value in x]


def unique(x: Iterable[Name]) -> Names:
    """``unique``: first occurrences, in order."""
    return list(dict.fromkeys(x))


def union(x: Sequence[Name], y: Sequence[Name]) -> Names:
    """``union(x, y)``: ``unique(c(x, y))``."""
    return unique([*x, *y])


def intersect(x: Sequence[Name], y: Sequence[Name]) -> Names:
    """``intersect(x, y)``: the distinct values of ``x`` found in ``y``, in the order of ``x``."""
    members = set(y)
    return unique(value for value in x if value in members)


def setdiff(x: Sequence[Name], y: Sequence[Name]) -> Names:
    """``setdiff(x, y)``: the distinct values of ``x`` not found in ``y``."""
    members = set(y)
    return unique(value for value in x if value not in members)


# ---- the two helpers of the package --------------------------------------------------------------
def removeNode(to_remove: Sequence[Name], x: Names) -> Names:
    """``removeNode``: ``x`` without the elements found in ``to.remove``."""
    if len(to_remove) > 0:
        members = set(to_remove)
        if any(value in members for value in x):
            return [value for value in x if value not in members]
    return x


def remove_by_name_vector(x: Sequence[Name], vec: NamedVector) -> NamedVector:
    """``removeByName(x, l)`` for a vector: ``l[!(names(l) %in% x)]``.

    Without names, ``names(l) %in% x`` is ``logical(0)`` and the result is empty.
    """
    if vec.names is None:
        return NamedVector(np.empty(0), None)
    members = set(x)
    keep = [name not in members for name in vec.names]
    return NamedVector(
        vec.values[np.array(keep, dtype=bool)],
        [n for n, k in zip(vec.names, keep, strict=True) if k],
    )


def remove_by_name_list(x: Sequence[Name], lst: RList) -> RList:
    """``removeByName(x, l)`` for a list of neighbourhoods."""
    if lst.names is None:
        return RList([], None)
    members = set(x)
    keep = [name not in members for name in lst.names]
    return RList(
        [item for item, k in zip(lst.items, keep, strict=True) if k],
        [name for name, k in zip(lst.names, keep, strict=True) if k],
    )


def concat_lists(a: RList, b: RList) -> RList:
    """``c(a, b)`` for two lists. The result has names when either has them; an element of a
    list without names gets the empty name."""
    if (a.names is None and b.names is None) or len(a.items) + len(b.items) == 0:
        return RList([*a.items, *b.items], None)  # nothing to name: no names
    names_a: Names = list(a.names) if a.names is not None else [""] * len(a.items)
    names_b: Names = list(b.names) if b.names is not None else [""] * len(b.items)
    return RList([*a.items, *b.items], [*names_a, *names_b])


# ---- indexing ---------------------------------------------------------------------------------
def subset_logical(vec: NamedVector, mask: Sequence[bool | None]) -> NamedVector:
    """``vec[mask]`` for a logical ``mask`` of the same length.

    Where the mask is ``NA`` R does not skip the element: it returns a missing value, and a
    missing name with it. This is how an entry named ``NA`` enters a neighbourhood when a
    correlation is missing.
    """
    if len(mask) != len(vec):
        raise ValueError("mask and vector of different length")
    values: list[float] = []
    names: Names = []
    for i, flag in enumerate(mask):
        if flag is None:
            values.append(NA_REAL)
            names.append(None)
        elif flag:
            values.append(float(vec.values[i]))
            names.append(None if vec.names is None else vec.names[i])
    return NamedVector(np.array(values, dtype=np.float64), None if vec.names is None else names)


def compare(values: FloatArray, op: str, threshold: float) -> list[bool | None]:
    """``values < threshold`` or ``values > threshold``: ``NA`` where the value is missing."""
    with np.errstate(invalid="ignore"):
        if op == "<":
            result = values < threshold
        elif op == ">":
            result = values > threshold
        else:
            raise ValueError(op)
    missing = np.isnan(values) | np.isnan(threshold)
    return [None if m else bool(r) for r, m in zip(result, missing, strict=True)]


def get_by_names(vec: NamedVector, names: Sequence[Name]) -> NamedVector:
    """``vec[names]`` for a character index: a name that is not there gives a missing value
    under a missing name. A vector without names gives missing values without names."""
    if vec.names is None:
        return NamedVector(np.full(len(names), NA_REAL), None)
    index = vec.first_index()
    out_values = np.empty(len(names), dtype=np.float64)
    out_names: Names = []
    for k, name in enumerate(names):
        pos = index.get(name) if name is not None else None
        if pos is None:
            out_values[k] = NA_REAL
            out_names.append(None)
        else:
            out_values[k] = vec.values[pos]
            out_names.append(name)
    return NamedVector(out_values, out_names)


def set_by_names(vec: NamedVector, names: Sequence[Name], new_values: FloatArray) -> NamedVector:
    """``vec[names] <- new_values`` for a character index.

    A name that is not in the vector is not an error: R appends a new element with that name
    (``stringSubscript`` in ``subscript.c``). A missing name never matches, so every missing
    name appends an element of its own; a new name given twice refers to the same new element.
    """
    if len(names) != len(new_values):
        raise ValueError("names and values of different length")
    if vec.names is None:
        raise UnportedState("assignment by name into a vector without names")
    values = list(vec.values)
    out_names: Names = list(vec.names)
    index = dict(vec.first_index())
    appended: dict[str, int] = {}
    for name, value in zip(names, new_values, strict=True):
        pos = index.get(name) if name is not None else None
        if pos is None and name is not None:
            pos = appended.get(name)
        if pos is None:
            values.append(float(value))
            out_names.append(name)
            if name is not None:
                appended[name] = len(values) - 1
        else:
            values[pos] = float(value)
    return NamedVector(np.array(values, dtype=np.float64), out_names)


def set_names(vec: NamedVector | None, names: Sequence[Name]) -> NamedVector | None:
    """``names(vec) <- names``. ``vec`` is ``None`` for ``NULL``.

    Fewer names than elements: the rest are missing. More names than elements: an error.
    Names for ``NULL`` are an error even when there are none (``names`` is a character
    vector here, never ``NULL``).
    """
    if vec is None:
        raise RError("attempt to set an attribute on NULL")
    if len(names) > len(vec):
        raise RError(
            f"'names' attribute [{len(names)}] must be the same length as the vector [{len(vec)}]"
        )
    padded: Names = [*names, *([None] * (len(vec) - len(names)))]
    return NamedVector(vec.values, padded)
