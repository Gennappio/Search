"""Unit cases for the verbatim port (``search.r2014``): one function, one call, R's answer.

    r2014_units.py run <out dir> [--scale K] [--report report.csv] [--summary report.md]
    r2014_units.py check <out dir> [--report report.csv] [--summary report.md]
    r2014_units.py golden <out dir> <tests/data/r2014>

The whole-run comparison (``compare_r2014.py``) shows that ``SEARCH()`` and
``SEARCH.randomise()`` return what R returns. It does not say which function is right when the
whole is: two errors can cancel, and a helper that is never reached with a difficult input is
never tested. Here every function of the package is called on its own, in R and in the port,
with inputs chosen for it: the ordinary ones, and those at the edges (nothing, one element, no
names, missing names, names given twice, missing values, ties).

A case is ``{"id", "fn", "group", "args", "note"}``. ``fn`` is the name of a function of the
package or, in the group ``RULES``, the R expression a rule of the port stands for. ``args``
are tagged values (see ``num``, ``chr_``, ``mat``, ``lst``) that both sides read exactly:
numbers travel as hexadecimal floating point. ``run`` writes the cases, has
``r/r2014_units.R`` answer them, runs the port and compares; ``check`` compares the port again
with answers R gave before, without R; ``golden`` packs cases and answers for
``tests/unit/test_r2014_units.py``, which needs no R.

Everything is derived from fixed random seeds: a case keeps its number between runs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import struct
import subprocess
import sys
import tempfile
import warnings
from collections import Counter
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

import r2014_cases  # noqa: E402
import r2014_lib as lib  # noqa: E402

from search.r2014 import (  # noqa: E402
    association,
    core,
    glm,
    platform,
    randomise,
    rio,
    rsem,
    summarise,
)
from search.r2014.rsem import (  # noqa: E402
    NA_REAL,
    NamedVector,
    NotReproducible,
    RError,
    RList,
    RMatrix,
    UnportedState,
)

Json = dict[str, Any]
REPO = Path(__file__).resolve().parents[3]
PACKAGE = "package"
RULES = "rules of R"

REGRESSION = ("nodes.glm", "nodes.fastLm")
"""The two functions whose numbers come out of other linear-algebra routines than R's: their
coefficients are compared with a tolerance, everything else bit for bit."""


# ---- tagged values ------------------------------------------------------------------------------
def _is_na(x: float) -> bool:
    return math.isnan(x) and struct.unpack("<Q", struct.pack("<d", x))[0] & 0xFFFFFFFF == 1954


def enc_number(x: float) -> str:
    """A double as text both languages read back exactly: hexadecimal floating point."""
    x = float(x)
    if math.isnan(x):
        return "NA" if _is_na(x) else "NaN"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    return x.hex()


def dec_number(text: str) -> float:
    if text == "NA":
        return NA_REAL
    return float.fromhex(text) if "0x" in text else float(text)


def num(values: Sequence[float], names: Sequence[str | None] | None = None) -> Json:
    return {
        "t": "num",
        "v": [enc_number(v) for v in values],
        "names": None if names is None else list(names),
    }


def chr_(values: Sequence[str | None]) -> Json:
    return {"t": "chr", "v": list(values), "names": None}


def lgl(value: bool) -> Json:
    return {"t": "lgl", "v": [bool(value)]}


def lst(items: Sequence[Json], names: Sequence[str | None] | None = None) -> Json:
    return {"t": "list", "items": list(items), "names": None if names is None else list(names)}


def mat(
    values: Any,
    rownames: Sequence[str | None] | None = None,
    colnames: Sequence[str | None] | None = None,
) -> Json:
    array = np.asarray(values, dtype=np.float64)
    return {
        "t": "mat",
        "dim": [int(array.shape[0]), int(array.shape[1])],
        "v": [enc_number(v) for v in array.T.ravel()],  # by columns, as R stores a matrix
        "rownames": None if rownames is None else list(rownames),
        "colnames": None if colnames is None else list(colnames),
    }


NULL: Json = {"t": "null"}


def scalar(value: Any) -> Json:
    if value is None:
        return NULL
    if isinstance(value, bool):
        return lgl(value)
    if isinstance(value, str):
        return chr_([value])
    return num([float(value)])


def port_value(v: Json) -> Any:
    """A tagged value as the object the port works with."""
    kind = v["t"]
    if kind == "null":
        return None
    if kind == "num":
        return NamedVector(np.array([dec_number(s) for s in v["v"]], dtype=np.float64), v["names"])
    if kind == "chr":
        return list(v["v"])
    if kind == "lgl":
        return [bool(x) for x in v["v"]]
    if kind == "mat":
        nrow, ncol = v["dim"]
        values = np.array([dec_number(s) for s in v["v"]], dtype=np.float64)
        return RMatrix(values.reshape(ncol, nrow).T.copy(), v["rownames"], v["colnames"])
    if kind == "list":
        return RList([port_value(item) for item in v["items"]], v["names"])
    if kind == "df":
        columns: dict[str, list[Any]] = {}
        for col in v["columns"]:
            if col["type"] == "numeric":
                columns[col["name"]] = [dec_number(s) for s in col["v"]]
            else:
                columns[col["name"]] = list(col["v"])
        return glm.Covariates(columns)
    raise ValueError(kind)


def plain(v: Json) -> Any:
    """A tagged value of one element as a Python scalar."""
    if v["t"] == "null":
        return None
    if v["t"] == "num":
        return dec_number(v["v"][0])
    return v["v"][0]


_WHOLE = {"max_iter", "n_boot", "n_boot_start", "n_boot_NULL", "batches_NULL", "iteration"}


def _keywords(args: Json, skip: Sequence[str] = ()) -> Json:
    """The scalar arguments of a case with the port's names: dots become underscores."""
    out: Json = {}
    for name, value in args.items():
        if name in skip:
            continue
        key = name.replace(".", "_")
        out[key] = plain(value)
        if key in _WHOLE and isinstance(out[key], float) and out[key] == int(out[key]):
            out[key] = int(out[key])
    return out


# ---- what the port returns, described as R's answer is ------------------------------------------
def describe(value: Any) -> Json:
    if value is None:
        return {"t": "null"}
    if isinstance(value, NamedVector):
        return {"t": "num", "v": value.values, "names": value.names}
    if isinstance(value, RMatrix):
        return {
            "t": "mat",
            "dim": [value.nrow, value.ncol],
            "v": value.values.T.ravel(),
            "rownames": value.rownames,
            "colnames": value.colnames,
        }
    if isinstance(value, RList):
        return {
            "t": "list",
            "items": [describe(item) for item in value.items],
            "names": value.names,
        }
    if isinstance(value, core.SearchResult):
        fields: Json = {
            "removed": value.removed,
            "correlationLists": value.correlationLists,
            "correlationMatrix": value.correlationMatrix,
            "D": value.D,
            "O": value.O,
            "network.growth": value.network_growth,
            "other.data": value.other_data,
        }
        if value.nodes_degree or getattr(value, "_has_degree", False):
            fields["nodes.degree"] = value.nodes_degree
        return describe(fields)
    if isinstance(value, summarise.IterationSummary):
        return describe({"real": value.real, "NULL": value.null})
    if isinstance(value, association.AssociationData):
        return describe({"CL": value.CL, "other.data": value.other_data})
    if isinstance(value, dict):
        return {
            "t": "list",
            "items": [describe(item) for item in value.values()],
            "names": list(value) if value else None,
        }
    if isinstance(value, str):
        return {"t": "chr", "v": [value], "names": None}
    if isinstance(value, (bool, int, float)):
        return {"t": "num", "v": np.array([float(value)]), "names": None}
    if isinstance(value, (list, tuple)):
        if len(value) == 0:
            return {"t": "empty"}
        if all(item is None or isinstance(item, str) for item in value):
            return {"t": "chr", "v": list(value), "names": None}
        if all(isinstance(item, (bool, int, float)) for item in value):
            return {"t": "num", "v": np.array([float(item) for item in value]), "names": None}
        return {"t": "list", "items": [describe(item) for item in value], "names": None}
    raise TypeError(type(value))


def _numbers(v: Any) -> np.ndarray:
    if isinstance(v, np.ndarray):
        return v
    return np.array([dec_number(s) for s in v], dtype=np.float64)


def _zero_length(v: Json) -> bool:
    if v["t"] == "empty":
        return True
    if v["t"] in ("num", "chr"):
        return len(v["v"]) == 0
    return v["t"] == "list" and len(v["items"]) == 0


def differences(r: Json, py: Json, what: str = "value", tolerance: float = 0.0) -> list[str]:
    """Where a value R returned and the one the port returned are not the same. Names must be
    the same strings in the same order, missing where R has them missing; numbers the same
    double, bit for bit, with NA and NaN kept apart (within ``tolerance``, relative, when it is
    not zero).

    One thing is not looked at: the type of a vector without elements. R tells ``logical(0)``
    from ``character(0)`` and from ``list()``; the port keeps an empty Python list for all
    three. Whether it has names is compared.
    """
    if _zero_length(r) and _zero_length(py):
        r_names, py_names = r.get("names"), py.get("names")
        if (r_names is None) != (py_names is None):
            return [f"{what}: empty, names R {r_names!r} / port {py_names!r}"]
        return []
    if r["t"] != py["t"]:
        return [f"{what}: R has {_kind(r)}, the port has {_kind(py)}"]
    kind = r["t"]
    out: list[str] = []
    if kind == "null":
        return out
    if kind in ("num", "chr") and r["names"] != py["names"]:
        out.append(f"{what}: names R {lib._short(r['names'])} / port {lib._short(py['names'])}")
    if kind == "chr":
        if list(r["v"]) != list(py["v"]):
            out.append(f"{what}: R {lib._short(r['v'])} / port {lib._short(py['v'])}")
        return out
    if kind == "num":
        out += _number_differences(what, _numbers(r["v"]), _numbers(py["v"]), tolerance)
        return out
    if kind == "mat":
        if list(r["dim"]) != list(py["dim"]):
            return [f"{what}: dim R {r['dim']} / port {py['dim']}"]
        for key in ("rownames", "colnames"):
            if r[key] != py[key]:
                out.append(f"{what}: {key} R {lib._short(r[key])} / port {lib._short(py[key])}")
        if tolerance > 0 and what.endswith("other.data") and not out:
            return _coefficient_differences(what, r, py, tolerance)
        out += _number_differences(what, _numbers(r["v"]), _numbers(py["v"]), tolerance)
        return out
    if kind == "list":
        if len(r["items"]) != len(py["items"]):
            return [f"{what}: {len(r['items'])} elements in R, {len(py['items'])} in the port"]
        if r["names"] != py["names"]:
            return [f"{what}: names R {lib._short(r['names'])} / port {lib._short(py['names'])}"]
        labels = r["names"] or [str(i + 1) for i in range(len(r["items"]))]
        for label, r_item, py_item in zip(labels, r["items"], py["items"], strict=True):
            out += differences(r_item, py_item, f"{what}${label}", tolerance)
        return out
    return [f"{what}: R returned {_kind(r)}, which the comparison does not know"]


def _kind(v: Json) -> str:
    if v["t"] in ("num", "chr"):
        return f"{v['t']}[{len(v['v'])}]"
    if v["t"] == "list":
        return f"list[{len(v['items'])}]"
    return str(v["t"])


def _number_differences(what: str, r: np.ndarray, py: np.ndarray, tolerance: float) -> list[str]:
    if r.shape != py.shape:
        return [f"{what}: {r.size} numbers in R, {py.size} in the port"]
    n_bad, largest = lib.same_numbers(r, py, tolerance)
    if n_bad == 0:
        return []
    return [f"{what}: {n_bad} of {r.size} numbers differ, largest difference {largest:.3g}"]


def _coefficient_differences(what: str, r: Json, py: Json, tolerance: float) -> list[str]:
    """The table of all coefficients of a regression function, as ``compare_r2014`` compares
    it: where a seed is regressed on itself only the estimates are looked at."""
    nrow, ncol = py["dim"]
    d = lib.Differences()
    d.coefficients(
        what,
        _numbers(r["v"]).reshape(ncol, nrow).T.copy(),
        _numbers(py["v"]).reshape(ncol, nrow).T.copy(),
        py["rownames"] or [],
        py["colnames"] or [],
        tolerance,
    )
    return d.items


# ---- the port, one function at a time -----------------------------------------------------------
def _results_object(v: Json) -> core.SearchResult:
    fields = dict(zip(v["names"], v["items"], strict=True))

    def per_iteration(name: str) -> dict[str, Any]:
        item = fields.get(name)
        if item is None or item["t"] != "list" or not item["items"]:
            return {}
        return {str(k): port_value(x) for k, x in zip(item["names"], item["items"], strict=True)}

    result = core.SearchResult(
        removed=port_value(fields["removed"]),
        correlationLists=port_value(fields["correlationLists"]),
        correlationMatrix=per_iteration("correlationMatrix"),
        D=per_iteration("D"),
        O=per_iteration("O"),
        network_growth=per_iteration("network.growth"),
        other_data=per_iteration("other.data"),
        nodes_degree=per_iteration("nodes.degree"),
    )
    result._has_degree = "nodes.degree" in fields  # type: ignore[attr-defined]
    return result


def _association(function: Callable[..., Any]) -> Callable[[Json, str], Any]:
    def call(a: Json, work: str) -> Any:
        covariates = port_value(a["covariates"]) if "covariates" in a else None
        return function(
            seeds=port_value(a["seeds"]),
            x=port_value(a["x"]),
            Blist=port_value(a["Blist"]),
            covariates=covariates,
            **_keywords(a, skip=("seeds", "x", "Blist", "covariates")),
        )

    return call


def _search(a: Json, work: str) -> Any:
    return core.SEARCH(
        x=port_value(a["x"]),
        Alist=port_value(a["Alist"]),
        Blist=port_value(a["Blist"]),
        **_keywords(a, skip=("x", "Alist", "Blist")),
    )


def _write_data_files(a: Json, work: str) -> Any:
    return randomise.write_data_files(
        _results_object(a["results.obj"]), f"{work}/", int(plain(a["iteration"]))
    )


def _randomise_inner(a: Json, work: str) -> Any:
    keywords = _keywords(a, skip=("x", "Alist", "Blist", "iteration", "sampling", "sampling.seed"))
    return randomise.SEARCH_randomise_inner(
        int(plain(a["iteration"])),
        port_value(a["x"]),
        Alist=port_value(a["Alist"]),
        Blist=port_value(a["Blist"]),
        covariates=None,
        sampling=plain(a["sampling"]),
        sampling_seed=plain(a["sampling.seed"]),
        output_directory=work,
        sample_kind=platform.current().sample_kind,
        **keywords,
    )


def _randomise(a: Json, work: str) -> Any:
    keywords = _keywords(a, skip=("data", "seed.nodes", "universe.nodes", "directory"))
    directory = plain(a["directory"]) if "directory" in a else None
    return randomise.SEARCH_randomise(
        data=port_value(a["data"]),
        seed_nodes=port_value(a["seed.nodes"]),
        universe_nodes=port_value(a["universe.nodes"]),
        output_directory=work if directory is None else f"{work}/{directory}",
        **keywords,
    )


def _write_gz(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as handle:
        handle.write(text.encode("utf-8"))


def _summarise(a: Json, work: str) -> Any:
    files = a["files"]
    for name, item in zip(files["names"] or [], files["items"], strict=True):
        _write_gz(Path(work) / name, item["v"][0])
    return summarise.summarise_iterations(
        output_directory=f"{work}/real/",
        output_directory_NULL=f"{work}/null/",
        **_keywords(a, skip=("files",)),
    )


def _remove_by_name(a: Json, work: str) -> Any:
    x, target = port_value(a["x"]) or [], port_value(a["l"])
    if isinstance(target, RList):
        return rsem.remove_by_name_list(x, target)
    return rsem.remove_by_name_vector(x, target)


def _matrix_format(a: Json, work: str) -> Any:
    o_matrix, correlation = core.get_matrix_format(
        port_value(a["x"]), port_value(a["all.features"])
    )
    return {"O": o_matrix, "subnetworks.correlation": correlation}


def _nodes_degree(a: Json, work: str) -> Any:
    v = a["correlation.matrix.list"]
    return core.estimate_nodes_degree(
        {str(k): port_value(x) for k, x in zip(v["names"] or [], v["items"], strict=True)}
    )


def _write_table(a: Json, work: str) -> Any:
    x = port_value(a["x"])
    if isinstance(x, NamedVector):
        return rio.table_of_vector(x)
    if isinstance(x, RMatrix):
        return rio.table_of_matrix(x)
    return rio.table_of_strings(x)


def _read_table(a: Json, work: str) -> Any:
    path = Path(work) / "table.txt.gz"
    _write_gz(path, plain(a["text"]))
    return rio.read_connectivity(str(path))


def _values(v: Json) -> np.ndarray:
    return np.asarray(port_value(v).values)


def _threshold(a: Json, absolute: bool, op: str) -> Any:
    v = port_value(a["v"])
    values = np.abs(v.values) if absolute else v.values
    return rsem.subset_logical(v, rsem.compare(values, op, plain(a["t"])))


def _colon(a: Json, work: str) -> Any:
    end = plain(a["b"])
    return summarise._colon(int(plain(a["a"])), None if end is None else int(end))


PORT: dict[str, Callable[[Json, str], Any]] = {
    # the functions of the package
    "getO": lambda a, w: core.getO(port_value(a["x"]), port_value(a["y"])),
    "getOmax": lambda a, w: core.getOmax(port_value(a["x"]), port_value(a["y"])),
    "removeNode": lambda a, w: rsem.removeNode(
        port_value(a["to.remove"]) or [], port_value(a["x"])
    ),
    "removeByName": _remove_by_name,
    "get.matrix.format": _matrix_format,
    "estimate.nodes.degree": _nodes_degree,
    "nodes.correlation": _association(association.nodes_correlation),
    "nodes.glm": _association(glm.nodes_glm),
    "nodes.fastLm": _association(glm.nodes_fastLm),
    "SEARCH": _search,
    "write.data.files": _write_data_files,
    ".SEARCH.randomise.inner": _randomise_inner,
    "SEARCH.randomise": _randomise,
    "summarise.iterations": _summarise,
    # the rules of R
    "v": lambda a, w: port_value(a["v"]),
    "v[names] <- values": lambda a, w: rsem.set_by_names(
        port_value(a["v"]), port_value(a["names"]), _values(a["values"])
    ),
    "v[names]": lambda a, w: rsem.get_by_names(port_value(a["v"]), port_value(a["names"])),
    "v[v < t]": lambda a, w: _threshold(a, False, "<"),
    "v[v > t]": lambda a, w: _threshold(a, False, ">"),
    "v[abs(v) > t]": lambda a, w: _threshold(a, True, ">"),
    "names(v) <- names": lambda a, w: rsem.set_names(
        port_value(a["v"]), port_value(a["names"]) or []
    ),
    "union(a, b)": lambda a, w: rsem.union(port_value(a["a"]), port_value(a["b"])),
    "intersect(a, b)": lambda a, w: rsem.intersect(port_value(a["a"]), port_value(a["b"])),
    "setdiff(a, b)": lambda a, w: rsem.setdiff(port_value(a["a"]), port_value(a["b"])),
    "unique(a)": lambda a, w: rsem.unique(port_value(a["a"])),
    "a %in% b": lambda a, w: rsem.is_in(port_value(a["a"]), port_value(a["b"])),
    "c(a, b)": lambda a, w: rsem.concat_lists(port_value(a["a"]), port_value(a["b"])),
    "sort(v, decreasing = TRUE)": lambda a, w: randomise._sort_decreasing(port_value(a["v"])),
    'as.numeric(sprintf("%.3f", v))': lambda a, w: NamedVector(
        randomise._three_decimals(_values(a["v"])), None
    ),
    "a:b": _colon,
    "write.table(x)": _write_table,
    "read.table(file)": _read_table,
}
WRITES_FILES = ("write.data.files", ".SEARCH.randomise.inner", "SEARCH.randomise")

NOT_PORTED = (
    "make.default.plots",
    "make.networkplot",
    "make.scatter.plot",
    "get.adjacency.matrices",
    "get.program.defaults",
)
"""The functions of the package the port does not have: the figures and their two helpers."""


def run_port(case: Json) -> Json:
    """The port's answer to a case: ``value`` as ``describe`` gives it, or ``error`` with R's
    message, or ``unported`` when the port refuses the input; ``files`` when the function
    writes any."""
    out: Json = {"error": None, "value": None}
    with tempfile.TemporaryDirectory() as work, warnings.catch_warnings():
        warnings.simplefilter("ignore", NotReproducible)
        try:
            out["value"] = describe(PORT[case["fn"]](case["args"], work))
        except RError as exc:
            out["error"] = exc.message
        except UnportedState as exc:
            out["unported"] = str(exc)
        if case["fn"] in WRITES_FILES:
            out["files"] = randomise.read_output_files(work)
    return out


def compare(case: Json, r: Json, port: Json, tolerance: float = 0.0) -> list[str]:
    """Everything in which the port's answer to a case is not R's."""
    if "unported" in port:
        return [f"the port refuses this input: {port['unported']}"]
    if r.get("error") is not None or port["error"] is not None:
        out = (
            []
            if r.get("error") == port["error"]
            else [f"error: R {r.get('error')!r} / port {port['error']!r}"]
        )
    else:
        out = differences(r["value"], port["value"], "value", tolerance)
    if "files" in r or "files" in port:
        r_files, py_files = r.get("files") or {}, port.get("files") or {}
        if sorted(r_files) != sorted(py_files):
            out.append(f"files written: R {sorted(r_files)} / port {sorted(py_files)}")
        for name in sorted(set(r_files) & set(py_files)):
            if r_files[name] != py_files[name]:
                out.append(f"file {name}: {lib._first_line_diff(r_files[name], py_files[name])}")
    return out


def tolerance_of(case: Json, regression: float) -> float:
    return regression if case["fn"] in REGRESSION else 0.0


# ---- cases --------------------------------------------------------------------------------------
POOL = [f"g{i}" for i in range(1, 9)]
Rng = np.random.Generator


class _Unit:
    """The cases of one function, numbered in the order they are added."""

    def __init__(self, fn: str, slug: str, group: str = PACKAGE) -> None:
        self.fn, self.slug, self.group = fn, slug, group
        self.cases: list[Json] = []

    def add(self, args: Json, note: str = "", **extra: Any) -> None:
        self.cases.append(
            {
                "id": f"{self.slug}/{len(self.cases) + 1:03d}",
                "fn": self.fn,
                "group": self.group,
                "note": note,
                "args": args,
                **extra,
            }
        )


def _number(rng: Rng) -> float:
    kind = int(rng.integers(5))
    if kind == 0:
        return float(rng.integers(-3, 4))
    if kind == 1:
        return float(np.round(rng.normal(), 2))
    return float(rng.normal())


def _some_numbers(rng: Rng, n: int, p_missing: float = 0.0) -> list[float]:
    out = [_number(rng) for _ in range(n)]
    for i in range(n):
        if rng.random() < p_missing:
            out[i] = NA_REAL if rng.random() < 0.7 else float("nan")
    return out


def _some_names(
    rng: Rng, n: int, pool: Sequence[str] = POOL, p_dup: float = 0.0, p_na: float = 0.0
) -> list[str | None]:
    names: list[str | None]
    if n > len(pool) or rng.random() < p_dup:
        names = [str(pool[int(i)]) for i in rng.integers(0, len(pool), size=n)]
    else:
        names = [str(pool[int(i)]) for i in rng.permutation(len(pool))[:n]]
    for i in range(n):
        if rng.random() < p_na:
            names[i] = None
    return names


def _neighbourhood(
    rng: Rng,
    pool: Sequence[str] = POOL,
    p_nonames: float = 0.0,
    p_dup: float = 0.0,
    p_na: float = 0.0,
    p_missing: float = 0.0,
) -> Json:
    k = int(rng.integers(0, len(pool) + 1))
    values = _some_numbers(rng, k, p_missing)
    if rng.random() < p_nonames:
        return num(values)
    return num(values, _some_names(rng, k, pool, p_dup, p_na))


def _neighbourhoods(rng: Rng, n: int, prefix: str = "s", **kw: Any) -> Json:
    return lst(
        [_neighbourhood(rng, **kw) for _ in range(n)], [f"{prefix}{i + 1}" for i in range(n)]
    )


def _ones(names: Sequence[str | None] | None, n: int = 0) -> Json:
    return num([1.0] * (n if names is None else len(names)), names)


def cases_getO(n: int) -> list[Json]:
    u = _Unit("getO", "getO")
    rng = np.random.default_rng(1001)
    a, b = ["g1", "g2", "g3"], ["g2", "g3", "g4"]
    u.add({"x": _ones(a), "y": _ones(a)}, "the same members: 1")
    u.add({"x": _ones(a), "y": _ones(["g4", "g5"])}, "no member in common: 0")
    u.add({"x": _ones(a), "y": _ones(b)}, "two of four")
    u.add({"x": _ones(a), "y": _ones(["g2"])}, "one inside the other")
    u.add({"x": _ones(["g1", "g1", "g2"]), "y": _ones(["g1"])}, "a member given twice")
    u.add({"x": _ones([]), "y": _ones([])}, "two empty neighbourhoods: 0 / 0")
    u.add({"x": _ones(None, 2), "y": _ones(None, 3)}, "no names on either side: 0 / 0")
    u.add({"x": _ones(None, 0), "y": _ones(None, 0)}, "empty and without names")
    u.add({"x": _ones(None, 2), "y": _ones(a)}, "names on one side only")
    u.add({"x": _ones(a), "y": _ones([])}, "one empty")
    u.add({"x": _ones(["g1", None]), "y": _ones([None, "g2"])}, "a missing name on both sides")
    u.add({"x": _ones(["g1", None]), "y": _ones(["g1"])}, "a missing name on one side")
    u.add({"x": _ones([None, None]), "y": _ones([None])}, "only missing names")
    for _ in range(n):
        kw = {"p_nonames": 0.1, "p_dup": 0.25, "p_na": 0.1}
        u.add({"x": _neighbourhood(rng, **kw), "y": _neighbourhood(rng, **kw)}, "random")
    return u.cases


def cases_getOmax(n: int) -> list[Json]:
    u = _Unit("getOmax", "getOmax")
    rng = np.random.default_rng(1002)
    same = lst([_ones(["g1", "g2"])] * 3, ["s1", "s2", "s3"])
    u.add({"x": same, "y": same}, "three identical lists: the last one gets 0", defect="R2")
    two = lst([_ones(["g1", "g2"]), _ones(["g1", "g2"])], ["s1", "s2"])
    u.add({"x": two, "y": two}, "two identical lists", defect="R2")
    u.add({"x": lst([]), "y": same}, "no list in x: NULL")
    u.add({"x": lst([]), "y": lst([])}, "no list at all")
    u.add({"x": same, "y": lst([])}, "no list in y")
    one = lst([_ones(["g1", "g2", "g3"])], ["s1"])
    u.add({"x": one, "y": one}, "one list compared with itself: never compared")
    u.add({"x": one, "y": lst([*one["items"], _ones(["g2", "g3"])], ["s1", "s2"])}, "one new list")
    u.add({"x": same, "y": one}, "x longer than y")
    bare = lst([_ones(None, 2), _ones(None, 1), _ones(None, 3)])
    u.add({"x": bare, "y": bare}, "lists without names: 0 / 0 for every pair")
    mixed = lst([_ones(None, 2), _ones(["g1"]), _ones([]), _ones(["g1", "g2"])])
    u.add({"x": mixed, "y": mixed}, "lists with and without names")
    for _ in range(n):
        kw = {"p_nonames": 0.08, "p_dup": 0.1, "p_na": 0.05}
        x = _neighbourhoods(rng, int(rng.integers(0, 6)), **kw)
        style = rng.random()
        if style < 0.3:
            y, note = x, "random, y is x (keep.Alist = FALSE)"
        elif style < 0.7:
            old = _neighbourhoods(rng, int(rng.integers(0, 5)), prefix="o", **kw)
            y = lst([*x["items"], *old["items"]], [*x["names"], *old["names"]])
            note = "random, y is c(x, old lists) (keep.Alist = TRUE)"
        else:
            y, note = _neighbourhoods(rng, int(rng.integers(0, 6)), **kw), "random, unrelated y"
        u.add({"x": x, "y": y}, note)
    return u.cases


def cases_removeNode(n: int) -> list[Json]:
    u = _Unit("removeNode", "removeNode")
    rng = np.random.default_rng(1003)
    x = chr_(["g1", "g2", "g3", "g4"])
    u.add({"to.remove": chr_(["g2"]), "x": x}, "one element")
    u.add({"to.remove": chr_([]), "x": x}, "nothing to remove")
    u.add({"to.remove": NULL, "x": x}, "NULL to remove")
    u.add({"to.remove": chr_(["g9"]), "x": x}, "an element that is not there")
    u.add({"to.remove": chr_(["g4", "g1", "g3", "g2"]), "x": x}, "everything")
    u.add({"to.remove": chr_(["g1"]), "x": chr_([])}, "nothing to remove from")
    u.add({"to.remove": chr_(["g1"]), "x": chr_(["g1", "g2", "g1"])}, "an element twice in x")
    u.add({"to.remove": chr_([None]), "x": chr_(["g1", None, "g2"])}, "a missing value")
    u.add({"to.remove": chr_(["g1", "g1"]), "x": x}, "an element twice in to.remove")
    for _ in range(n):
        u.add(
            {
                "to.remove": chr_(_some_names(rng, int(rng.integers(0, 5)), p_dup=0.2, p_na=0.1)),
                "x": chr_(_some_names(rng, int(rng.integers(0, 8)), p_dup=0.15, p_na=0.08)),
            },
            "random",
        )
    return u.cases


def cases_removeByName(n: int) -> list[Json]:
    u = _Unit("removeByName", "removeByName")
    rng = np.random.default_rng(1004)
    vector = num([0.5, 0.25, 0.75], ["g1", "g2", "g3"])
    lists = lst([_ones(["g1"]), _ones(["g2", "g3"]), _ones([])], ["g1", "g2", "g3"])
    for target, kind in ((vector, "vector"), (lists, "list")):
        u.add({"x": chr_(["g2"]), "l": target}, f"{kind}: one name")
        u.add({"x": chr_([]), "l": target}, f"{kind}: no name")
        u.add({"x": chr_(["g9"]), "l": target}, f"{kind}: a name that is not there")
        u.add({"x": chr_(["g1", "g2", "g3"]), "l": target}, f"{kind}: every name")
    u.add({"x": chr_(["g1"]), "l": num([1.0, 2.0])}, "a vector without names: nothing is kept")
    u.add({"x": chr_([]), "l": num([1.0, 2.0])}, "a vector without names, nothing to remove")
    u.add({"x": chr_(["g1"]), "l": lst([_ones(["g1"])])}, "a list without names")
    u.add({"x": chr_(["g1"]), "l": num([1.0, 2.0, 3.0], ["g1", "g2", "g1"])}, "a name twice")
    u.add({"x": chr_([None]), "l": num([1.0, 2.0], ["g1", None])}, "a missing name")
    u.add({"x": chr_(["g1"]), "l": num([1.0, 2.0], ["g1", None])}, "a missing name that stays")
    for _ in range(n):
        x = chr_(_some_names(rng, int(rng.integers(0, 5)), p_dup=0.2, p_na=0.1))
        if rng.random() < 0.5:
            k = int(rng.integers(0, 8))
            names = None if rng.random() < 0.08 else _some_names(rng, k, p_dup=0.15, p_na=0.1)
            u.add({"x": x, "l": num(_some_numbers(rng, k, 0.1), names)}, "random vector")
        else:
            k = int(rng.integers(0, 6))
            items = [_neighbourhood(rng) for _ in range(k)]
            names = None if rng.random() < 0.08 else _some_names(rng, k, p_dup=0.15, p_na=0.1)
            u.add({"x": x, "l": lst(items, names)}, "random list")
    return u.cases


def cases_get_matrix_format(n: int) -> list[Json]:
    u = _Unit("get.matrix.format", "get_matrix_format")
    rng = np.random.default_rng(1005)
    features = chr_(POOL)

    def call(x: Json, note: str, all_features: Json = features, **extra: Any) -> None:
        u.add({"x": x, "all.features": all_features}, note, **extra)

    full = lst(
        [
            num([1.0, 0.6, -0.7], ["g1", "g2", "g5"]),
            num([0.6, 1.0], ["g1", "g2"]),
            num([1.0, 0.9, 0.8, -0.55], ["g3", "g4", "g5", "g6"]),
        ],
        ["g1", "g2", "g3"],
    )
    call(full, "three lists")
    call(lst([num([1.0, 0.6], ["g1", "g2"])], ["g1"]), "one list")
    call(lst([num([1.0], ["g1"])], ["g1"]), "one list of one member: the names are dropped")
    call(
        lst([num([1.0], ["g1"]), num([0.7], ["g1"])], ["g1", "g2"]),
        "two lists with the same single member: a column without names",
    )
    call(lst([num([], []), num([1.0, 0.5], ["g2", "g3"])], ["g1", "g2"]), "an empty list")
    call(lst([num([]), num([1.0, 0.5], ["g2", "g3"])], ["g1", "g2"]), "empty, without names")
    call(lst([num([], []), num([], [])], ["g1", "g2"]), "only empty lists: no column")
    call(lst([num([1.0, 0.5], ["g1", "g9"])], ["g1"]), "a member that is not a feature")
    call(lst([num([1.0, 0.5], ["g1", None])], ["g1"]), "a member with a missing name")
    call(
        lst([num([1.0, NA_REAL], ["g1", "g2"]), num([NA_REAL, 1.0], ["g2", "g3"])], ["g1", "g3"]),
        "a column left with missing values only is dropped",
    )
    call(
        lst([num([1.0, float("nan")], ["g1", "g2"]), num([0.5], ["g1"])], ["g1", "g3"]),
        "NaN counts as missing",
    )
    call(lst([num([1.0, 0.5, 0.25], ["g1", "g2", "g1"])], ["g1"]), "a member twice in a list")
    call(lst([num([1.0, 0.5], ["g2", "g3"])], ["zz"]), "a list whose name is not a feature")
    call(full, "features in another order", chr_(list(reversed(POOL))))
    call(lst([num([1.0, 0.5])], ["g1"]), "a list of values without names")
    call(
        lst([num([1.0], ["g1"]), num([1.0], ["g2"])]),
        "lists without names",
        outside="lists without names: R stops at the first step, the port does not model it",
    )
    call(
        lst([num([1.0], ["g1"]), num([1.0], ["g2"])], ["g1", "g1"]),
        "two lists with the same name",
        outside="two lists with the same name",
    )
    call(
        lst([num([1.0], ["g1"]), num([1.0], ["g2"])], ["g1", None]),
        "a list with a missing name",
        outside="a list with a missing name",
    )
    call(
        full,
        "a feature given twice",
        chr_([*POOL, "g1"]),
        outside="duplicated row names in the data matrix",
    )
    for _ in range(n):
        k = int(rng.integers(1, 6))
        seeds = _some_names(rng, k)
        order = [POOL[int(i)] for i in rng.permutation(len(POOL))]
        kw = {"pool": POOL if rng.random() < 0.9 else [*POOL, "g9"], "p_missing": 0.12}
        items = [_neighbourhood(rng, p_dup=0.05, **kw) for _ in range(k)]
        call(lst(items, seeds), "random", chr_(order))
    return u.cases


def _random_matrix_with_missing(rng: Rng, p: float) -> np.ndarray:
    k, m = int(rng.integers(1, 6)), int(rng.integers(0, 7))
    values = rng.normal(size=(k, m))
    missing = rng.random((k, m)) < p
    values[missing] = NA_REAL
    values[(rng.random((k, m)) < 0.04) & ~missing] = float("nan")
    return values


def cases_estimate_nodes_degree(n: int) -> list[Json]:
    u = _Unit("estimate.nodes.degree", "estimate_nodes_degree")
    rng = np.random.default_rng(1006)
    na = NA_REAL

    def call(matrices: Sequence[Json], note: str) -> None:
        names = [str(i + 1) for i in range(len(matrices))]
        u.add({"correlation.matrix.list": lst(matrices, names)}, note)

    seeds, members = ["g1", "g2", "g3"], ["g1", "g2", "g3", "g4"]
    full = mat([[1, 0.6, na, 0.7], [0.6, 1, na, na], [na, na, 1, 0.8]], seeds, members)
    call([full], "three seeds, one of them selected by no other")
    call([mat([[1, 0.5], [na, 0.6]], ["g1", "g2"], ["g1", "g5"])], "one seed among the members")
    call([mat([[0.5, 0.6]], ["g1"], ["g4", "g5"])], "no seed among the members")
    call([mat([[1.0], [0.7]], None, None)], "a column without any name")
    call([mat(np.empty((2, 0)), ["g1", "g2"], None)], "no column")
    call([mat([[1.0]], ["g1"], ["g1"])], "one seed, itself")
    call([mat([[na, na], [na, na]], ["g1", "g2"], ["g1", "g2"])], "nothing but missing values")
    call([mat([[float("nan"), na], [na, 1]], ["g1", "g2"], ["g1", "g2"])], "NaN is counted")
    call([full, mat([[1, na], [0.6, 1]], ["g1", "g2"], ["g2", "g1"])], "two iterations")
    call([], "no iteration")
    for _ in range(n):
        matrices = []
        for _ in range(int(rng.integers(1, 4))):
            values = _random_matrix_with_missing(rng, 0.4)
            rows = _some_names(rng, values.shape[0])
            cols = _some_names(rng, values.shape[1])
            matrices.append(mat(values, rows, cols if cols else None))
        call(matrices, "random")
    return u.cases


def _expression(rng: Rng, **kw: Any) -> tuple[np.ndarray, list[str], list[str]]:
    """A small expression matrix as the whole-run cases make them; a missing value is NA."""
    x, rownames, colnames = r2014_cases._random_matrix(rng, **kw)
    x = np.where(np.isnan(x), NA_REAL, x)
    return x, [str(r) for r in rownames], [str(c) for c in colnames]


def _random_universe(rng: Rng, rownames: Sequence[str]) -> list[str]:
    style = rng.random()
    if style < 0.4:
        return list(rownames)
    order = [str(rownames[int(i)]) for i in rng.permutation(len(rownames))]
    if style < 0.7:
        return order
    if style < 0.9:
        return order[: max(2, len(order) - int(rng.integers(1, 4)))]
    return [*order, "not_in_matrix"]


METHODS = ("pearson", "spearman", "kendall")
USES = ("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")


def cases_nodes_correlation(n: int) -> list[Json]:
    u = _Unit("nodes.correlation", "nodes_correlation")
    rng = np.random.default_rng(1007)
    genes, samples = ["g1", "g2", "g3", "g4", "g5"], [f"s{i + 1}" for i in range(6)]
    values = np.array(
        [
            [5.0, 5.0, 5.0, 5.0, 5.0, 5.01],
            [9.3, 8.95, 7.3, 6.73, 7.38, 8.04],
            [5.67, 7.78, 6.75, 7.27, 7.46, 7.68],
            [8.41, 9.04, 7.87, 9.37, 7.33, 8.35],
            [8.9, 8.09, 7.26, 7.08, 7.54, 8.22],
        ]
    )
    x = mat(values, genes, samples)

    def call(note: str, matrix: Json = x, **kw: Any) -> None:
        extra = {k: kw.pop(k) for k in ("defect",) if k in kw}
        args: Json = {
            "seeds": chr_(kw.pop("seeds", ["g3"])),
            "x": matrix,
            "Blist": chr_(kw.pop("Blist", genes)),
        }
        args.update({name.replace("_", "."): scalar(value) for name, value in kw.items()})
        u.add(args, note, **extra)

    call("the defaults: Spearman above 0.5")
    for method in METHODS:
        call(f"{method}, threshold 0", correlation_method=method, correlation_threshold=0.0)
    call("two seeds", seeds=["g3", "g5"])
    call("the universe in another order", Blist=["g5", "g1", "g4", "g3", "g2"])
    call("a universe smaller than the matrix", Blist=["g3", "g4", "g5"])
    call("a universe with a gene the matrix does not have", Blist=[*genes, "zz"])
    call("a threshold nothing passes", correlation_threshold=1.0)
    call("a negative threshold: everything passes", correlation_threshold=-0.1)
    call("a seed outside the universe", seeds=["g3"], Blist=["g1", "g2", "g4"])
    call("a seed the matrix does not have", seeds=["zz"])
    call("a universe of one gene", Blist=["g3"])
    call("no gene of the universe in the matrix", Blist=["zz"])
    call("a correlation type that is not supported", correlation_type="partial")
    call("a method that does not exist", correlation_method="cosine")
    call("an abbreviated method", correlation_method="pear")
    call("a rule for missing values that does not exist", correlation_use="sometimes")
    holed = values.copy()
    holed[1, 2] = NA_REAL
    for use in USES:
        for method in METHODS:
            call(
                f"a missing value, {method}, use = {use}",
                mat(holed, genes, samples),
                correlation_method=method,
                correlation_use=use,
                correlation_threshold=0.3,
                defect="R7" if use == "everything" else None,
            )
    flat = values.copy()
    flat[0] = 5.0
    call("a gene with the same value in every sample", mat(flat, genes, samples))
    call("a seed with the same value in every sample", mat(flat, genes, samples), seeds=["g1"])
    gaps = values.copy()
    gaps[:, 0] = NA_REAL
    gaps[2, 1:] = [1.0, NA_REAL, NA_REAL, NA_REAL, 2.0]
    call(
        "pairs of two samples only",
        mat(gaps, genes, samples),
        correlation_use="pairwise.complete.obs",
        correlation_method="pearson",
    )
    empty = np.full_like(values, NA_REAL)
    for use in ("complete.obs", "na.or.complete", "pairwise.complete.obs"):
        call(
            f"nothing but missing values, use = {use}",
            mat(empty, genes, samples),
            correlation_use=use,
        )
    for k in range(n):
        m, rownames, colnames = _expression(rng, p_missing=0.3)
        if rng.random() < 0.1:
            m[int(rng.integers(m.shape[0])), int(rng.integers(m.shape[1]))] = float("nan")
        universe = _random_universe(rng, rownames)
        pick = rng.permutation(len(rownames))[: int(rng.integers(1, 4))]
        call(
            "random",
            mat(m, rownames, colnames),
            seeds=[rownames[int(i)] for i in pick],
            Blist=universe,
            correlation_threshold=float(rng.choice([0.0, 0.3, 0.5, 0.5, 0.8, 1.0, -0.1])),
            correlation_method=METHODS[k % 3],
            correlation_use=USES[(k // 3) % 5],
        )
    return u.cases


def _covariates(rng: Rng, samples: Sequence[str]) -> Json:
    n = len(samples)
    return {
        "t": "df",
        "rownames": list(samples),
        "columns": [
            {
                "name": "age",
                "type": "numeric",
                "v": [enc_number(float(np.round(rng.normal(60, 10), 1))) for _ in range(n)],
            },
            {
                "name": "batch",
                "type": "character",
                "v": [str(rng.choice(["a", "b", "c"])) for _ in range(n)],
            },
        ],
    }


FORMULAS = ("age", "batch", "age + batch", "x.node:batch", "age*batch", "batch + x.node:age")
GLM_FAMILIES = (
    "gaussian",
    "gaussian",
    "Gamma",
    "inverse.gaussian",
    "poisson",
    "quasi",
    "quasipoisson",
)
ADJUSTMENTS = ("BH", "bonferroni", "holm", "BY", "hochberg", "hommel", "fdr", "none")


def _regression_cases(function: str, slug: str, n: int, seed: int) -> list[Json]:
    u = _Unit(function, slug)
    rng = np.random.default_rng(seed)
    fast = function == "nodes.fastLm"

    def matrix() -> tuple[np.ndarray, list[str], list[str]]:
        # enough samples and no small whole numbers, so that no model fits its seed exactly;
        # no gene with the same value in every sample (R defects R11 and R12)
        while True:
            m, rownames, colnames = _expression(
                rng,
                min_samples=12,
                max_samples=18,
                styles=("continuous", "blocks"),
                p_missing=0.0,
                p_constant=0.0,
            )
            if all(np.ptp(row) > 0 for row in m):
                return m, rownames, colnames

    def call(note: str, m: np.ndarray, rownames: list[str], colnames: list[str], **kw: Any) -> None:
        extra = {k: kw.pop(k) for k in ("defect",) if k in kw}
        args: Json = {
            "seeds": chr_(kw.pop("seeds")),
            "x": mat(m, rownames, colnames),
            "Blist": chr_(kw.pop("Blist", rownames)),
        }
        if "covariates" in kw:
            args["covariates"] = kw.pop("covariates")
        args.update({name.replace("_", "."): scalar(value) for name, value in kw.items()})
        u.add(args, note, **extra)

    m, rownames, colnames = matrix()
    call("the defaults", m, rownames, colnames, seeds=rownames[:2])
    call(
        "a threshold of 1: every gene",
        m,
        rownames,
        colnames,
        seeds=rownames[:1],
        significance_threshold=1.0,
    )
    call(
        "a threshold nothing passes: the seed vanishes from the result",
        m, rownames, colnames, seeds=rownames[:2], significance_threshold=0.0, defect="R6",
    )  # fmt: skip
    call(
        "only the seed passes: its neighbourhood loses the name",
        m, rownames, colnames, seeds=rownames[:2], significance_threshold=1e-12, defect="R5",
    )  # fmt: skip
    call(
        "a family that is not supported",
        m,
        rownames,
        colnames,
        seeds=rownames[:1],
        glm_family="binomial",
    )
    call("a universe of one gene", m, rownames, colnames, seeds=rownames[:1], Blist=rownames[:1])
    call("a seed the matrix does not have", m, rownames, colnames, seeds=["zz"])
    covariates = _covariates(rng, colnames)
    call(
        "covariates without a formula are ignored",
        m, rownames, colnames, seeds=rownames[:1], covariates=covariates,
    )  # fmt: skip
    for formula in FORMULAS:
        call(
            f"formula.string = {formula}",
            m, rownames, colnames, seeds=rownames[:1], covariates=covariates,
            formula_string=formula, significance_threshold=0.2,
        )  # fmt: skip
    for adjustment in ADJUSTMENTS:
        call(
            f"adjust.p = {adjustment}",
            m, rownames, colnames, seeds=rownames[:2], adjust_p=adjustment,
            significance_threshold=0.2,
        )  # fmt: skip
    if fast:
        flat = m.copy()
        flat[len(flat) - 1] = 8.0
        call(
            "a gene with the same value in every sample gets a coefficient",
            flat, rownames, colnames, seeds=rownames[:1], significance_threshold=1.0, defect="R12",
        )  # fmt: skip
    else:
        for family in sorted(set(GLM_FAMILIES)):
            call(
                f"glm.family = {family}",
                m,
                rownames,
                colnames,
                seeds=rownames[:2],
                glm_family=family,
            )
        call(
            "negative values under Gamma",
            m - 9.0,
            rownames,
            colnames,
            seeds=rownames[:1],
            glm_family="Gamma",
        )
    # a batch seen in one sample only, and a gene without a value in that sample: its model has
    # one coefficient fewer, and the tables of the genes cannot be bound
    holed = m.copy()
    holed[3, -1] = NA_REAL
    lone = {
        "t": "df",
        "rownames": list(colnames),
        "columns": [
            {
                "name": "batch",
                "type": "character",
                "v": [*(["a", "b"] * len(colnames))[: len(colnames) - 1], "c"],
            }
        ],
    }
    call(
        "a model with one coefficient fewer for the fourth gene",
        holed,
        rownames,
        colnames,
        seeds=rownames[:1],
        covariates=lone,
        formula_string="batch",
    )
    for _ in range(n):
        m, rownames, colnames = matrix()
        if rng.random() < 0.1:
            m[int(rng.integers(m.shape[0])), int(rng.integers(m.shape[1]))] = NA_REAL
        pick = rng.permutation(len(rownames))[: int(rng.integers(1, 4))]
        kw: Json = {
            "seeds": [rownames[int(i)] for i in pick],
            "Blist": [rownames[int(i)] for i in rng.permutation(len(rownames))],
            "significance_threshold": float(rng.choice([0.01, 0.05, 0.05, 0.2, 0.5, 1e-5])),
        }
        if not fast:
            kw["glm_family"] = str(rng.choice(GLM_FAMILIES))
        if rng.random() < 0.3:
            kw["adjust_p"] = str(rng.choice(ADJUSTMENTS))
        if rng.random() < 0.4:
            kw["covariates"] = _covariates(rng, colnames)
            kw["formula_string"] = str(rng.choice(FORMULAS))
        call("random", m, rownames, colnames, **kw)
    return u.cases


def _search_arguments(rng: Rng, **fixed: Any) -> Json:
    args: Json = {
        "scale": bool(rng.random() < 0.85),
        "D.lower": float(rng.choice([0.0, 0.0, 0.1, 0.25, 0.5])),
        "D.higher": float(rng.choice([0.0, 0.1, 0.24, 0.5])),
        "O.lower": float(rng.choice([0.0, 0.01, 0.1, 0.3, 0.6])),
        "keep.Alist": bool(rng.random() < 0.6),
        "max.iter": int(rng.choice([1, 2, 3, 4])),
        "sd.threshold": float(rng.choice([0.0, 0.0, 0.1, 0.5])),
        "verbose": False,
        "correlation.threshold": float(rng.choice([0.3, 0.5, 0.5, 0.7])),
        "correlation.method": str(rng.choice(METHODS)),
    }
    args.update(fixed)
    return args


def _tagged(args: Json) -> Json:
    return {name: scalar(value) for name, value in args.items()}


def _search_case(rng: Rng, p_missing: float = 0.05, **fixed: Any) -> Json:
    m, rownames, colnames = _expression(rng, p_missing=p_missing)
    pick = rng.permutation(len(rownames))[: int(rng.integers(1, 5))]
    return {
        "x": mat(m, rownames, colnames),
        "Alist": chr_([rownames[int(i)] for i in pick]),
        "Blist": chr_(_random_universe(rng, rownames)),
        **_tagged(_search_arguments(rng, **fixed)),
    }


def _port_search(args: Json) -> core.SearchResult | None:
    try:
        result = _search(args, "")
    except (RError, UnportedState):
        return None
    return result if isinstance(result, core.SearchResult) else None


def _shows_r3(result: core.SearchResult) -> bool:
    # the seed removed on Omax is counted again from its list: D gets an element for it with
    # a missing value, and the seeds added include a missing name
    return any(None in growth for growth in result.network_growth.values())


def _shows_r4(result: core.SearchResult) -> bool:
    removed = {name for name in result.removed if name is not None}
    return any(removed & set(item.names or []) for item in result.correlationLists.items)


def cases_SEARCH(n: int) -> list[Json]:
    u = _Unit("SEARCH", "SEARCH")
    rng = np.random.default_rng(1010)
    genes, samples = ["g1", "g2", "g3", "g4", "g5"], [f"s{i + 1}" for i in range(6)]
    values = np.array(
        [
            [5.0, 5.0, 5.0, 5.0, 5.0, 5.01],
            [9.3, 8.95, 7.3, 6.73, 7.38, 8.04],
            [5.67, 7.78, 6.75, 7.27, 7.46, 7.68],
            [8.41, 9.04, 7.87, 9.37, 7.33, 8.35],
            [8.9, 8.09, 7.26, 7.08, 7.54, 8.22],
        ]
    )
    x = mat(values, genes, samples)

    def call(note: str, matrix: Json = x, **kw: Any) -> None:
        extra = {k: kw.pop(k) for k in ("defect",) if k in kw}
        args: Json = {
            "x": matrix,
            "Alist": chr_(kw.pop("Alist", ["g3", "g4", "g5"])),
            "Blist": chr_(kw.pop("Blist", genes)),
            "verbose": lgl(False),
        }
        args.update({name.replace("_", "."): scalar(value) for name, value in kw.items()})
        u.add(args, note, **extra)

    call("the defaults")
    call(
        "the sd prefilter removes the gene that carries the name, not the flat one",
        sd_threshold=0.1,
        defect="R1",
    )
    call("one seed: D is divided by zero", Alist=["g3"])
    call("no seed", Alist=[])
    call("max.iter = 0: the loop is not entered", max_iter=0)
    call("three iterations", max_iter=3, D_higher=0.3)
    call("without scaling", scale=False, max_iter=2)
    call("an association function that does not exist", nodes_association_function="nodes.nothing")
    call("a seed the matrix does not have", Alist=["g3", "zz"])
    call("seeds that are not in the universe", Blist=["g1", "g2"])
    call(
        "a universe of one gene and the sd prefilter", Alist=["g3"], Blist=["g3"], sd_threshold=0.1
    )
    holed = values.copy()
    holed[1, 2] = NA_REAL
    call(
        "a missing value: an entry named NA, and the run stops when the matrices are built",
        mat(holed, genes, samples), Alist=["g3", "g4"], D_higher=-1.0, max_iter=2, defect="R7",
    )  # fmt: skip
    # cases found by trying random inputs on the port until the symptom shows; R confirms it
    for defect, shows, fixed in (
        ("R3", _shows_r3, {"keep.Alist": False}),
        ("R4", _shows_r4, {"max.iter": 3}),
    ):
        found = 0
        for _ in range(20000):
            args = _search_case(rng, p_missing=0.0, **fixed)
            result = _port_search(args)
            if result is not None and shows(result):
                u.add(args, f"found by search: shows {defect}", defect=defect)
                found += 1
                if found == 3:
                    break
    for _ in range(n):
        u.add(_search_case(rng), "random")
    return u.cases


def _fraction(rng: Rng) -> float:
    kind = rng.random()
    if kind < 0.6:
        d = int(rng.integers(1, 12))
        return int(rng.integers(0, d + 1)) / d
    if kind < 0.75:
        return float(rng.choice([0.123455, 0.123445, 0.000005, 0.5, 1.0, 0.0, 2.675e-3]))
    if kind < 0.85:
        return -1.0 / int(rng.integers(1, 9))
    return float(rng.random())


def _results_object_case(rng: Rng) -> Json:
    def per_iteration(make: Callable[[], Json], k: int) -> Json:
        if k == 0:
            return lst([])  # list(): no names
        return lst([make() for _ in range(k)], [str(i + 1) for i in range(k)])

    def correlation_matrix() -> Json:
        nrow = int(rng.integers(1, 5))
        ncol = int(rng.choice([0, 1, 1, 2, 3, 4, 5, 6]))
        values = rng.normal(size=(nrow, ncol))
        values[rng.random((nrow, ncol)) < 0.35] = NA_REAL
        if ncol == 1 and rng.random() < 0.6:
            return mat(values, None, None)  # what get.matrix.format leaves of a single column
        return mat(values, _some_names(rng, nrow), _some_names(rng, ncol) if ncol else None)

    def overlap_matrix() -> Json:
        k = int(rng.integers(1, 5))
        values = np.triu(np.array([[_fraction(rng) for _ in range(k)] for _ in range(k)]))
        names = _some_names(rng, k)
        return mat(np.abs(values), names, names)

    def connectivity() -> Json:
        k = int(rng.integers(1, 9))
        values = [_fraction(rng) for _ in range(k)]
        for i in range(k):
            roll = rng.random()
            if roll < 0.05:
                values[i] = float("nan")
            elif roll < 0.08:
                values[i] = float("inf") if rng.random() < 0.5 else float("-inf")
            elif roll < 0.11:
                values[i] = NA_REAL
        return num(values, _some_names(rng, k))

    def growth() -> Json:
        return chr_([str(v) for v in _some_names(rng, int(rng.integers(0, 4)))])

    def degree() -> Json:
        k = int(rng.integers(0, 6))
        if k == 1 and rng.random() < 0.7:
            return num([_fraction(rng)])  # one common node: the names are dropped
        return num([abs(_fraction(rng)) for _ in range(k)], _some_names(rng, k))

    def coefficients() -> Json:
        return mat(
            rng.normal(size=(2, 4)) * 10.0 ** rng.integers(-8, 8, size=(2, 4)),
            ["g1.(Intercept)", "g1.x.node"],
            ["g1.Estimate", "g1.P", "g2.Estimate", "g2.P"],
        )

    n_matrix = int(rng.choice([0, 1, 1, 2]))
    n_d = n_matrix + int(rng.integers(0, 2))
    other = per_iteration(coefficients, n_d) if rng.random() < 0.2 else lst([])
    return lst(
        [
            chr_([str(v) for v in _some_names(rng, int(rng.integers(0, 3)))]),
            _neighbourhoods(rng, int(rng.integers(0, 4))),
            per_iteration(correlation_matrix, n_matrix),
            per_iteration(connectivity, n_d),
            per_iteration(overlap_matrix, n_matrix),
            per_iteration(growth, n_d),
            other,
            per_iteration(degree, n_matrix),
        ],
        [
            "removed",
            "correlationLists",
            "correlationMatrix",
            "D",
            "O",
            "network.growth",
            "other.data",
            "nodes.degree",
        ],
    )


def _object(**fields: Json) -> Json:
    base: Json = {
        "removed": chr_([]),
        "correlationLists": lst([]),
        "correlationMatrix": lst([]),
        "D": lst([]),
        "O": lst([]),
        "network.growth": lst([]),
        "other.data": lst([]),
        "nodes.degree": lst([]),
    }
    base.update({name.replace("_", "."): value for name, value in fields.items()})
    return lst(list(base.values()), list(base))


def cases_write_data_files(n: int) -> list[Json]:
    u = _Unit("write.data.files", "write_data_files")
    rng = np.random.default_rng(1011)
    na = NA_REAL

    def call(obj: Json, note: str, iteration: int = 1) -> None:
        u.add({"results.obj": obj, "iteration": num([float(iteration)])}, note)

    one = ["1"]
    call(_object(), "nothing to write")
    d = num([0.123455, 0.5, 0.123445, 1.0, 0.5, 0.0], ["g1", "g2", "g3", "g4", "g5", "g6"])
    call(_object(D=lst([d], one)), "D is rounded to five decimals and sorted, ties in order")
    call(_object(D=lst([d, d], ["1", "2"])), "two iterations of D", iteration=7)
    odd = num([float("nan"), 0.5, na, float("inf"), -0.25], ["g1", "g2", "g3", "g4", "g5"])
    call(_object(D=lst([odd], one)), "D with NaN, NA and Inf: missing values are dropped")
    call(_object(D=lst([num([0.5, 0.25], ["g1", None])], one)), "D with a missing name")
    call(_object(D=lst([num([0.5, 0.25], ["g1", "g1"])], one)), "D with a name twice")
    call(_object(D=lst([num([0.5, 0.25])], one)), "D without names")
    cm = mat([[1, 0.61234, na], [0.61234, 1, -0.5555]], ["g1", "g2"], ["g1", "g2", "g3"])
    call(_object(correlationMatrix=lst([cm], one)), "a matrix goes through three decimals")
    row = mat([[1, 0.61234, na]], ["g1"], ["g1", "g2", "g3"])
    call(_object(correlationMatrix=lst([row], one)), "a matrix of one row")
    call(
        _object(correlationMatrix=lst([mat([[1.0], [0.7]], None, None)], one)),
        "a column without names is written as a row",
    )
    call(
        _object(correlationMatrix=lst([mat([[1.0], [0.7]], ["g1", "g2"], ["g1"])], one)),
        "a named column of two rows",
    )
    call(_object(correlationMatrix=lst([mat([[1.0]], ["g1"], ["g1"])], one)), "one by one")
    call(
        _object(correlationMatrix=lst([mat(np.empty((2, 0)), ["g1", "g2"], None)], one)),
        "a matrix without columns",
    )
    call(
        _object(correlationMatrix=lst([mat(np.empty((1, 0)), ["g1"], None)], one)),
        "one row and no column",
    )
    o = mat([[1, 1 / 3, 0], [0, 1, 2 / 3], [0, 0, 1]], ["g1", "g2", "g3"], ["g1", "g2", "g3"])
    call(_object(O=lst([o], one)), "the O matrix")
    call(_object(O=lst([mat([[float("nan")]], ["g1"], ["g1"])], one)), "an O of NaN")
    call(
        _object(network_growth=lst([chr_(["g4", "g5"]), chr_([])], ["1", "2"])),
        "seeds added, and an iteration that adds none",
    )
    degree = num([2 / 3, 1 / 3, 2 / 3], ["g1", "g2", "g3"])
    call(_object(nodes_degree=lst([degree], one)), "node degrees")
    call(_object(nodes_degree=lst([num([1.0])], one)), "one node degree without a name")
    call(_object(nodes_degree=lst([num([])], one)), "no node degree")
    table = mat(
        [[1.23456789012345678e-7, 0.01, 123456.789, 1e-300], [1e15, 1e16, 0.1 + 0.2, 100000.0]],
        ["g1.(Intercept)", "g1.x.node"],
        ["g1.Estimate", "g1.P", "g2.Estimate", "g2.P"],
    )
    call(_object(other_data=lst([table], one)), "coefficients are written with 15 digits")
    for _ in range(n):
        call(_results_object_case(rng), "random", iteration=int(rng.choice([1, 2, 12])))
    return u.cases


def _randomise_arguments(rng: Rng) -> Json:
    args = _search_arguments(rng)
    args["max.iter"] = int(rng.choice([1, 2, 3]))
    args["sampling"] = str(rng.choice(["bootstrap", "bootstrap", "jacknife"]))
    args["sampling.seed"] = float(rng.integers(1, 1000))
    return args


def cases_randomise_inner(n: int) -> list[Json]:
    u = _Unit(".SEARCH.randomise.inner", "SEARCH_randomise_inner")
    rng = np.random.default_rng(1012)

    def call(note: str, **fixed: Any) -> None:
        m, rownames, colnames = _expression(rng, p_missing=0.03)
        pick = rng.permutation(len(rownames))[: int(rng.integers(1, 5))]
        args = _randomise_arguments(rng)
        args["keep.Alist"] = args.pop("keep.Alist")
        args.update(fixed)
        u.add(
            {
                "iteration": num([float(args.pop("iteration", int(rng.integers(1, 6))))]),
                "x": mat(m, rownames, colnames),
                "Alist": chr_([rownames[int(i)] for i in pick]),
                "Blist": chr_(_random_universe(rng, rownames)),
                **_tagged(args),
            },
            note,
        )

    call("bootstrap", sampling="bootstrap")
    call("jacknife: one sample fewer, none twice", sampling="jacknife")
    call("a sampling that does not exist", sampling="permutation")
    call("the seed given as text", **{"sampling.seed": "17"})
    call("replicate 1 with seed 1", iteration=1, **{"sampling.seed": 1.0})
    call("replicate 0", iteration=0, **{"sampling.seed": 5.0})
    for _ in range(n):
        call("random")
    return u.cases


def cases_SEARCH_randomise(n: int) -> list[Json]:
    u = _Unit("SEARCH.randomise", "SEARCH_randomise")
    rng = np.random.default_rng(1013)

    def call(note: str, matrix: Json | None = None, **fixed: Any) -> None:
        extra = {k: fixed.pop(k) for k in ("defect",) if k in fixed}
        m, rownames, colnames = _expression(rng, p_missing=0.03)
        pick = rng.permutation(len(rownames))[: int(rng.integers(1, 5))]
        args = _randomise_arguments(rng)
        args["keep.seeds"] = args.pop("keep.Alist")
        args["n.boot"] = int(rng.choice([1, 2, 3]))
        args["n.boot.start"] = int(rng.choice([1, 1, 1, 4]))
        seeds = fixed.pop("seeds", [rownames[int(i)] for i in pick])
        universe = fixed.pop("universe", _random_universe(rng, rownames))
        args.update(fixed)
        u.add(
            {
                "data": mat(m, rownames, colnames) if matrix is None else matrix,
                "seed.nodes": chr_(seeds),
                "universe.nodes": chr_(universe),
                **_tagged(args),
            },
            note,
            **extra,
        )

    genes, samples = ["g1", "g2", "g3", "g4"], ["s1", "s2", "s3", "s4", "s5"]
    small = np.array([[1.0, 2, 3, 4, 5], [2.0, 1, 4, 3, 6], [5.0, 4, 3, 2, 1], [1.0, 3, 2, 5, 4]])
    x = mat(small, genes, samples)
    call("three replicates", x, seeds=["g1", "g2"], universe=genes, **{"n.boot": 3})
    call("a matrix of one row", mat(small[:1], genes[:1], samples), seeds=["g1"], universe=genes)
    call(
        "a matrix of one column",
        mat(small[:, :1], genes, samples[:1]),
        seeds=["g1"],
        universe=genes,
    )
    call("no seed in the matrix", x, seeds=["zz"], universe=genes)
    call(
        "seeds given twice, and one the matrix does not have",
        x,
        seeds=["g2", "g1", "g2", "zz"],
        universe=genes,
    )
    call(
        "n.boot = 0: replicates 1 and 0",
        x,
        seeds=["g1", "g2"],
        universe=genes,
        **{"n.boot": 0, "n.boot.start": 1},
    )
    call(
        "replicates 4 to 6",
        x,
        seeds=["g1", "g2"],
        universe=genes,
        **{"n.boot": 3, "n.boot.start": 4},
    )
    call("jacknife", x, seeds=["g1", "g2"], universe=genes, sampling="jacknife")
    call("a sampling that does not exist", x, seeds=["g1", "g2"], universe=genes, sampling="none")
    call("the seed given as text", x, seeds=["g1", "g2"], universe=genes, **{"sampling.seed": "3"})
    call(
        "a blank in the output directory: the files stay uncompressed",
        x, seeds=["g1", "g2"], universe=genes, directory="with blank", defect="R10",
    )  # fmt: skip
    call("the defaults of the package", x, seeds=["g1", "g2"], universe=genes, **{
        "D.lower": 0.1, "D.higher": 0.24, "O.lower": 0.01, "sd.threshold": 0.1, "max.iter": 1,
    })  # fmt: skip
    for _ in range(n):
        call("random")
    return u.cases


CONNECTIVITY_VALUES = (
    "0", "0.1", "0.1", "0.2", "0.25", "0.33333", "0.5", "0.5", "0.66667", "1", "-0.2", "NA", "Inf",
)  # fmt: skip


def _connectivity_text(rows: Sequence[tuple[str, str]]) -> str:
    return "\tx\n" + "".join(f"{name}\t{value}\n" for name, value in rows)


def _real(replicate: int, iteration: int = 1) -> str:
    return f"real/output/connectivity__iteration_{replicate}_{iteration}.txt.gz"


def _null(batch: int, replicate: int, iteration: int = 1) -> str:
    return (
        f"null/batch_{batch}_0_FALSE_0_0/output/"
        f"connectivity__iteration_{replicate}_{iteration}.txt.gz"
    )


def cases_summarise_iterations(n: int) -> list[Json]:
    u = _Unit("summarise.iterations", "summarise_iterations")
    rng = np.random.default_rng(1014)

    def call(files: dict[str, str], note: str, **kw: Any) -> None:
        extra = {k: kw.pop(k) for k in ("defect",) if k in kw}
        args: Json = {
            "files": lst([chr_([text]) for text in files.values()], list(files)),
            "max.iter": num([float(kw.pop("max_iter", 1))]),
            "n.boot": num([float(kw.pop("n_boot", 1))]),
        }
        for name in ("n_boot_NULL", "batches_NULL"):
            if name in kw:
                args[name.replace("_", ".")] = scalar(kw.pop(name))
        u.add(args, note, **extra)

    def table(**values: str) -> str:
        return _connectivity_text(list(values.items()))

    a = table(g1="0.5", g2="0.2", g3="0.1")
    b = table(g1="0.9", g2="0.4", g3="0.1")
    low = table(g1="0.1", g2="0.2", g3="0.1")
    nulls = {_null(1, 1): low, _null(2, 1): table(g1="0.7", g2="0.1", g3="0.2")}
    both = {"n_boot_NULL": 1, "batches_NULL": 2}
    call({_real(1): a, _real(2): b, **nulls}, "two replicates, two null runs", n_boot=2, **both)
    call(
        {_real(1): a, _real(3): b, **nulls},
        "replicate 2 has no file: replicate 3 is read and never used",
        n_boot=3, defect="R9", **both,
    )  # fmt: skip
    call(
        {_real(1): a, _null(1, 1): low, _null(1, 3): b},
        "null replicate 2 has no file: the same for the null runs",
        n_boot_NULL=3, batches_NULL=1, defect="R9",
    )  # fmt: skip
    call(
        {_real(1): a, _null(1, 1): table(g1="0.1", g9="0.9"), _null(2, 1): low},
        "a gene of a null run that no real run has",
        defect="R8", **both,
    )  # fmt: skip
    call({_real(1): a}, "no null run asked for: batches.NULL = NULL")
    call({_real(1): a}, "null runs asked for, none there", **both)
    call({_real(1): a, _null(1, 1): low}, "a single null file", n_boot_NULL=1, batches_NULL=1)
    call({_real(1): table(g1="0.5"), **nulls}, "a single gene", **both)
    call({**nulls}, "no real file at all", **both)
    call({_real(1): a, **nulls}, "an iteration no replicate reached", max_iter=2, **both)
    call(
        {_real(1): a, _real(1, 2): table(g4="0.3", g1="0.2"), **nulls,
         _null(1, 1, 2): table(g1="0.3", g4="0.3"), _null(2, 1, 2): table(g1="0.1", g4="0.5")},
        "two iterations: the genes of the first stay in the table of the second",
        max_iter=2, **both,
    )  # fmt: skip
    tenth = table(g1="0.1", g2="0.1")
    call(
        {**{_real(i): tenth for i in range(1, 11)}, **{_null(1, i): tenth for i in range(1, 4)}},
        "the mean of ten values of 0.1 against a null value of 0.1",
        n_boot=10, n_boot_NULL=3, batches_NULL=1,
    )  # fmt: skip
    call(
        {_real(1): a, _real(2): table(g3="0.3", g4="0.8"), **nulls},
        "replicates with different genes",
        n_boot=2, **both,
    )  # fmt: skip
    call(
        {_real(1): table(g1="0.5", g2="0.5", g3="0.5"), _real(2): table(g3="0.5", g1="0.5"),
         **nulls},
        "equal means keep the order of the genes",
        n_boot=2, **both,
    )  # fmt: skip
    call({_real(1): table(g1="NA", g2="0.2", g3="Inf"), **nulls}, "NA and Inf in a file", **both)
    numbered = {"007": "0.5", "12": "0.2", "100000": "0.1"}
    numbered_null = {
        _null(1, 1): _connectivity_text([("7", "0.1"), ("12", "0.3"), ("100000", "0.1")]),
        _null(2, 1): _connectivity_text([("7", "0.6"), ("12", "0.1"), ("100000", "0.2")]),
    }
    call(
        {_real(1): _connectivity_text(list(numbered.items())), **numbered_null},
        "gene names that look like whole numbers: 007 comes back as 7",
        defect="R13", **both,
    )  # fmt: skip
    decimal = [("1.50", "0.5"), ("100000", "0.2"), ("2e3", "0.1")]
    call(
        {_real(1): _connectivity_text(decimal), _real(2): _connectivity_text(decimal)},
        "gene names that look like decimal numbers",
        n_boot=2, defect="R13",
    )  # fmt: skip
    call(
        {_real(1): _connectivity_text([("1.0", "0.5"), ("1", "0.2")])},
        "two names that become the same number",
        defect="R13",
    )
    call(
        {_real(1): _connectivity_text([("NA", "0.5"), ("g1", "0.2")])},
        "a gene called NA",
        defect="R13",
    )
    call(
        {_real(1): _connectivity_text([("T", "0.5"), ("FALSE", "0.2")])},
        "gene names that are all logical words",
        defect="R13",
    )
    call(
        {_real(1): _connectivity_text([("T", "0.5"), ("g1", "0.2")])},
        "one logical word among names",
    )
    call({_real(1): "\tx\n", **nulls}, "a file with no gene", **both)
    call({_real(1): "", **nulls}, "an empty file", **both)
    genes = ["g1", "g2", "g3", "g4", "g5", "g6"]

    def random_table(pool: Sequence[str]) -> tuple[str, list[str]]:
        k = int(rng.integers(2, len(pool) + 1))
        chosen = [str(pool[int(i)]) for i in rng.permutation(len(pool))[:k]]
        rows = [(g, str(rng.choice(CONNECTIVITY_VALUES))) for g in chosen]
        return _connectivity_text(rows), chosen

    for _ in range(n):
        max_iter = int(rng.choice([1, 1, 2]))
        n_boot = int(rng.integers(1, 5))
        batches, n_null = int(rng.integers(1, 4)), int(rng.integers(1, 4))
        files: dict[str, str] = {}
        p_real, p_null = float(rng.choice([1.0, 1.0, 0.8])), float(rng.choice([1.0, 0.85]))
        seen: list[str] = []  # the genes of the real runs so far: the rows of the summary
        for iteration in range(1, max_iter + 1):
            for i in range(1, n_boot + 1):
                if rng.random() < p_real:
                    files[_real(i, iteration)], chosen = random_table(genes)
                    seen += [g for g in chosen if g not in seen]
            # mostly null runs whose genes the real runs have; otherwise R defect R8 stops it
            null_pool = seen if len(seen) >= 2 and rng.random() < 0.85 else [*genes, "g7"]
            for batch in range(1, batches + 1):
                for i in range(1, n_null + 1):
                    if rng.random() < p_null:
                        files[_null(batch, i, iteration)], _ = random_table(null_pool)
        call(
            files, "random",
            max_iter=max_iter, n_boot=n_boot, n_boot_NULL=n_null, batches_NULL=batches,
        )  # fmt: skip
    return u.cases


# ---- the rules of R -----------------------------------------------------------------------------
def cases_rules(n: int) -> list[Json]:
    rng = np.random.default_rng(1020)
    out: list[Json] = []
    na, nan, inf = NA_REAL, float("nan"), float("inf")

    def unit(fn: str, slug: str) -> _Unit:
        return _Unit(fn, slug, RULES)

    def named(k: int, p_dup: float = 0.15, p_na: float = 0.1, p_missing: float = 0.1) -> Json:
        return num(_some_numbers(rng, k, p_missing), _some_names(rng, k, p_dup=p_dup, p_na=p_na))

    def maybe_named(k: int) -> Json:
        return num(_some_numbers(rng, k, 0.1)) if rng.random() < 0.15 else named(k)

    def names(k: int) -> Json:
        return chr_(_some_names(rng, k, p_dup=0.25, p_na=0.12))

    u = unit("v", "numbers_travel_exactly")
    hard = [0.1, 1 / 3, 9.118641810491681, 5e-324, 1.7976931348623157e308, -0.0, 0.0, 1e22, 1e23]
    u.add({"v": num([*hard, na, nan, inf, -inf])}, "numbers hard to carry as text")
    u.add(
        {"v": num(list(rng.normal(size=300) * 10.0 ** rng.integers(-30, 30, size=300)))}, "random"
    )
    u.add({"v": num([0.5, na], ["g1", None])}, "a missing name")
    out += u.cases

    u = unit("v[names] <- values", "assign_by_name")
    v = num([1.0, 2.0, 3.0], ["g1", "g2", "g3"])
    u.add({"v": v, "names": chr_(["g2"]), "values": num([9.0])}, "a name that is there")
    u.add({"v": v, "names": chr_(["zz"]), "values": num([9.0])}, "a new name is appended")
    u.add(
        {
            "v": v,
            "names": chr_([None, None, "q", "q", "g1"]),
            "values": num([10.0, 11, 12, 13, 14]),
        },
        "every missing name is an element of its own; a new name twice is one element",
    )
    u.add({"v": v, "names": chr_([]), "values": num([])}, "no name")
    u.add(
        {
            "v": num([1.0, 2.0, 3.0], ["g1", "g1", None]),
            "names": chr_(["g1", None]),
            "values": num([8.0, 9.0]),
        },
        "a name the vector has twice, and one it has missing",
    )
    u.add(
        {"v": num([1.0, 2.0]), "names": chr_(["g1"]), "values": num([9.0])},
        "a vector without names",
        outside="assignment by name into a vector without names",
    )
    for _ in range(n):
        k = int(rng.integers(0, 5))
        u.add(
            {
                "v": named(int(rng.integers(0, 6))),
                "names": names(k),
                "values": num(_some_numbers(rng, k, 0.1)),
            },
            "random",
        )
    out += u.cases

    u = unit("v[names]", "read_by_name")
    u.add(
        {"v": v, "names": chr_(["g1", "zz", None])}, "a name that is not there, and a missing one"
    )
    u.add({"v": num([1.0, 2.0]), "names": chr_(["g1"])}, "a vector without names")
    u.add({"v": v, "names": chr_([])}, "no name")
    u.add(
        {"v": num([1.0, 2.0, 3.0], ["g1", "g1", None]), "names": chr_(["g1", None])}, "a name twice"
    )
    for _ in range(n):
        u.add(
            {"v": maybe_named(int(rng.integers(0, 6))), "names": names(int(rng.integers(0, 6)))},
            "random",
        )
    out += u.cases

    thresholds = [0.0, 0.5, -1.0, 1.0, na, nan]
    for fn, slug in (
        ("v[v < t]", "subset_below"),
        ("v[v > t]", "subset_above"),
        ("v[abs(v) > t]", "subset_abs_above"),
    ):
        u = unit(fn, slug)
        w = num([0.9, na, 0.2, nan, -0.7, 0.5], ["g1", "g2", "g3", "g4", "g5", "g6"])
        for t in (0.5, -1.0, 2.0, na):
            u.add({"v": w, "t": num([t])}, "a missing value is kept as a missing element")
        u.add({"v": num([0.9, na, 0.2]), "t": num([0.5])}, "a vector without names")
        u.add({"v": num([], []), "t": num([0.5])}, "an empty vector")
        u.add({"v": num([0.5, 0.5], ["g1", "g2"]), "t": num([0.5])}, "values on the threshold")
        for _ in range(n):
            u.add(
                {
                    "v": maybe_named(int(rng.integers(0, 7))),
                    "t": num([float(rng.choice(thresholds))]),
                },
                "random",
            )
        out += u.cases

    u = unit("names(v) <- names", "set_names")
    plain_v = num([1.0, 2.0, 3.0])
    u.add({"v": plain_v, "names": chr_(["g1", "g2", "g3"])}, "as many names as elements")
    u.add({"v": plain_v, "names": chr_(["g1"])}, "fewer names: the others are missing")
    u.add({"v": plain_v, "names": chr_(["g1", "g2", "g3", "g4"])}, "more names than elements")
    u.add({"v": plain_v, "names": chr_([])}, "no name")
    u.add({"v": NULL, "names": chr_([])}, "NULL and no name")
    u.add({"v": NULL, "names": chr_(["g1"])}, "a name for NULL")
    u.add({"v": num([]), "names": chr_([])}, "an empty vector and no name")
    u.add({"v": plain_v, "names": chr_(["g1", None, "g1"])}, "a name twice and a missing one")
    for _ in range(n):
        k = int(rng.integers(0, 6))
        u.add(
            {"v": num(_some_numbers(rng, k)), "names": names(int(rng.integers(0, k + 2)))}, "random"
        )
    out += u.cases

    for fn, slug in (
        ("union(a, b)", "union"),
        ("intersect(a, b)", "intersect"),
        ("setdiff(a, b)", "setdiff"),
        ("a %in% b", "in"),
    ):
        u = unit(fn, slug)
        u.add({"a": chr_(["g1", "g2", "g1"]), "b": chr_(["g2", "g3"])}, "an element twice")
        u.add({"a": chr_(["g1", None]), "b": chr_([None, "g2"])}, "a missing value on both sides")
        u.add({"a": chr_([]), "b": chr_(["g1"])}, "an empty first set")
        u.add({"a": chr_(["g1"]), "b": chr_([])}, "an empty second set")
        u.add({"a": chr_([]), "b": chr_([])}, "two empty sets")
        for _ in range(n):
            u.add(
                {"a": names(int(rng.integers(0, 7))), "b": names(int(rng.integers(0, 7)))}, "random"
            )
        out += u.cases

    u = unit("unique(a)", "unique")
    u.add({"a": chr_(["g2", "g1", "g2", None, None])}, "first occurrences, in order")
    u.add({"a": chr_([])}, "empty")
    for _ in range(n):
        u.add({"a": names(int(rng.integers(0, 8)))}, "random")
    out += u.cases

    u = unit("c(a, b)", "concatenate_lists")
    one_named = lst([_ones(["g1"])], ["s1"])
    bare = lst([_ones(["g2"]), _ones(["g3"])])
    u.add({"a": one_named, "b": bare}, "names on one side: the other elements get an empty name")
    u.add({"a": bare, "b": bare}, "no names on either side")
    u.add({"a": lst([]), "b": one_named}, "an empty list first")
    u.add({"a": one_named, "b": lst([])}, "an empty list second")
    u.add({"a": lst([]), "b": lst([])}, "two empty lists")
    for _ in range(n):
        a = _neighbourhoods(rng, int(rng.integers(0, 4)))
        b = _neighbourhoods(rng, int(rng.integers(0, 4)), prefix="o")
        if rng.random() < 0.3:
            a = lst(a["items"])
        if rng.random() < 0.3:
            b = lst(b["items"])
        u.add({"a": a, "b": b}, "random")
    out += u.cases

    u = unit("sort(v, decreasing = TRUE)", "sort_decreasing")
    ties = num([0.5, 1.0, 0.5, 0.0, 1.0, 0.5], ["g1", "g2", "g3", "g4", "g5", "g6"])
    u.add({"v": ties}, "equal values keep their order")
    u.add(
        {"v": num([0.5, na, nan, inf, -inf, 0.0], ["g1", "g2", "g3", "g4", "g5", "g6"])},
        "missing values are dropped",
    )
    u.add({"v": num([0.5, 1.0, 0.5])}, "without names")
    u.add({"v": num([], [])}, "empty")
    u.add({"v": num([0.0, -0.0, 0.0], ["g1", "g2", "g3"])}, "zero and minus zero")
    for _ in range(n):
        k = int(rng.integers(0, 9))
        values = [float(rng.choice([0.0, 0.25, 0.5, 1.0, -1.0, na, nan, inf])) for _ in range(k)]
        u.add(
            {"v": num(values, None if rng.random() < 0.15 else _some_names(rng, k, p_dup=0.1))},
            "random",
        )
    out += u.cases

    u = unit('as.numeric(sprintf("%.3f", v))', "three_decimals")
    u.add(
        {"v": num([0.0005, 0.0015, 0.0025, 2.5, 1e10, -0.0004, 0.9995, 1 / 3, 2 / 3])},
        "halves and carries",
    )
    u.add({"v": num([na, nan, inf, -inf])}, "missing and infinite values")
    u.add({"v": num([])}, "empty")
    u.add(
        {"v": num(list(rng.integers(-5000, 5000, size=400) / 10000 + 0.00005))},
        "values on a half, as stored",
    )
    u.add({"v": num(list(rng.normal(size=400)))}, "random")
    out += u.cases

    u = unit("a:b", "colon")
    for a, b in ((1, 3), (1, 1), (1, 0), (3, 1), (1, -1), (0, 2)):
        u.add({"a": num([float(a)]), "b": num([float(b)])}, f"{a}:{b}")
    u.add({"a": num([1.0]), "b": NULL}, "1:NULL")
    out += u.cases

    u = unit("write.table(x)", "write_table")
    u.add({"x": num([0.5, 0.25], ["g1", "g2"])}, "a vector with names")
    u.add({"x": num([0.5, 0.25])}, "a vector without names: the rows are numbered")
    u.add({"x": num([0.5, 0.25], ["g1", "g1"])}, "a name twice")
    u.add({"x": num([0.5, 0.25], ["g1", None])}, "a missing name")
    u.add({"x": num([0.5, 0.25], [None, None])}, "only missing names")
    u.add({"x": num([0.5, 0.25], ["", ""])}, "empty names")
    u.add({"x": num([0.5, 0.25], ["g1", ""])}, "one empty name")
    u.add({"x": num([], [])}, "an empty vector with names")
    u.add({"x": num([])}, "an empty vector")
    u.add(
        {
            "x": num(
                [na, nan, inf, -inf, 1e5, 123456.0, 1e-5, 0.0001, 1e15, 1e16, 0.1 + 0.2],
                [f"g{i}" for i in range(11)],
            )
        },
        "numbers of every size",
    )
    u.add({"x": chr_(["g1", "g2"])}, "strings")
    u.add({"x": chr_([])}, "no string")
    u.add({"x": chr_(["g1", None])}, "a missing string")
    u.add({"x": mat([[1, 0.5], [na, 1]], ["g1", "g2"], ["g1", "g2"])}, "a matrix")
    u.add({"x": mat([[1, 0.5], [na, 1]], None, None)}, "a matrix without names")
    u.add({"x": mat([[1, 0.5], [na, 1]], ["g1", "g2"], None)}, "a matrix with row names only")
    u.add({"x": mat([[1, 0.5], [na, 1]], None, ["g1", "g2"])}, "a matrix with column names only")
    u.add({"x": mat(np.empty((2, 0)), ["g1", "g2"], None)}, "a matrix without columns")
    u.add({"x": mat(np.empty((1, 0)), ["g1"], None)}, "one row and no column")
    u.add({"x": mat([[1.0, 0.5]], ["g1"], ["g1", "g2"])}, "a matrix of one row")
    u.add({"x": mat([[1.0], [0.5]], ["g1", "g2"], ["g1"])}, "a matrix of one column")
    for _ in range(n):
        k = int(rng.integers(1, 8))
        values = list(rng.normal(size=k) * 10.0 ** rng.integers(-12, 12, size=k))
        if rng.random() < 0.5:
            u.add({"x": num(values, _some_names(rng, k))}, "random vector")
        else:
            m = rng.normal(size=(k, 3)) * 10.0 ** rng.integers(-12, 12, size=(k, 3))
            m[rng.random((k, 3)) < 0.2] = na
            u.add({"x": mat(m, _some_names(rng, k), ["c1", "c2", "c3"])}, "random matrix")
    out += u.cases

    u = unit("read.table(file)", "read_table")

    def text(rows: Sequence[tuple[str, str]]) -> Json:
        return chr_([_connectivity_text(rows)])

    u.add({"text": text([("g1", "0.5"), ("g2", "0.25")])}, "two genes")
    u.add({"text": text([("g1", "0.5")])}, "one gene")
    u.add({"text": chr_(["\tx\n"])}, "no gene")
    u.add({"text": chr_([""])}, "an empty file")
    u.add(
        {"text": text([("g1", "NA"), ("g2", "Inf"), ("g3", "-Inf"), ("g4", "NaN")])},
        "NA, NaN and Inf",
    )
    u.add({"text": text([("g1", "0.5"), ("g1", "0.25")])}, "a gene twice")
    u.add({"text": text([("007", "0.5"), ("12", "0.25")])}, "names that are whole numbers")
    u.add({"text": text([("1.50", "0.5"), ("100000", "0.25")])}, "names that are decimal numbers")
    u.add({"text": text([("1e3", "0.5"), ("1000", "0.25")])}, "two names that are the same number")
    u.add({"text": text([("T", "0.5"), ("F", "0.25")])}, "names that are logical words")
    u.add({"text": text([("TRUE", "0.5"), ("true", "0.25")])}, "a lowercase logical word is a name")
    u.add({"text": text([("NA", "0.5"), ("g1", "0.25")])}, "a gene called NA")
    u.add({"text": text([("0x10", "0.5"), ("16", "0.25")])}, "a name in hexadecimal")
    u.add(
        {"text": text([("3000000000", "0.5"), ("1", "0.25")])},
        "a whole number too large for an integer",
    )
    u.add({"text": text([("-1", "0.5"), ("+1", "0.25")])}, "signed numbers as names")
    u.add(
        {"text": text([("g1", "9.118641810491681"), ("g2", "0.1234567890123456789")])},
        "values of 17 digits and more",
    )
    u.add(
        {"text": text([("g1", "1e-5"), ("g2", "1E5"), ("g3", "5e-324"), ("g4", "1e400")])},
        "exponents",
    )
    for _ in range(n):
        k = int(rng.integers(1, 7))
        if rng.random() < 0.25:
            pool = ["1", "02", "3.0", "4e1", "5", "0.5", "1000000", "7"]
        else:
            pool = POOL
        rows = [
            (str(name), str(rng.choice([*CONNECTIVITY_VALUES, repr(float(rng.normal()))])))
            for name in _some_names(rng, k, pool)
        ]
        u.add({"text": text(rows)}, "random")
    out += u.cases
    return out


RANDOM = {
    "getO": 60,
    "getOmax": 80,
    "removeNode": 40,
    "removeByName": 50,
    "get.matrix.format": 80,
    "estimate.nodes.degree": 50,
    "nodes.correlation": 150,
    "nodes.glm": 40,
    "nodes.fastLm": 24,
    "SEARCH": 60,
    "write.data.files": 60,
    ".SEARCH.randomise.inner": 16,
    "SEARCH.randomise": 24,
    "summarise.iterations": 80,
    "rules": 30,
}
"""Random cases per function in the set kept for the tests, on top of the systematic ones."""


def all_cases(scale: float = 1.0) -> list[Json]:
    def n(name: str) -> int:
        return max(1, round(RANDOM[name] * scale))

    cases = [
        *cases_getO(n("getO")),
        *cases_getOmax(n("getOmax")),
        *cases_removeNode(n("removeNode")),
        *cases_removeByName(n("removeByName")),
        *cases_get_matrix_format(n("get.matrix.format")),
        *cases_estimate_nodes_degree(n("estimate.nodes.degree")),
        *cases_nodes_correlation(n("nodes.correlation")),
        *_regression_cases("nodes.glm", "nodes_glm", n("nodes.glm"), 1008),
        *_regression_cases("nodes.fastLm", "nodes_fastLm", n("nodes.fastLm"), 1009),
        *cases_SEARCH(n("SEARCH")),
        *cases_write_data_files(n("write.data.files")),
        *cases_randomise_inner(n(".SEARCH.randomise.inner")),
        *cases_SEARCH_randomise(n("SEARCH.randomise")),
        *cases_summarise_iterations(n("summarise.iterations")),
        *cases_rules(n("rules")),
    ]
    ids = [case["id"] for case in cases]
    assert len(set(ids)) == len(ids)
    return [{k: v for k, v in case.items() if v is not None} for case in cases]


# ---- running ------------------------------------------------------------------------------------
def answer_in_r(cases_file: Path, out_file: Path) -> None:
    script = REPO / "validation" / "r_parity" / "r" / "r2014_units.R"
    subprocess.run(
        ["Rscript", str(script), str(REPO), str(cases_file), str(out_file)],
        check=True,
        env=os.environ.copy(),
    )


def check(
    cases: list[Json], answers: Json, regression_tolerance: float
) -> tuple[list[Json], dict[str, Counter[str]]]:
    """Every case against R's answer: one row per case, and the counts per function."""
    by_id = {r["id"]: r for r in answers["results"]}
    target = platform.matching(
        answers["info"]["r_version"], int(answers["info"]["long_double_bytes"])
    )
    if target is None:
        raise SystemExit(f"the port has no platform for {answers['info']['r_version']}")
    rows: list[Json] = []
    counts: dict[str, Counter[str]] = {}
    with platform.use(target):
        for case in cases:
            r = by_id[case["id"]]
            count = counts.setdefault(case["fn"], Counter())
            if r.get("skipped"):
                outcome, detail = "skipped", r["skipped"]
            else:
                port = run_port(case)
                found = compare(case, r, port, tolerance_of(case, regression_tolerance))
                if "outside" in case:
                    refused = "unported" in port
                    outcome = "refused, as stated" if refused else "DIFFERENT"
                    detail = case["outside"] if refused else "the port was expected to refuse"
                elif found:
                    outcome, detail = "DIFFERENT", "; ".join(found[:3])
                else:
                    outcome = "equal, R stops" if r.get("error") is not None else "equal"
                    detail = r.get("error") or ""
            count[outcome] += 1
            rows.append(
                {
                    "case": case["id"],
                    "function": case["fn"],
                    "group": case["group"],
                    "what": case.get("note", ""),
                    "outcome": outcome,
                    "detail": detail,
                }
            )
    return rows, counts


def summary_table(info: Json, counts: dict[str, Counter[str]], groups: dict[str, str]) -> str:
    """The counts per function as a Markdown table."""
    columns = ("equal", "equal, R stops", "refused, as stated", "DIFFERENT", "skipped")
    lines = [
        f"Unit cases of `search.r2014` against {info['r_version']} on {info['platform']}, "
        f"{info['date']}.",
        "",
        "| function or rule of R | cases | equal | of which R stops | refused by the port, "
        "as stated | differ | skipped |",
        "|---|---|---|---|---|---|---|",
    ]

    def row(label: str, count: Counter[str]) -> str:
        equal = count["equal"] + count["equal, R stops"]
        cells = [sum(count.values()), equal, *(count[c] for c in columns[1:])]
        return f"| {label} | " + " | ".join(str(c) for c in cells) + " |"

    for group in (PACKAGE, RULES):
        subtotal: Counter[str] = Counter()
        for fn, count in counts.items():
            if groups[fn] == group:
                subtotal.update(count)
                lines.append(row(f"`{fn}`", count))
        lines.append(row(f"**{group}**", subtotal))
    return "\n".join(lines) + "\n"


def report(cases: list[Json], answers: Json, args: argparse.Namespace) -> int:
    rows, counts = check(cases, answers, args.tolerance)
    print(json.dumps(answers["info"]))
    total: Counter[str] = Counter()
    for fn, count in counts.items():
        total.update(count)
        print(f"{fn:34s} {sum(count.values()):6d}  {dict(count)}")
    print(f"{'all':34s} {sum(total.values()):6d}  {dict(total)}")
    for row in rows:
        if row["outcome"] == "DIFFERENT":
            print(f"  {row['case']} ({row['what']}): {row['detail']}")
    if args.report:
        with args.report.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    if args.summary:
        groups = {case["fn"]: case["group"] for case in cases}
        args.summary.write_text(summary_table(answers["info"], counts, groups))
    return total["DIFFERENT"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "check"):
        command = sub.add_parser(name)
        command.add_argument("out_dir", type=Path)
        command.add_argument("--report", type=Path)
        command.add_argument("--summary", type=Path)
        command.add_argument("--tolerance", type=float, default=r2014_cases.REGRESSION_TOLERANCE)
        if name == "run":
            command.add_argument(
                "--scale", type=float, default=1.0, help="multiplies the random cases"
            )
    golden = sub.add_parser("golden")
    golden.add_argument("out_dir", type=Path)
    golden.add_argument("data_dir", type=Path)
    args = parser.parse_args()
    cases_file = args.out_dir / "cases.json"
    answers_file = args.out_dir / "answers.json"
    if args.command == "run":
        args.out_dir.mkdir(parents=True, exist_ok=True)
        cases_file.write_text(json.dumps(all_cases(args.scale)))
        answer_in_r(cases_file, answers_file)
    cases = json.loads(cases_file.read_text())
    answers = json.loads(answers_file.read_text())
    if args.command in ("run", "check"):
        sys.exit(1 if report(cases, answers, args) else 0)
    target = platform.matching(
        answers["info"]["r_version"], int(answers["info"]["long_double_bytes"])
    )
    slug = "r31_x86" if target is platform.R_3_1_X86 else "r44_arm64"
    path = args.data_dir / f"units_{slug}.json.gz"
    content = {
        "info": answers["info"],
        "cases": cases,
        "answers": {r["id"]: r for r in answers["results"]},
    }
    with gzip.GzipFile(path, "wb", compresslevel=9, mtime=0) as handle:
        handle.write(json.dumps(content, separators=(",", ":")).encode())
    print(f"{path}: {len(cases)} cases, {path.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
