"""Verbatim port against the R package: reading the cases and the oracle, running the port, and
comparing the two exactly.

"Exactly" means: the same names in the same order, missing names where R has them, and for
every number the same double, bit for bit, with NA and NaN kept apart. Quantities that cannot be
bit-identical (the coefficients and p-values of the regression functions, which come out of
different linear-algebra and distribution routines) are compared with a stated tolerance.
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from search.errors import SearchError
from search.r2014 import core, randomise, summarise
from search.r2014.rsem import (
    NA_REAL,
    NamedVector,
    RError,
    RMatrix,
    UnportedState,
    is_na_real,
)

Json = dict[str, Any]


# ---- from JSON --------------------------------------------------------------------------------
def numbers(values: list[str]) -> np.ndarray:
    return np.array([NA_REAL if v == "NA" else float(v) for v in values], dtype=np.float64)


def case_matrix(case: Json, datasets: dict[str, RMatrix]) -> RMatrix:
    d = case["data"]
    if "dataset" in d:
        base = datasets[d["dataset"]]
        values = base.values
        rownames = list(base.rownames or [])
        colnames = list(base.colnames or [])
        if d.get("rows") is not None:
            values = values[d["rows"]]
            rownames = [rownames[i] for i in d["rows"]]
        if d.get("cols") is not None:
            values = values[:, d["cols"]]
            colnames = [colnames[i] for i in d["cols"]]
        x = RMatrix(values.copy(), rownames, colnames)
    else:
        values = np.array(
            [NA_REAL if v is None or v == "NA" else float(v) for v in d["values"]], dtype=np.float64
        ).reshape(d["nrow"], d["ncol"])
        x = RMatrix(values, list(d["rownames"]), list(d["colnames"]) if d.get("colnames") else None)
    for r, c in case.get("na_cells", []):
        x.values[r, c] = NA_REAL
    return x


def python_args(r_args: Json) -> Json:
    """R argument names to the port's: dots become underscores."""
    return {name.replace(".", "_"): value for name, value in r_args.items()}


def case_covariates(case: Json) -> Any:
    cv = case.get("covariates")
    if cv is None:
        return None
    from search.r2014.glm import Covariates

    columns: dict[str, list[Any]] = {}
    for col in cv["columns"]:
        if col["type"] == "numeric":
            columns[col["name"]] = [NA_REAL if v is None else float(v) for v in col["values"]]
        else:
            columns[col["name"]] = [None if v is None else str(v) for v in col["values"]]
    return Covariates(columns)


# ---- running the port -------------------------------------------------------------------------
@dataclass
class PortRun:
    result: Any = None
    error: str | None = None
    unported: str | None = None
    extra: Json = field(default_factory=dict)


def run_search(case: Json, datasets: dict[str, RMatrix]) -> PortRun:
    x = case_matrix(case, datasets)
    kwargs = {**python_args(case.get("args", {})), **python_args(case.get("assoc_args", {}))}
    try:
        universe = case["universe"] if case.get("universe") is not None else list(x.rownames or [])
        res = core.SEARCH(
            x=x, Alist=case["seeds"], Blist=universe, covariates=case_covariates(case), **kwargs
        )
        try:
            res.nodes_degree = core.estimate_nodes_degree(res.correlationMatrix)
        except SearchError:
            res.nodes_degree = {}
        return PortRun(result=res)
    except RError as exc:
        return PortRun(error=exc.message)
    except UnportedState as exc:
        return PortRun(unported=str(exc))


def run_randomise(case: Json, datasets: dict[str, RMatrix], work: Path | None = None) -> PortRun:
    x = case_matrix(case, datasets)
    kwargs = {**python_args(case.get("args", {})), **python_args(case.get("assoc_args", {}))}
    if case.get("sample_kind"):
        kwargs["sample_kind"] = case["sample_kind"]
    out = PortRun()
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) if work is None else work
        real_dir = case.get("real_dir", "real/")
        universe = case["universe"] if case.get("universe") is not None else list(x.rownames or [])
        runs = [
            {
                "dir": real_dir,
                "seeds": case["seeds"],
                "n.boot": case["args"].get("n.boot"),
                "n.boot.start": case["args"].get("n.boot.start"),
            },
            *case.get("null_runs", []),
        ]
        summaries: Json = {}
        try:
            for r in runs:
                a = dict(kwargs)
                if r.get("n.boot") is not None:
                    a["n_boot"] = r["n.boot"]
                if r.get("n.boot.start") is not None:
                    a["n_boot_start"] = r["n.boot.start"]
                current = r["dir"]
                summaries[current] = randomise.SEARCH_randomise(
                    data=x,
                    seed_nodes=r["seeds"],
                    universe_nodes=universe,
                    covariates=case_covariates(case),
                    output_directory=f"{root}/{current}",
                    **a,
                )
        except RError as exc:
            out.error = exc.message
            out.extra["error_run"] = current
        except UnportedState as exc:
            out.unported = str(exc)
        out.extra["summaries"] = summaries
        if out.error is None and out.unported is None and case.get("summarise"):
            s = case["summarise"]
            try:
                out.extra["summarise"] = summarise.summarise_iterations(
                    output_directory=f"{root}/{real_dir}",
                    max_iter=s["max.iter"],
                    n_boot=s["n.boot"],
                    output_directory_NULL=f"{root}/null/",
                    n_boot_NULL=s["n.boot.NULL"],
                    batches_NULL=s["batches.NULL"],
                )
            except RError as exc:
                out.extra["summarise_error"] = exc.message
            except UnportedState as exc:
                out.unported = str(exc)
        out.extra["files"] = randomise.read_output_files(root)
    return out


# ---- exact comparison -------------------------------------------------------------------------
def same_numbers(r: np.ndarray, py: np.ndarray, tolerance: float = 0.0) -> tuple[int, float]:
    """Number of positions that differ and the largest absolute difference. With tolerance 0
    the two doubles must have the same bits (and the same kind of missing value)."""
    r = np.ascontiguousarray(r, dtype=np.float64).ravel()
    py = np.ascontiguousarray(py, dtype=np.float64).ravel()
    if r.shape != py.shape:
        return max(len(r), len(py)), float("inf")
    both_missing = np.isnan(r) & np.isnan(py) & (is_na_real(r) == is_na_real(py))
    with np.errstate(invalid="ignore"):
        diff = np.where(both_missing, 0.0, np.abs(r - py))
    diff = np.where(np.isnan(diff), np.inf, diff)
    diff = np.where(np.isinf(r) & np.isinf(py) & (np.sign(r) == np.sign(py)), 0.0, diff)
    if tolerance == 0.0:
        bad = ~((r.view(np.int64) == py.view(np.int64)) | both_missing)
    else:
        # relative to the size of the number, and absolute below 1. Numbers beyond 1e12 are the
        # t statistics of a seed regressed on itself: a division by rounding noise on both sides
        scale = np.maximum(1.0, np.maximum(np.abs(r), np.abs(py)))
        with np.errstate(invalid="ignore"):
            noise = (np.abs(r) > 1e12) & (np.abs(py) > 1e12) & (np.sign(r) == np.sign(py))
            bad = (diff > tolerance * scale) & ~noise
        diff = np.where(noise, 0.0, diff / np.where(np.isfinite(scale), scale, 1.0))
    return int(bad.sum()), float(diff.max()) if len(diff) else 0.0


class Differences:
    """Collects what differs in one case, quantity by quantity."""

    def __init__(self) -> None:
        self.items: list[str] = []
        self.compared = 0
        self.values = 0
        self.max_diff = 0.0

    def check(self, what: str, ok: bool, detail: str = "") -> bool:
        self.compared += 1
        if not ok:
            self.items.append(f"{what}: {detail}" if detail else what)
        return ok

    def names(self, what: str, r: list[Any] | None, py: list[Any] | None) -> bool:
        self.values += len(r or [])
        return self.check(what, r == py, f"R {_short(r)} / port {_short(py)}")

    def numbers(self, what: str, r: np.ndarray, py: np.ndarray, tolerance: float = 0.0) -> bool:
        n_bad, largest = same_numbers(r, py, tolerance)
        self.values += len(np.ravel(r))
        if n_bad == 0 and np.isfinite(largest):
            self.max_diff = max(self.max_diff, largest)
        return self.check(
            what, n_bad == 0, f"{n_bad} of {len(np.ravel(r))} values, largest {largest:.3g}"
        )

    def vector(self, what: str, r: Json, py: NamedVector | None, tolerance: float = 0.0) -> None:
        if py is None:
            self.check(what, bool(r.get("is_null")), "the port has NULL, R has a vector")
            return
        r_names = r["names"] if r["has_names"] else None
        self.names(f"{what}: names", r_names, py.names)
        self.numbers(f"{what}: values", numbers(r["values"]), py.values, tolerance)

    def coefficient_table(self, what: str, r: Json, py: RMatrix, tolerance: float) -> None:
        """``other.data``: every coefficient of every model."""
        if not self.check(
            f"{what}: dim",
            list(r["dim"]) == [py.nrow, py.ncol],
            f"R {r['dim']} / port {[py.nrow, py.ncol]}",
        ):
            return
        self.names(f"{what}: rownames", r["rownames"], py.rownames)
        self.names(f"{what}: colnames", r["colnames"], py.colnames)
        r_values = numbers(r["values"]).reshape(py.ncol, py.nrow).T.copy()
        self.coefficients(
            what, r_values, py.values.copy(), py.rownames or [], py.colnames or [], tolerance
        )

    def coefficients(
        self,
        what: str,
        r_values: np.ndarray,
        py_values: np.ndarray,
        rownames: list[Any],
        colnames: list[Any],
        tolerance: float,
    ) -> None:
        """Rows are ``<seed>.<term>``, columns ``<gene>.<statistic>``. Where a seed is regressed
        on itself the fit is exact, and its standard errors, t statistics and p-values are
        divisions by rounding noise: there only the estimates are compared, as the 0 or 1 they
        stand for."""
        for j, column in enumerate(colnames):
            is_estimate = column.endswith(".Estimate")
            pieces = 2 if column.endswith((".t.value", ".p.value")) else 1
            gene = column[: -len(".Estimate")] if is_estimate else column.rsplit(".", pieces)[0]
            for i, row in enumerate(rownames):
                if row.startswith(gene + "."):
                    if is_estimate:
                        r_values[i, j] = round(r_values[i, j], 9)
                        py_values[i, j] = round(py_values[i, j], 9)
                    else:
                        r_values[i, j] = py_values[i, j] = 0.0
        self.numbers(f"{what}: values", r_values, py_values, tolerance)

    def matrix(self, what: str, r: Json, py: RMatrix, tolerance: float = 0.0) -> None:
        if not self.check(
            f"{what}: dim",
            list(r["dim"]) == [py.nrow, py.ncol],
            f"R {r['dim']} / port {[py.nrow, py.ncol]}",
        ):
            return
        self.names(f"{what}: rownames", r["rownames"] if r["has_rownames"] else None, py.rownames)
        self.names(f"{what}: colnames", r["colnames"] if r["has_colnames"] else None, py.colnames)
        r_values = (
            numbers(r["values"]).reshape(py.ncol, py.nrow).T if py.nrow * py.ncol else np.empty(0)
        )
        self.numbers(f"{what}: values", r_values, py.values, tolerance)


def _short(v: Any) -> str:
    text = repr(v)
    return text if len(text) <= 160 else text[:157] + "..."


def compare_search(r: Json, run: PortRun, tolerance: float = 0.0) -> Differences:
    """Everything ``SEARCH()`` returns, plus the node degrees computed from it."""
    d = Differences()
    if r.get("error") is not None or run.error is not None:
        d.check("error", r.get("error") == run.error, f"R {r.get('error')!r} / port {run.error!r}")
        return d
    res: core.SearchResult = run.result
    d.names("removed", r["removed"], res.removed)
    cl = r["correlationLists"]
    d.names(
        "correlationLists: names",
        cl["names"] if cl["has_names"] else None,
        res.correlationLists.names,
    )
    if d.check(
        "correlationLists: length",
        len(cl["items"]) == len(res.correlationLists),
        f"R {len(cl['items'])} / port {len(res.correlationLists)}",
    ):
        for k, (r_item, py_item) in enumerate(
            zip(cl["items"], res.correlationLists.items, strict=True)
        ):
            d.vector(f"correlationLists[[{k + 1}]]", r_item, py_item, tolerance)
    if d.check(
        "D: iterations", list(r["D"]) == list(res.D), f"R {list(r['D'])} / port {list(res.D)}"
    ):
        for it in res.D:
            d.vector(f"D[[{it}]]", r["D"][it], res.D[it])
    if d.check(
        "network.growth: iterations",
        list(r["network.growth"]) == list(res.network_growth),
        f"R {list(r['network.growth'])} / port {list(res.network_growth)}",
    ):
        for it in res.network_growth:
            d.names(
                f"network.growth[[{it}]]", r["network.growth"][it]["values"], res.network_growth[it]
            )
    for key, port in (
        ("O", res.O),
        ("correlationMatrix", res.correlationMatrix),
        ("other.data", res.other_data),
    ):
        if d.check(
            f"{key}: iterations",
            list(r[key]) == list(port),
            f"R {list(r[key])} / port {list(port)}",
        ):
            for it in port:
                if key == "other.data" and tolerance > 0:
                    d.coefficient_table(f"{key}[[{it}]]", r[key][it], port[it], tolerance)
                else:
                    d.matrix(
                        f"{key}[[{it}]]", r[key][it], port[it], tolerance if key != "O" else 0.0
                    )
    if d.check(
        "nodes.degree: iterations",
        list(r["nodes.degree"]) == list(res.nodes_degree),
        f"R {list(r['nodes.degree'])} / port {list(res.nodes_degree)}",
    ):
        for it in res.nodes_degree:
            d.vector(f"nodes.degree[[{it}]]", r["nodes.degree"][it], res.nodes_degree[it])
    return d


def _parse_table(text: str) -> tuple[list[str], list[str], np.ndarray]:
    lines = [line for line in text.split("\n") if line != ""]
    columns = lines[0].split("\t")[1:]
    rows: list[str] = []
    values: list[list[float]] = []
    for line in lines[1:]:
        fields = line.split("\t")
        rows.append(fields[0])
        values.append([NA_REAL if f == "NA" else float(f) for f in fields[1:]])
    return rows, columns, np.array(values, dtype=np.float64).reshape(len(rows), len(columns))


def compare_randomise(r: Json, run: PortRun, tolerance: float = 0.0) -> Differences:
    """What ``SEARCH.randomise`` returns and writes, and what ``summarise.iterations`` returns.

    The files are compared character for character once uncompressed. With a tolerance (the
    regression functions) the files of coefficients are compared number by number instead.
    """
    d = Differences()
    d.check("error", r.get("error") == run.error, f"R {r.get('error')!r} / port {run.error!r}")
    if r.get("error") is not None:
        d.check(
            "error: run",
            r.get("error_run") == run.extra.get("error_run"),
            f"R {r.get('error_run')} / port {run.extra.get('error_run')}",
        )
    r_files: Json = r.get("files") or {}
    py_files: Json = run.extra.get("files") or {}
    d.names("files written", sorted(r_files), sorted(py_files))
    for name in sorted(set(r_files) & set(py_files)):
        if tolerance > 0 and "other_data" in name:
            r_rows, r_cols, r_values = _parse_table(r_files[name])
            py_rows, py_cols, py_values = _parse_table(py_files[name])
            if d.check(f"file {name}: names", (r_rows, r_cols) == (py_rows, py_cols)):
                d.coefficients(f"file {name}", r_values, py_values, py_rows, py_cols, tolerance)
        else:
            d.values += r_files[name].count("\n")
            d.check(
                f"file {name}",
                r_files[name] == py_files[name],
                _first_line_diff(r_files[name], py_files[name]),
            )
    d.check(
        "returned summaries",
        _plain(r.get("summaries")) == _plain(run.extra.get("summaries")),
        f"R {_short(r.get('summaries'))} / port {_short(run.extra.get('summaries'))}",
    )
    d.check(
        "summarise.iterations: error",
        r.get("summarise_error") == run.extra.get("summarise_error"),
        f"R {r.get('summarise_error')!r} / port {run.extra.get('summarise_error')!r}",
    )
    r_sum = r.get("summarise")
    py_sum = run.extra.get("summarise")
    if (
        d.check("summarise.iterations: returned", (r_sum is None) == (py_sum is None))
        and r_sum is not None
    ):
        for key, port in (("real", py_sum.real), ("null", py_sum.null)):
            if d.check(
                f"summarise {key}: iterations",
                list(r_sum[key]) == list(port),
                f"R {list(r_sum[key])} / port {list(port)}",
            ):
                for it in port:
                    d.matrix(f"summarise {key}[[{it}]]", r_sum[key][it], port[it])
    return d


def _plain(v: Any) -> Any:
    """Summaries as nested lists of floats, whatever container they came in. An empty list of
    R is written ``[]`` whether or not it had names, so every empty container is the same."""
    if isinstance(v, dict):
        return {k: _plain(x) for k, x in v.items()} or None
    if isinstance(v, (list, tuple)):
        return [_plain(x) for x in v] or None
    if v is None:
        return None
    return float(v)


def _first_line_diff(a: str, b: str) -> str:
    la, lb = a.split("\n"), b.split("\n")
    for i, (x, y) in enumerate(zip(la, lb, strict=False)):
        if x != y:
            return f"line {i + 1}: R {x[:120]!r} / port {y[:120]!r}"
    return f"{len(la)} lines in R, {len(lb)} in the port"
