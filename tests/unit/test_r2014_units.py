"""Unit tests of the verbatim port (``search.r2014``): every function of the 2014 package, on
its own, against the answer the same function gave in R.

``test_r2014_golden.py`` compares whole runs. A whole run can be right while a function is
wrong, when the error never reaches the result or two errors cancel, and it reaches the
helpers only with the inputs a run happens to produce. Here every function of the package is
called directly, in R and in the port, on inputs chosen for it: ordinary ones and those at the
edges (nothing, one element, no names, missing names, a name twice, missing values, ties).
The rules of R the port had to write down (what ``v[names] <- values`` does with a name that
is not there, how ``write.table`` writes a number) are tested the same way, each as the R
expression it stands for.

The cases and R's answers are in ``tests/data/r2014/units_r44_arm64.json.gz``, written by
``validation/r_parity/run_r2014_units.sh`` from the unmodified package in R 4.4.3; the tests
need no R. One test is one case: its name says the function and the number of the case, and
``pytest -k getOmax`` runs those of one function. What is compared and how is in
``docs/r2014_tests.md``.
"""

from __future__ import annotations

import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "validation" / "r_parity" / "py"))

import r2014_units as units  # noqa: E402

from search.r2014 import association, core, platform, randomise, rio, summarise  # noqa: E402
from search.r2014.rformat import format_real  # noqa: E402
from search.r2014.rsem import NA_REAL, NamedVector, RList, RMatrix  # noqa: E402

Json = dict[str, Any]

with gzip.open(REPO / "tests" / "data" / "r2014" / "units_r44_arm64.json.gz", "rt") as _handle:
    GOLDEN: Json = json.load(_handle)
CASES: list[Json] = GOLDEN["cases"]
ANSWERS: dict[str, Json] = GOLDEN["answers"]
PLATFORM = platform.R_4_ARM64

REGRESSION_TOLERANCE = 1e-6
"""Coefficients and p-values of ``nodes.glm`` and ``nodes.fastLm``, relative: their linear
algebra is not R's and differs between machines. Everything else must be the same bit for bit."""

MINIMUM_CASES = {
    "getO": 50,
    "getOmax": 60,
    "removeNode": 30,
    "removeByName": 40,
    "get.matrix.format": 60,
    "estimate.nodes.degree": 40,
    "nodes.correlation": 150,
    "nodes.glm": 50,
    "nodes.fastLm": 30,
    "SEARCH": 50,
    "write.data.files": 60,
    ".SEARCH.randomise.inner": 15,
    "SEARCH.randomise": 25,
    "summarise.iterations": 80,
}


def _port_answer(case: Json) -> Json:
    with platform.use(PLATFORM):
        return units.run_port(case)


def _differences(case: Json) -> list[str]:
    """Everything in which the port's answer to a case is not the one recorded from R."""
    return units.compare(
        case,
        ANSWERS[case["id"]],
        _port_answer(case),
        units.tolerance_of(case, REGRESSION_TOLERANCE),
    )


def _cases_of(fn: str) -> list[Json]:
    return [case for case in CASES if case["fn"] == fn and "outside" not in case]


# ---- one test per case --------------------------------------------------------------------------
def test_the_answers_come_from_the_r_the_port_stands_for() -> None:
    info = GOLDEN["info"]
    assert platform.matching(info["r_version"], int(info["long_double_bytes"])) is PLATFORM
    assert not [a["id"] for a in ANSWERS.values() if a.get("skipped")], "R skipped some cases"


@pytest.mark.parametrize("case", CASES, ids=[case["id"] for case in CASES])
def test_the_function_returns_what_it_returns_in_r(case: Json) -> None:
    """The value returned, the files written, or the message of the error R stops with."""
    if "outside" in case:
        # an input the port does not model: it must refuse it, not answer something else
        assert "unported" in _port_answer(case), case["outside"]
        return
    found = _differences(case)
    assert not found, f"{case['fn']} ({case['note']}): " + "; ".join(found)


# ---- the set of cases is complete ---------------------------------------------------------------
def _functions_of_the_package() -> set[str]:
    names: set[str] = set()
    for source in (REPO / "reference" / "r_sources" / "SEARCH" / "R").glob("*.R"):
        names.update(re.findall(r"^([A-Za-z.][\w.]*) <- function", source.read_text(), re.M))
    return names


def test_every_function_of_the_package_has_cases_or_is_listed_as_not_ported() -> None:
    package = _functions_of_the_package()
    assert len(package) == 19
    tested = {case["fn"] for case in CASES if case["group"] == units.PACKAGE}
    assert tested | set(units.NOT_PORTED) == package
    assert not tested & set(units.NOT_PORTED)
    assert tested == set(MINIMUM_CASES) == set(units.PORT) - _rules()


def _rules() -> set[str]:
    return {case["fn"] for case in CASES if case["group"] == units.RULES}


def test_every_rule_of_r_the_port_wrote_down_has_cases() -> None:
    assert _rules() == set(units.PORT) - set(MINIMUM_CASES)
    assert len(_rules()) == 18


def test_every_function_has_enough_cases_and_some_in_which_r_stops() -> None:
    count = Counter(case["fn"] for case in CASES)
    assert not {fn: count[fn] for fn, least in MINIMUM_CASES.items() if count[fn] < least}
    stops = Counter(case["fn"] for case in CASES if ANSWERS[case["id"]].get("error") is not None)
    # the functions that can stop are tested where they do; the others never raise an error
    can_stop = set(MINIMUM_CASES) - {
        "getO",
        "getOmax",
        "removeNode",
        "removeByName",
        "estimate.nodes.degree",
    }
    assert not {fn for fn in can_stop if stops[fn] < 3}
    assert sum(stops.values()) >= 150


def test_the_stored_cases_are_the_ones_the_generator_makes() -> None:
    """The file is not edited by hand: its systematic cases are, one by one, those of
    ``r2014_units.all_cases``, and it has as many random ones for every function."""

    def systematic(cases: list[Json]) -> list[tuple[str, str, str]]:
        return [(c["id"], c["fn"], c["note"]) for c in cases if not _by_chance(c)]

    def _by_chance(case: Json) -> bool:
        return bool(case["note"].startswith(("random", "found by search")))

    made = units.all_cases()
    assert systematic(made) == systematic(CASES)
    assert Counter(c["fn"] for c in made) == Counter(c["fn"] for c in CASES)


# ---- the defects of the package are in the cases, as R shows them -------------------------------
def _fields(value: Json) -> Json:
    return dict(zip(value["names"], value["items"], strict=True))


def _neighbourhoods(r: Json) -> list[Json]:
    value = _fields(r["value"])
    lists = value["CL"] if "CL" in value else value["correlationLists"]
    return list(lists["items"])


def _r1(case: Json, r: Json) -> bool:
    # g1 is flat and g3 is not: the prefilter removed g3 and kept g1
    names = _fields(r["value"])["D"]["items"][0]["names"]
    return "g1" in names and "g3" not in names


def _r2(case: Json, r: Json) -> bool:
    # identical lists overlap completely, and the last one gets 0
    values = [units.dec_number(v) for v in r["value"]["v"]]
    return values[0] == 1.0 and values[-1] == 0.0


def _r3(case: Json, r: Json) -> bool:
    # a missing name among the seeds added
    return any(None in added["v"] for added in _fields(r["value"])["network.growth"]["items"])


def _r4(case: Json, r: Json) -> bool:
    # a gene removed from the universe is still a member of a neighbourhood
    removed = set(_fields(r["value"])["removed"]["v"])
    return any(removed & set(item["names"] or []) for item in _neighbourhoods(r))


def _r5(case: Json, r: Json) -> bool:
    return any(len(item["v"]) == 1 and item["names"] is None for item in _neighbourhoods(r))


def _r6(case: Json, r: Json) -> bool:
    return len(_neighbourhoods(r)) < len(case["args"]["seeds"]["v"])


def _r7(case: Json, r: Json) -> bool:
    if case["fn"] == "SEARCH":
        return bool(r["error"] == "subscript out of bounds")
    return any(None in (item["names"] or []) for item in _neighbourhoods(r))


def _r8(case: Json, r: Json) -> bool:
    return bool(r["error"] == "subscript out of bounds")


def _r9(case: Json, r: Json) -> bool:
    # a column of the summary that no replicate filled, with more files than columns filled
    table = _fields(_fields(r["value"])["real"])["1"]
    nrow = table["dim"][0]
    columns = [table["v"][j * nrow : (j + 1) * nrow] for j in range(table["dim"][1])]
    empty_real = any(set(column) == {"NA"} for column in columns[:-2])
    null = _fields(r["value"])["NULL"]
    if not null["items"]:
        return empty_real
    table = null["items"][0]
    nrow = table["dim"][0]
    columns = [table["v"][j * nrow : (j + 1) * nrow] for j in range(table["dim"][1])]
    return empty_real or any(set(column) == {"NA"} for column in columns)


def _r10(case: Json, r: Json) -> bool:
    return bool(r["files"]) and all(name.endswith(".txt") for name in r["files"])


def _r12(case: Json, r: Json) -> bool:
    # the gene with one value in every sample has a coefficient
    table = _fields(r["value"])["other.data"]
    gene = case["args"]["x"]["rownames"][-1]
    column = table["colnames"].index(f"{gene}.Estimate")
    row = table["rownames"].index(f"{case['args']['seeds']['v'][0]}.x.node")
    return bool(np.isfinite(units.dec_number(table["v"][column * table["dim"][0] + row])))


def _r13(case: Json, r: Json) -> bool:
    # the genes of the summary are not called as they are in the files
    if r["error"] is not None:
        return True
    files = case["args"]["files"]
    written = {
        line.split("\t")[0]
        for name, item in zip(files["names"], files["items"], strict=True)
        if name.startswith("real/")
        for line in item["v"][0].split("\n")[1:]
        if line
    }
    return set(_fields(_fields(r["value"])["real"])["1"]["rownames"]) != written


SHOWN_BY = {
    "R1": _r1, "R2": _r2, "R3": _r3, "R4": _r4, "R5": _r5, "R6": _r6, "R7": _r7, "R8": _r8,
    "R9": _r9, "R10": _r10, "R12": _r12, "R13": _r13,
}  # fmt: skip
WITHOUT_A_CASE = {
    "R11": "a seed with one value in every sample under a regression function: R's own answer "
    "is rounding noise, so there is none to record (the port warns: test_r2014_defects.py)",
}


def test_every_defect_of_the_package_is_accounted_for() -> None:
    assert set(SHOWN_BY) | set(WITHOUT_A_CASE) == {f"R{i}" for i in range(1, 14)}
    assert {case["defect"] for case in CASES if "defect" in case} == set(SHOWN_BY)


@pytest.mark.parametrize("defect", sorted(SHOWN_BY, key=lambda name: int(name[1:])))
def test_the_defect_is_in_r_s_answer_and_the_port_reproduces_it(defect: str) -> None:
    """The cases written for a defect show it in what R returned, and the port returns the
    same: the defect is reproduced, not only tolerated."""
    cases = [case for case in CASES if case.get("defect") == defect]
    assert cases
    for case in cases:
        assert SHOWN_BY[defect](case, ANSWERS[case["id"]]), f"{case['id']} does not show {defect}"
        assert not _differences(case), case["id"]


# ---- the cases notice a defect that is repaired -------------------------------------------------
def _getOmax_as_meant(x: RList, y: RList) -> NamedVector | None:
    """Every list against every other one, not only against those after it."""
    if len(x) == 0:
        return None
    overlap = core._overlap_block(x.items, y.items)
    best = [
        max((float(overlap[i, j]) for j in range(len(y)) if j != i), default=0.0)
        for i in range(len(x))
    ]
    return NamedVector(np.array(best), None)


def _subset_without_missing(vec: NamedVector, mask: list[bool | None]) -> NamedVector:
    """A missing condition leaves the element out, as one would expect."""
    keep = [flag is True for flag in mask]
    names = None if vec.names is None else [n for n, k in zip(vec.names, keep, strict=True) if k]
    return NamedVector(vec.values[np.array(keep, dtype=bool)], names)


def _collect_by_replicate(files: dict[str, RMatrix], features: list[str | None]) -> RMatrix:
    """Every file in a column of its own, whatever its replicate number."""
    row_of = {name: i for i, name in enumerate(features)}
    values = np.full((len(features), len(files)), NA_REAL)
    for j, key in enumerate(sorted(files, key=int)):
        table = files[key]
        for name, value in zip(table.rownames or [], table.values[:, 0], strict=True):
            values[row_of[name], j] = value
    return RMatrix(values, list(features), [str(j + 1) for j in range(len(files))])


def _always_gzip(text: str, path: str) -> None:
    with gzip.open(path + ".gz", "wb") as handle:
        handle.write(text.encode("utf-8"))


def _sort_ties_reversed(vec: NamedVector) -> NamedVector:
    keep = ~np.isnan(vec.values)
    values = vec.values[keep]
    names = None if vec.names is None else [n for n, k in zip(vec.names, keep, strict=True) if k]
    order = np.argsort(values, kind="stable")[::-1]
    return NamedVector(values[order], None if names is None else [names[i] for i in order])


REPAIRS: dict[str, tuple[Any, str, Any, str]] = {
    "R2: getOmax compares every pair": (core, "getOmax", _getOmax_as_meant, "getOmax"),
    "R7: a missing correlation is left out": (
        association, "subset_logical", _subset_without_missing, "nodes.correlation",
    ),
    "R9: every file is a column": (
        summarise, "_collect", _collect_by_replicate, "summarise.iterations",
    ),
    "R10: the files are always compressed": (
        randomise, "write_and_gzip", _always_gzip, "SEARCH.randomise",
    ),
    "R13: gene names are read as they are": (
        rio, "_row_names_as_r_reads_them", lambda tokens: tokens, "summarise.iterations",
    ),
    "write.table: NaN is written as NaN": (rio, "_cell", format_real, "write.data.files"),
    "sort: equal values in the opposite order": (
        randomise, "_sort_decreasing", _sort_ties_reversed, "write.data.files",
    ),
}  # fmt: skip


@pytest.mark.parametrize("repair", sorted(REPAIRS))
def test_the_cases_notice_a_repair(repair: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """With a defect of the package repaired in the port, or a rule of R replaced by the
    reasonable one, cases of the function fail: the cases hold the port to what R does, not to
    what would be right."""
    module, name, replacement, fn = REPAIRS[repair]
    cases = _cases_of(fn)
    assert not [case["id"] for case in cases if _differences(case)]
    monkeypatch.setattr(module, name, replacement)
    noticed = [case["id"] for case in cases if _safely_differs(case)]
    assert noticed, f"no case of {fn} notices: {repair}"


def _safely_differs(case: Json) -> bool:
    try:
        return bool(_differences(case))
    except Exception:  # a repair may also make the port fail outright
        return True


# ---- the comparison itself ----------------------------------------------------------------------
def _next_double(text: str) -> str:
    return units.enc_number(float(np.nextafter(units.dec_number(text), np.inf)))


def test_the_comparison_notices_one_bit_one_name_and_one_kind_of_missing_value() -> None:
    case = next(c for c in CASES if c["id"] == "nodes_correlation/001")
    r = ANSWERS[case["id"]]
    port = _port_answer(case)
    assert units.compare(case, r, port) == []

    def changed(edit: Any) -> list[str]:
        other = json.loads(json.dumps(r))
        edit(_fields(other["value"])["CL"]["items"][0])
        return units.compare(case, other, port)

    def one_bit(item: Json) -> None:
        item["v"][0] = _next_double(item["v"][0])

    def one_name(item: Json) -> None:
        item["names"][0] = item["names"][0] + "x"

    def missing_name(item: Json) -> None:
        item["names"][0] = None

    def not_a_number(item: Json) -> None:
        item["v"][0] = "NaN"

    def no_names(item: Json) -> None:
        item["names"] = None

    for edit in (one_bit, one_name, missing_name, not_a_number, no_names):
        assert changed(edit), edit.__name__
    assert units.compare(case, {**r, "error": "subscript out of bounds"}, port)


def test_na_and_nan_are_not_the_same_and_minus_zero_is_not_zero() -> None:
    def vector(*values: float) -> Json:
        return units.describe(NamedVector(np.array(values), None))

    assert units.differences(vector(NA_REAL), vector(NA_REAL)) == []
    assert units.differences(vector(NA_REAL), vector(float("nan")))
    assert units.differences(vector(0.0), vector(-0.0))
    assert units.differences(vector(0.1 + 0.2), vector(0.3))
    assert units.differences(vector(0.1 + 0.2), vector(0.3), tolerance=1e-9) == []


def test_numbers_reach_r_and_come_back_unchanged() -> None:
    """The cases give R their numbers as hexadecimal text and read its answers from 17
    digits. The rule ``v`` returns what it was given: what came back is what was sent."""
    for case in _cases_of("v"):
        sent = case["args"]["v"]
        back = ANSWERS[case["id"]]["value"]
        assert len(sent["v"]) == len(back["v"])
        assert units.differences(back, units.describe(units.port_value(sent))) == []
