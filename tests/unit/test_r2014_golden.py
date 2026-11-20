"""The verbatim port of SEARCH 0.0.2 (``search.r2014``) against answers recorded from R.

For this port, unlike the engine, R is the specification: the files under ``tests/data/r2014``
hold what the unmodified R package returned on a selection of cases, in R 4.4.3 on macOS arm64
and in R 3.1.0 on Linux x86-64, and the port must return the same, to the last bit for the
correlation functions and within a stated tolerance for the regression ones. The files are
written by ``validation/r_parity/py/make_r2014_golden.py`` from a run of
``validation/r_parity/run_r2014.sh``; the tests themselves need no R.
"""

from __future__ import annotations

import gzip
import json
import sys
import warnings
from pathlib import Path
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "validation" / "r_parity" / "py"))

import check_r2014_primitives  # noqa: E402
import compare_r2014  # noqa: E402
import r2014_lib as lib  # noqa: E402

from search.r2014 import platform  # noqa: E402
from search.r2014.rsem import NotReproducible  # noqa: E402

DATA = REPO / "tests" / "data" / "r2014"
PLATFORMS = {"r44_arm64": platform.R_4_ARM64, "r31_x86": platform.R_3_1_X86}
# Cases and, among them, cases in which R stops, that each file must hold. R 3.1.0 was run on
# the named cases only (it runs emulated, in a container).
MINIMUM = {"r44_arm64": (400, 50), "r31_x86": (8, 1)}


def _load(name: str) -> dict[str, Any]:
    with gzip.open(DATA / name, "rt") as handle:
        return dict(json.load(handle))


@pytest.fixture(scope="module")
def datasets() -> dict[str, Any]:
    return compare_r2014.load_datasets()


@pytest.mark.parametrize("slug", sorted(PLATFORMS))
def test_building_blocks_equal_r(slug: str) -> None:
    """Random draws, reading and writing numbers, round, sd, cor, p.adjust, rowMeans: every
    value the same as R's, bit for bit."""
    ref = _load(f"primitives_{slug}.json.gz")
    assert (
        platform.matching(ref["info"]["r_version"], int(ref["info"]["long_double_bytes"]))
        is PLATFORMS[slug]
    )
    with platform.use(PLATFORMS[slug]):
        rows = check_r2014_primitives.check(ref)
    assert len(rows) > 30
    different = [
        f"{r['what']}: {r['different']} of {r['compared']}" for r in rows if r["different"]
    ]
    assert not different, different


@pytest.mark.parametrize("slug", sorted(PLATFORMS))
def test_package_cases_equal_r(slug: str, datasets: dict[str, Any]) -> None:
    """Whole runs of ``SEARCH``, ``SEARCH.randomise`` and ``summarise.iterations``: what is
    returned, what is written, and where R stops with which message."""
    golden = _load(f"golden_{slug}.json.gz")
    failures: list[str] = []
    stops = 0
    with platform.use(PLATFORMS[slug]), warnings.catch_warnings():
        warnings.simplefilter("ignore", NotReproducible)
        for case in golden["cases"]:
            r = golden["answers"][case["id"]]
            # the correlation path is exact; the regression functions are compared to six
            # digits here, because their linear algebra differs between machines (the study
            # in docs/replication.md, on one machine, compares them to nine)
            tolerance = 1e-6 if case.get("tolerance") else 0.0
            if case["kind"] == "randomise":
                case = {**case, "sample_kind": r.get("sample_kind", case.get("sample_kind"))}
                run = lib.run_randomise(case, datasets)
                diff = lib.compare_randomise(r, run, tolerance=tolerance)
            else:
                run = lib.run_search(case, datasets)
                diff = lib.compare_search(r, run, tolerance)
            stops += r.get("error") is not None
            if run.unported is not None:
                failures.append(f"{case['id']}: not ported: {run.unported}")
            elif diff.items:
                failures.append(f"{case['id']}: {diff.items[0]}")
    n_cases, n_stops = MINIMUM[slug]
    assert len(golden["cases"]) >= n_cases
    assert stops >= n_stops, "the selection should include cases in which R stops with an error"
    assert not failures, failures[:10]
