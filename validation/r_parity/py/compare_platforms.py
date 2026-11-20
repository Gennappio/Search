"""The port against itself: the same cases computed as two installations of R would.

    compare_platforms.py <out/r2014> [--limit N]

``search.r2014`` reproduces R 4.4.3 on arm64 and R 3.1.0 on x86-64. Here it runs every random
case of the correlation path as each of the two, with the same random draws, and counts the
cases in which the two results differ: how much a result of the package depends on the
processor and on the version of R alone. This is a prediction made with the port, checked
against R only as far as the port itself is; ``compare_r_versions.py`` is the measurement in R,
on the cases both installations ran.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

import compare_r2014  # noqa: E402
import r2014_lib as lib  # noqa: E402

from search.r2014 import platform  # noqa: E402
from search.r2014.core import SearchResult  # noqa: E402
from search.r2014.rsem import NamedVector, RMatrix  # noqa: E402


def snapshot(value: Any, numbers: bool) -> Any:
    """A result as plain data; without ``numbers``, only which genes are where."""
    if isinstance(value, NamedVector):
        return (value.names, value.values.tobytes() if numbers else len(value))
    if isinstance(value, RMatrix):
        return (
            value.rownames,
            value.colnames,
            value.values.tobytes() if numbers else value.values.shape,
        )
    if isinstance(value, SearchResult):
        return [
            value.removed,
            value.correlationLists.names,
            [snapshot(item, numbers) for item in value.correlationLists.items],
            {k: snapshot(v, numbers) for k, v in value.D.items()},
            {k: snapshot(v, numbers) for k, v in value.O.items()},
            value.network_growth,
        ]
    if isinstance(value, np.ndarray):
        return value.tobytes() if numbers else value.shape
    return value


def run(
    case: dict[str, Any], datasets: dict[str, RMatrix], which: platform.RPlatform
) -> lib.PortRun:
    with platform.use(which):
        if case["kind"] == "randomise":
            return lib.run_randomise({**case, "sample_kind": "Rounding"}, datasets)
        return lib.run_search(case, datasets)


def classify(case: dict[str, Any], a: lib.PortRun, b: lib.PortRun) -> str:
    if (a.error, a.unported) != (b.error, b.unported):
        return "one stops and the other does not"
    if case["kind"] != "randomise":
        if a.error is not None or snapshot(a.result, True) == snapshot(b.result, True):
            return "identical"
        return (
            "numbers only"
            if snapshot(a.result, False) == snapshot(b.result, False)
            else "other genes"
        )
    if a.extra.get("summarise_error") != b.extra.get("summarise_error"):
        return "one stops and the other does not"
    fa, fb = a.extra.get("files") or {}, b.extra.get("files") or {}
    sa, sb = a.extra.get("summarise"), b.extra.get("summarise")
    p_a = None if sa is None else {k: (m.rownames, m.values.tobytes()) for k, m in sa.real.items()}
    p_b = None if sb is None else {k: (m.rownames, m.values.tobytes()) for k, m in sb.real.items()}
    if fa == fb and p_a == p_b:
        return "identical"
    if fa != fb:
        return "the files differ"
    return "the files are the same, the summary (mean, p-value) differs"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=Path)
    parser.add_argument("--limit", type=int, default=0, help="cases per file, 0 for all")
    args = parser.parse_args()
    datasets = compare_r2014.load_datasets()
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    warnings.simplefilter("ignore")
    for case_file in sorted((args.dir / "cases").glob("*.json")):
        if not case_file.name.startswith(("search_", "randomise_")):
            continue
        cases = json.loads(case_file.read_text())
        for case in cases[: args.limit] if args.limit else cases:
            a = run(case, datasets, platform.R_4_ARM64)
            b = run(case, datasets, platform.R_3_1_X86)
            counts[case["family"]][classify(case, a, b)] += 1
    for family in sorted(counts):
        total = sum(counts[family].values())
        print(f"{family}: {total} cases")
        for outcome, n in counts[family].most_common():
            print(f"   {n:6d}  {outcome}")


if __name__ == "__main__":
    main()
