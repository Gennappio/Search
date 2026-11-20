"""Verbatim port (``search.r2014``) against the SEARCH 0.0.2 package run by R.

    compare_r2014.py cases <dir> [--fuzz N] [--seed S] [--chunks C] [--gbm]
    compare_r2014.py check <dir> [--report report.csv] [--show K]

``cases`` writes ``<dir>/cases/<family>_<k>.json``. ``r/r2014_oracle.R`` runs one such file and
writes ``<dir>/oracle/<family>_<k>/oracle.json``. ``check`` runs the port on every case that has
an answer from R and compares the two. ``run_r2014.sh`` does the three steps.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

import r2014_cases  # noqa: E402
import r2014_lib as lib  # noqa: E402

from search.r2014 import platform  # noqa: E402
from search.r2014.rsem import NotReproducible, RMatrix  # noqa: E402

REPO = Path(__file__).resolve().parents[3]


def _number(token: str) -> float:
    return float("nan") if token in ("NA", "") else float(token)


def load_datasets() -> dict[str, RMatrix]:
    """The matrices the named cases use, read the way R reads them: every number is the nearest
    double to its decimal text (the files have at most 15 significant digits, for which R's
    reader and Python's agree)."""
    out: dict[str, RMatrix] = {}
    for name in ("Breastdata1", "Breastdata2"):
        path = REPO / "tests" / "data" / "breast" / name / "mRNA_abundance.txt.gz"
        if not path.exists():
            continue
        with gzip.open(path, "rt") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            rows: list[str] = []
            values: list[list[float]] = []
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                rows.append(parts[0].strip('"'))
                values.append([_number(p) for p in parts[1:]])
        out[name] = RMatrix(np.array(values), rows, [h.strip('"') for h in header[1:]])
    gbm = os.environ.get("SEARCH_GBM_DATA", "")
    if gbm:
        path = Path(gbm) / "DATA" / "NORMALIZED" / "IT_GBM_RNA_Deseq2_vst_Combat.csv"
        if path.exists():
            with path.open(newline="") as handle:
                reader = csv.reader(handle)
                header = next(reader)
                rows, values = [], []
                for record in reader:
                    rows.append(record[0])
                    values.append([_number(p) for p in record[1:]])
            out["IT"] = RMatrix(np.array(values), rows, header[1:])
    return out


def r_literal(value: object) -> str:
    """A case as R source: the oracle reads it with ``source()`` and needs no JSON package,
    which the R of 2014 does not have."""
    if value is None:
        return "NULL"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, int):
        return f"{value}L" if abs(value) < 2**31 else repr(float(value))
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, dict):
        return "list(" + ", ".join(f"`{k}` = {r_literal(v)}" for k, v in value.items()) + ")"
    if isinstance(value, (list, tuple)):
        return "list(" + ", ".join(r_literal(v) for v in value) + ")"
    raise TypeError(type(value))


def load_oracle(path: Path) -> dict[str, Any]:
    """The answers of R. An empty list with names was once written as ``{"":}``; files from
    that version of the oracle are still read."""
    return dict(json.loads(path.read_text().replace('{"":}', "{}")))


def write_cases(args: argparse.Namespace) -> None:
    out = args.dir / "cases"
    out.mkdir(parents=True, exist_ok=True)
    for old in [*out.glob("*.json"), *out.glob("*.R")]:
        old.unlink()
    total = 0
    for family in r2014_cases.FAMILIES:
        cases = r2014_cases.cases_of(family, args.fuzz, args.seed, args.gbm)
        if family == "named":
            # one file per case: the long ones run side by side
            chunks = [[case] for case in cases]
        else:
            n = max(1, min(args.chunks, len(cases)))
            chunks = [cases[k::n] for k in range(n)]
        for k, chunk in enumerate(chunks):
            stem = f"named_{chunk[0]['id']}" if family == "named" else f"{family}_{k + 1:03d}"
            (out / f"{stem}.json").write_text(json.dumps(chunk))
            source = "cases <- list(\n" + ",\n".join(r_literal(case) for case in chunk) + "\n)\n"
            (out / f"{stem}.R").write_text(source)
        total += len(cases)
    print(f"{total} cases written to {out}")


def check(args: argparse.Namespace) -> None:
    datasets = load_datasets()
    rows: list[dict[str, object]] = []
    status: dict[str, Counter[str]] = defaultdict(Counter)
    totals: dict[str, list[float]] = defaultdict(lambda: [0, 0, 0.0])
    shown = 0
    info: dict[str, object] = {}
    port_as = ""
    warnings.simplefilter("ignore", NotReproducible)
    for case_file in sorted((args.dir / "cases").glob("*.json")):
        oracle_file = args.dir / args.oracle / case_file.stem / "oracle.json"
        if not oracle_file.exists():
            continue
        oracle = load_oracle(oracle_file)
        info = oracle["info"]
        by_id = {r["id"]: r for r in oracle["results"]}
        chosen = platform.matching(info["r_version"], int(info["long_double_bytes"]))
        if chosen is None:
            raise SystemExit(f"the port has no platform for {info}")
        token = platform._current.set(chosen)
        port_as = chosen.name
        for case in json.loads(case_file.read_text()):
            r = by_id.get(case["id"])
            family = case["family"]
            if r is None or r.get("skipped"):
                status[family]["skipped"] += 1
                continue
            if r.get("error") and str(r["error"]).startswith("SEARCH_GBM_DATA"):
                status[family]["skipped"] += 1
                continue
            tolerance = float(case.get("tolerance", 0.0))
            if case["kind"] == "randomise":
                # the sample() R used: before R 3.6.0 there is one kind only
                case = {**case, "sample_kind": r.get("sample_kind", case.get("sample_kind"))}
                run = lib.run_randomise(case, datasets)
                diff = lib.compare_randomise(r, run, tolerance=tolerance)
            else:
                run = lib.run_search(case, datasets)
                diff = lib.compare_search(r, run, tolerance)
            r_error = r.get("error") or r.get("summarise_error") or ""
            if run.unported is not None:
                outcome = "unported"
            elif diff.items:
                outcome = "differs"
            elif r.get("error") is not None:
                outcome = "equal, R stops"
            else:
                outcome = "equal"
            status[family][outcome] += 1
            totals[family][0] += diff.compared
            totals[family][1] += diff.values
            if outcome.startswith("equal"):
                totals[family][2] = max(totals[family][2], diff.max_diff)
            rows.append(
                {
                    "case": case["id"],
                    "family": family,
                    "outcome": outcome,
                    "quantities": diff.compared,
                    "values": diff.values,
                    "largest_difference": f"{diff.max_diff:.3g}",
                    "r_error": r_error.replace("\n", " ").strip(),
                    "detail": (run.unported or "; ".join(diff.items[:3])).replace("\n", " "),
                    "note": case.get("note", ""),
                }
            )
            if (
                outcome in ("differs", "unported")
                and shown < args.show
                and "constant seed" not in family
            ):
                shown += 1
                print(f"--- {case['id']} [{outcome}] {case.get('note', '')}")
                print("   args:", {**case.get("args", {}), **case.get("assoc_args", {})})
                for item in diff.items[:6] if outcome == "differs" else [run.unported]:
                    print("   ", item)
        platform._current.reset(token)
    print(f"R side: {info}; port run as: {port_as}")
    families: list[dict[str, object]] = []
    for family in sorted(status):
        compared, values, largest = totals[family]
        counts = status[family]
        families.append(
            {
                "family": family,
                "cases": sum(counts.values()) - counts["skipped"],
                "equal": counts["equal"],
                "equal, R stops": counts["equal, R stops"],
                "differs": counts["differs"],
                "not ported": counts["unported"],
                "quantities": int(compared),
                "values": int(values),
                "largest difference among the equal": f"{largest:.3g}",
            }
        )
        print(
            f"{family}: {dict(counts)}; {int(compared)} quantities, "
            f"{int(values)} values, largest difference {largest:.3g}"
        )
    if args.report and rows:
        write_reports(args.report, rows, families, {**info, "port_run_as": port_as})


def write_reports(
    report: Path,
    rows: list[dict[str, object]],
    families: list[dict[str, object]],
    info: dict[str, object],
) -> None:
    """Three files: ``<report>.csv`` with one row per family, ``<report>_cases.csv`` with the
    named cases and every random case that is not simply equal, and ``<report>.md`` with the
    family table for reading. The row of every random case is in the output folder only."""
    report.parent.mkdir(parents=True, exist_ok=True)

    def table(path: Path, content: list[dict[str, object]]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(content[0]))
            writer.writeheader()
            writer.writerows(content)

    table(report, families)

    def random_case(r: dict[str, object]) -> bool:
        return "random matrices" in str(r["family"]) or "constant seed" in str(r["family"])

    kept = [r for r in rows if not random_case(r)]
    kept += [
        r
        for r in rows
        if "random matrices" in str(r["family"]) and r["outcome"] in ("differs", "unported")
    ]
    table(report.with_name(report.stem + "_cases.csv"), kept)
    lines = [
        f"# Verbatim port against {info['r_version']} ({info['platform']})",
        "",
        f"Port run as: {info['port_run_as']}. R side: `sample()` kind {info['sample_kind']}, "
        f"long double of {info['long_double_bytes']} bytes, "
        f"collation {info.get('locale_collate', '?')}, run on {info['date']}.",
        "",
        "| " + " | ".join(families[0]) + " |",
        "|" + "---|" * len(families[0]),
        *("| " + " | ".join(str(v) for v in f.values()) + " |" for f in families),
        "",
    ]
    report.with_suffix(".md").write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    p_cases = sub.add_parser("cases")
    p_cases.add_argument("dir", type=Path)
    p_cases.add_argument("--fuzz", type=int, default=2000)
    p_cases.add_argument("--seed", type=int, default=20141)
    p_cases.add_argument("--chunks", type=int, default=8)
    p_cases.add_argument("--gbm", action="store_true")
    p_check = sub.add_parser("check")
    p_check.add_argument("dir", type=Path)
    p_check.add_argument("--oracle", default="oracle", help="folder of the R answers under <dir>")
    p_check.add_argument("--report", type=Path)
    p_check.add_argument("--show", type=int, default=5)
    args = parser.parse_args()
    if args.command == "cases":
        write_cases(args)
    else:
        check(args)


if __name__ == "__main__":
    main()
