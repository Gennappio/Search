"""R against R: the same cases answered by two installations of R, compared with each other.

    compare_r_versions.py <out/r2014> [--a oracle] [--b oracle_r31] [--report report.csv]

No port is involved. This measures how much a result of the SEARCH package depends on the R
that computes it. A case is *identical* when everything R returned and wrote is the same in
both; it differs *in numbers* when the same genes are in the same places with other values;
it differs *in results* when a gene enters or leaves a neighbourhood, a seed list or the
summary, or one R stops with an error and the other does not. Runs of ``SEARCH.randomise`` are
compared only where both drew the same samples (R changed ``sample()`` in version 3.6.0).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).parent))

import compare_r2014  # noqa: E402

Json = dict[str, Any]


def _names_only(value: Any) -> Any:
    """A result with every number taken out: what is left is which genes are where."""
    if isinstance(value, dict):
        return {k: _names_only(v) for k, v in value.items() if k not in ("values", "warnings")}
    if isinstance(value, list):
        return [_names_only(v) for v in value]
    return value


def _table_names(text: str) -> list[str]:
    lines = text.split("\n")
    return [lines[0], *(line.split("\t", 1)[0] for line in lines[1:])]


def classify(case: Json, a: Json, b: Json) -> tuple[str, str]:
    if case["kind"] == "randomise" and a.get("sample_kind") != b.get("sample_kind"):
        return "other random draws", ""
    drop = ("id", "warnings", "sample_kind")
    ra = {k: v for k, v in a.items() if k not in drop}
    rb = {k: v for k, v in b.items() if k not in drop}
    if ra == rb:
        return "identical", ""
    if (a.get("error"), a.get("summarise_error")) != (b.get("error"), b.get("summarise_error")):
        stop_a = a.get("error") or a.get("summarise_error")
        stop_b = b.get("error") or b.get("summarise_error")
        return "results", f"errors: {stop_a!r} / {stop_b!r}"
    if case["kind"] != "randomise":
        if _names_only(ra) != _names_only(rb):
            return "results", "different genes in D, in a neighbourhood or among the seeds"
        return "numbers", ""
    fa, fb = a.get("files") or {}, b.get("files") or {}
    if sorted(fa) != sorted(fb) or a.get("summaries") != b.get("summaries"):
        return "results", "different files or sizes"
    for name in fa:
        if fa[name] != fb[name] and "connectivity" not in name and "correlation" not in name:
            if _table_names(fa[name]) != _table_names(fb[name]) or "network_growth" in name:
                return "results", f"{name.rsplit('/', 1)[-1]} lists other genes"
    sa, sb = a.get("summarise"), b.get("summarise")
    if sa is not None and sb is not None:
        for it in sa["real"]:
            ma, mb = sa["real"][it], sb["real"].get(it)
            if mb is None or sorted(ma["rownames"] or []) != sorted(mb["rownames"] or []):
                return "results", "the summary lists other genes"
            rows_a = dict(
                zip(ma["rownames"] or [], ma["values"][-len(ma["rownames"] or []) :], strict=False)
            )
            rows_b = dict(
                zip(mb["rownames"] or [], mb["values"][-len(mb["rownames"] or []) :], strict=False)
            )
            if rows_a != rows_b:
                changed = sum(rows_a[g] != rows_b[g] for g in rows_a)
                return "results", f"the p-value of {changed} of {len(rows_a)} genes"
    return "numbers", ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=Path)
    parser.add_argument("--a", default="oracle")
    parser.add_argument("--b", default="oracle_r31")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    reasons: dict[str, Counter[str]] = defaultdict(Counter)
    rows: list[Json] = []
    info: dict[str, Any] = {}
    for case_file in sorted((args.dir / "cases").glob("*.json")):
        file_a = args.dir / args.a / case_file.stem / "oracle.json"
        file_b = args.dir / args.b / case_file.stem / "oracle.json"
        if not (file_a.exists() and file_b.exists()):
            continue
        oracle_a, oracle_b = compare_r2014.load_oracle(file_a), compare_r2014.load_oracle(file_b)
        info = {"a": oracle_a["info"], "b": oracle_b["info"]}
        by_a = {r["id"]: r for r in oracle_a["results"]}
        by_b = {r["id"]: r for r in oracle_b["results"]}
        for case in json.loads(case_file.read_text()):
            a, b = by_a.get(case["id"]), by_b.get(case["id"])
            if a is None or b is None or a.get("skipped") or b.get("skipped"):
                continue
            outcome, reason = classify(case, a, b)
            counts[case["family"]][outcome] += 1
            if reason:
                reasons[case["family"]][
                    reason.split(":")[0] if reason.startswith("errors") else reason.split(" of ")[0]
                ] += 1
            rows.append(
                {"case": case["id"], "family": case["family"], "outcome": outcome, "detail": reason}
            )
    print(json.dumps(info))
    for family in sorted(counts):
        print(f"{family}: {dict(counts[family])}")
        for reason, n in reasons[family].most_common(6):
            print(f"      {n:6d}  {reason}")
    if args.report and rows:
        with args.report.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
