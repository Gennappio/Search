"""Golden answers for the test suite: a selection of the cases of ``run_r2014.sh`` with what R
returned for them, so that the tests of the verbatim port need no R.

    make_r2014_golden.py <out/r2014> <tests/data/r2014> [--primitives name=primitives.json ...]

For every folder of R answers under ``<out/r2014>`` (``oracle`` for the R of the machine,
``oracle_r31`` for R 3.1.0) it writes ``golden_<platform>.json.gz``: the first cases of each
random family and a few named ones on Breastdata1. Nothing derived from the GBM cohorts is
included. Each ``--primitives`` file is cut down to ``primitives_<name>.json.gz``.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).parent))

import compare_r2014  # noqa: E402

RANDOM = {"search": 300, "randomise": 60, "regression": 60}
NAMED = (
    "breast_single_pass",
    "breast_sd_prefilter",
    "breast_universe_reversed",
    "breast_doc_not_immune",
    "breast_missing_values",
    "breast_one_seed",
    "breast_bootstrap_defaults",
    "breast_manual_correlation",
    "breast_blank_in_path",
    "breast_bootstrap_6iter",
    "names_whole_numbers",
    "names_decimal_numbers",
    "names_colliding_numbers",
    "names_logical",
    "names_all_logical",
    "names_with_na",
)
SLUG = {"oracle": "r44_arm64", "oracle_r31": "r31_x86"}


def write(path: Path, content: Any) -> None:
    with gzip.GzipFile(path, "wb", compresslevel=9, mtime=0) as handle:
        handle.write(json.dumps(content, separators=(",", ":")).encode())
    print(f"{path}: {path.stat().st_size / 1e6:.2f} MB")


def golden(run_dir: Path, oracle: str) -> dict[str, Any] | None:
    cases: list[Any] = []
    answers: dict[str, Any] = {}
    info: dict[str, Any] = {}
    for family, limit in RANDOM.items():
        pool: dict[str, tuple[Any, Any]] = {}
        for case_file in sorted((run_dir / "cases").glob(f"{family}_*.json")):
            oracle_file = run_dir / oracle / case_file.stem / "oracle.json"
            if not oracle_file.exists():
                continue
            answer = compare_r2014.load_oracle(oracle_file)
            info = answer["info"]
            by_id = {r["id"]: r for r in answer["results"]}
            for case in json.loads(case_file.read_text()):
                r = by_id.get(case["id"])
                if r is not None and not r.get("skipped") and "constant seed" not in case["family"]:
                    pool[case["id"]] = (case, r)
        for case_id in sorted(pool)[:limit]:
            cases.append(pool[case_id][0])
            answers[case_id] = pool[case_id][1]
    for name in NAMED:
        case_file = run_dir / "cases" / f"named_{name}.json"
        oracle_file = run_dir / oracle / f"named_{name}" / "oracle.json"
        if case_file.exists() and oracle_file.exists():
            answer = compare_r2014.load_oracle(oracle_file)
            cases.append(json.loads(case_file.read_text())[0])
            answers[name] = answer["results"][0]
    if not cases:
        return None
    return {"info": info, "cases": cases, "answers": answers}


def cut_primitives(ref: dict[str, Any]) -> dict[str, Any]:
    """A sample of every block, small enough for the repository."""

    def every(values: list[Any], step: int) -> list[Any]:
        return values[::step]

    out = dict(ref)
    out["rng"] = [c for c in ref["rng"] if int(c["n"]) <= 1000]
    out["format"] = {k: every(v, 20) for k, v in ref["format"].items()}
    out["read"] = {k: every(v, 7) for k, v in ref["read"].items()}
    out["round"] = {k: every(v, 25) for k, v in ref["round"].items()}
    out["cor"] = ref["cor"][: 15 * 12]
    out["padj"] = ref["padj"][: 8 * 40]
    rm = ref["row_means"]
    nrow, ncol, keep = int(rm["nrow"]), int(rm["ncol"]), 400
    rows = [i for i in range(0, nrow, nrow // keep)][:keep]
    out["row_means"] = {
        "nrow": len(rows),
        "ncol": ncol,
        "value": [rm["value"][j * nrow + i] for j in range(ncol) for i in rows],
        "mean": [rm["mean"][i] for i in rows],
    }
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--primitives", action="append", default=[], help="name=path")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for oracle, slug in SLUG.items():
        content = golden(args.run_dir, oracle)
        if content is not None:
            write(args.out_dir / f"golden_{slug}.json.gz", content)
    for item in args.primitives:
        name, _, path = item.partition("=")
        write(
            args.out_dir / f"primitives_{name}.json.gz",
            cut_primitives(json.loads(Path(path).read_text())),
        )


if __name__ == "__main__":
    main()
