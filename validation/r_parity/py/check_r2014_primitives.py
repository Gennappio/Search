"""The building blocks of the verbatim port against R: random draws, numbers to text and back,
``round``, ``sd``, ``cor``, ``p.adjust``, ``rowMeans``.

    check_r2014_primitives.py <primitives.json> [--report report.csv]

``primitives.json`` is written by ``r/r2014_primitives.R``. Every comparison is exact: the same
integers, the same characters, the same double bit for bit.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from search.r2014 import platform, rformat, rstats
from search.r2014.rng import RRandom
from search.r2014.rsem import NA_REAL, RError, is_na_real
from search.r2014.summarise import _row_means


def numbers(values: list[str]) -> np.ndarray:
    return np.array([NA_REAL if v == "NA" else float(v) for v in values], dtype=np.float64)


def differing(a: np.ndarray, b: np.ndarray) -> int:
    a = np.ascontiguousarray(a, dtype=np.float64).ravel()
    b = np.ascontiguousarray(b, dtype=np.float64).ravel()
    if a.shape != b.shape:
        return max(len(a), len(b))
    both_missing = np.isnan(a) & np.isnan(b) & (is_na_real(a) == is_na_real(b))
    return int((~((a.view(np.int64) == b.view(np.int64)) | both_missing)).sum())


def check(ref: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    def row(block: str, what: str, n: int, bad: int, note: str = "") -> None:
        rows.append({"block": block, "what": what, "compared": n, "different": bad, "note": note})

    # sample()
    for kind in sorted({c["kind"] for c in ref["rng"]}):
        n = bad = 0
        for c in (c for c in ref["rng"] if c["kind"] == kind):
            seed, size = int(c["seed"]), int(c["n"])
            boot = RRandom(seed, kind).sample(size, replace=True).astype(float)
            jack = (
                RRandom(seed, kind).sample(size, size - 1, replace=False).astype(float)
                if size > 1
                else np.empty(0)
            )
            gen = RRandom(seed, kind)
            unif = np.array([gen.unif_rand() for _ in range(3)])
            ok = (
                differing(boot[:200], numbers(c["boot_head"])) == 0
                and float((boot * np.arange(1, len(boot) + 1)).sum()) == float(c["boot_sum"])
                and differing(jack[:200], numbers(c["jack_head"])) == 0
                and float((jack * np.arange(1, len(jack) + 1)).sum()) == float(c["jack_sum"])
                and differing(unif, numbers(c["unif"])) == 0
            )
            n += 1
            bad += not ok
        row("random numbers", f"set.seed and sample(), kind {kind}: seeds x sizes", n, bad)

    # numbers to text
    values = numbers(ref["format"]["value"])
    text = ref["format"]["text"]
    bad = sum(rformat.format_real(float(v)) != t for v, t in zip(values, text, strict=True))
    row("numbers to text", "write.table of a number (15 significant digits)", len(text), bad)

    # text to numbers, by number of significant digits
    texts = ref["read"]["text"]
    read = numbers(ref["read"]["value"])
    by_digits: dict[int, list[int]] = {}
    for t, v in zip(texts, read, strict=True):
        digits = len(t.lstrip("-").split("e")[0].replace(".", "").lstrip("0"))
        tally = by_digits.setdefault(digits if digits <= 15 else min(digits, 20), [0, 0, 0])
        tally[0] += 1
        tally[1] += rformat.r_strtod(t) != v
        tally[2] += float(t) != v
    for digits in sorted(by_digits):
        n, bad, not_nearest = by_digits[digits]
        row(
            "text to numbers",
            f"as.numeric of a text with {digits} significant digits",
            n,
            bad,
            f"R read {not_nearest} of them into a double that is not the nearest",
        )

    # round()
    values = numbers(ref["round"]["value"])
    row(
        "round",
        "round(x, 5)",
        len(values),
        differing(rstats.r_round(values, 5), numbers(ref["round"]["five"])),
    )
    row(
        "round",
        "round(x, 3)",
        len(values),
        differing(rstats.r_round(values, 3), numbers(ref["round"]["three"])),
    )

    # sd() and cor()
    tally_cor: dict[tuple[str, str], list[int]] = {}
    sd_n = sd_bad = 0
    for c in ref["cor"]:
        g, n = int(c["g"]), int(c["n"])
        x = numbers(c["x"]).reshape(n, g).T.copy()
        seeds = [int(float(s)) for s in c["seeds"]]
        tally = tally_cor.setdefault((c["method"], c["use"]), [0, 0])
        try:
            got: np.ndarray | None = rstats.r_cor(x[seeds], x, c["method"], c["use"])
            error = None
        except RError as exc:
            got, error = None, exc.message
        if c.get("error") is not None or error is not None:
            tally[0] += 1
            tally[1] += c.get("error") != error
        else:
            assert got is not None
            expected = numbers(c["r"])
            tally[0] += len(expected)
            tally[1] += differing(got.T, expected)
        expected_sd = numbers(c["sd"])
        sd_n += len(expected_sd)
        sd_bad += differing(rstats.row_sd(x), expected_sd)
    for (method, use), (n, bad) in sorted(tally_cor.items()):
        row("cor", f"cor(method = {method}, use = {use})", n, bad)
    row("sd", "sd of a row", sd_n, sd_bad)

    # p.adjust()
    tally_p: dict[str, list[int]] = {}
    for c in ref["padj"]:
        expected = numbers(c["out"])
        tally = tally_p.setdefault(c["method"], [0, 0])
        tally[0] += len(expected)
        tally[1] += differing(rstats.p_adjust(numbers(c["p"]), c["method"]), expected)
    for method, (n, bad) in sorted(tally_p.items()):
        row("p.adjust", f"p.adjust(method = {method})", n, bad)

    # rowMeans()
    rm = ref["row_means"]
    matrix = numbers(rm["value"]).reshape(int(rm["ncol"]), int(rm["nrow"])).T
    expected = numbers(rm["mean"])
    got_mean = _row_means(matrix)
    row(
        "rowMeans",
        "rowMeans(x, na.rm = TRUE) of ten values with five decimals",
        len(expected),
        differing(got_mean, expected),
    )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("primitives", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    ref = json.loads(args.primitives.read_text())
    info = ref["info"]
    chosen = platform.matching(info["r_version"], int(info["long_double_bytes"]))
    if chosen is None:
        print("R side:", info)
        raise SystemExit("the port has no platform for this R")
    with platform.use(chosen):
        rows = check(ref)
    print("R side:", info, "| port as:", chosen.name)
    for r in rows:
        flag = "" if r["different"] == 0 else "   <-- differs"
        counts = f"{r['compared']:7d} compared, {r['different']:6d} different"
        print(f"{r['block']:16s} {r['what']:62s} {counts}{flag}")
    for c in ref["constant_seed"]:
        print("constant seed", c)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with args.report.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
