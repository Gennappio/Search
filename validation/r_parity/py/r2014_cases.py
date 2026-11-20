"""Cases for the comparison of the verbatim port (``search.r2014``) with the R package.

Random cases (``fuzz_*``) are small matrices with random arguments, built to reach the states
in which the package misbehaves: with a handful of samples many correlations sit exactly on the
threshold, missing values and constant rows are frequent, and the universe is given in an order
that is not the one of the matrix. Named cases run the package on Breastdata1, the data it ships
with, the way its manual does, and on the GBM IT cohort when it is available.

Every case belongs to a ``family``; the report counts by family:

- ``SEARCH()``, ``SEARCH.randomise()``: correlation; compared bit for bit, files included.
- ``regression``: ``nodes.glm`` and ``nodes.fastLm``; names, memberships, D and O compared
  exactly, coefficients and p-values with ``REGRESSION_TOLERANCE``.
- ``regression, constant seed``: a seed with the same value in every sample. R's own answer is
  rounding noise there (it changes when the samples are listed in another order), so these
  cases are run and counted but nothing can be expected of them.

Everything is derived from a fixed random seed, so a case keeps its number between runs.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np

Json = dict[str, Any]

REPO = Path(__file__).resolve().parents[3]

REGRESSION_TOLERANCE = 1e-9
"""Coefficients and p-values of ``nodes.glm`` and ``nodes.fastLm`` come out of different
linear-algebra and distribution routines in the two languages; they are compared relative to
their size with this tolerance."""


def _fmt(v: float) -> str | None:
    return None if np.isnan(v) else repr(float(v))


def _random_matrix(
    rng: np.random.Generator,
    min_samples: int = 3,
    max_samples: int = 10,
    styles: tuple[str, ...] = ("ties", "continuous", "blocks", "blocks_ties"),
    p_constant: float = 0.25,
    p_missing: float = 0.12,
) -> tuple[np.ndarray, list[str], list[str]]:
    g = int(rng.integers(4, 15))
    n = int(rng.integers(min_samples, max_samples + 1))
    style = str(rng.choice(styles))
    if style == "ties":
        x = rng.integers(0, 5, size=(g, n)).astype(float)
    elif style == "continuous":
        x = np.round(rng.normal(8, 1.5, size=(g, n)), 2)
    else:
        n_blocks = int(rng.integers(1, 4))
        signals = rng.normal(0, 1, size=(n_blocks, n))
        member = rng.integers(-1, n_blocks, size=g)  # -1: noise only
        sign = rng.choice([1.0, 1.0, -1.0], size=g)
        noise = rng.normal(0, rng.choice([0.2, 0.6, 1.2]), size=(g, n))
        x = (
            np.where(member[:, None] >= 0, sign[:, None] * signals[np.maximum(member, 0)], 0.0)
            + noise
            + 8
        )
        x = np.round(x, 0 if style == "blocks_ties" else 2)
    if rng.random() < p_constant:  # constant rows
        for r in rng.choice(g, size=int(rng.integers(1, 3)), replace=False):
            x[r] = np.round(rng.normal(8, 1), 1)
    if rng.random() < 0.15:  # low-variance rows, for the sd prefilter
        for r in rng.choice(g, size=int(rng.integers(1, 4)), replace=False):
            x[r] = 8 + np.round(rng.normal(0, 0.05, size=n), 2)
    if rng.random() < p_missing:  # missing values
        for _ in range(int(rng.integers(1, 4))):
            x[int(rng.integers(g)), int(rng.integers(n))] = np.nan
    order = rng.permutation(g)
    rownames = [f"g{i + 1}" for i in order]
    colnames = [f"s{j + 1}" for j in range(n)]
    return x, rownames, colnames


def _data(x: np.ndarray, rownames: list[str], colnames: list[str]) -> Json:
    return {
        "nrow": len(rownames),
        "ncol": len(colnames),
        "rownames": rownames,
        "colnames": colnames,
        "values": [_fmt(v) for v in x.ravel()],
    }


def _universe(rng: np.random.Generator, rownames: list[str], exotic: bool = True) -> list[str]:
    g = len(rownames)
    u = rng.random()
    if u < 0.45:
        return list(rownames)  # matrix order
    if u < 0.78:
        return [rownames[i] for i in rng.permutation(g)]  # another order
    if u < 0.95 or not exotic:
        keep = rng.choice(g, size=int(rng.integers(2, g + 1)), replace=False)
        return [rownames[i] for i in keep]  # a subset; the seeds are added by SEARCH
    return [*rownames, "absent1"] if rng.random() < 0.5 else [rownames[0]]


# ---- SEARCH() with correlations -----------------------------------------------------------------
def fuzz_search_cases(n_cases: int, seed: int = 20141, prefix: str = "search") -> list[Json]:
    rng = np.random.default_rng(seed)
    cases: list[Json] = []
    for k in range(n_cases):
        x, rownames, colnames = _random_matrix(rng)
        g = len(rownames)
        n_seeds = int(rng.choice([1, 2, 2, 3, 3, 4, 5]))
        seeds = [rownames[i] for i in rng.choice(g, size=min(n_seeds, g), replace=False)]
        universe = _universe(rng, rownames)
        if rng.random() < 0.02:
            seeds = [*seeds, "not_in_matrix"]
        method = str(rng.choice(["spearman"] * 6 + ["pearson"] * 2 + ["kendall"]))
        assoc: Json = {
            "correlation.threshold": float(rng.choice([0.5, 0.5, 0.5, 0.3, 0.7, 0.0, 1.0])),
            "correlation.method": method,
        }
        if rng.random() < 0.12:
            assoc["correlation.use"] = str(
                rng.choice(["pairwise.complete.obs", "complete.obs", "all.obs", "na.or.complete"])
            )
        args: Json = {
            "scale": bool(rng.random() < 0.85),
            "D.lower": float(rng.choice([0.0, 0.0, 0.1, 0.1, 0.3, 0.5])),
            "D.higher": float(rng.choice([0.0, 0.1, 0.24, 0.24, 0.5, 0.8])),
            "O.lower": float(rng.choice([0.0, 0.01, 0.01, 0.1, 0.3, 0.6])),
            "keep.Alist": bool(rng.random() < 0.6),
            "max.iter": int(rng.choice([1, 1, 2, 3, 4, 6])),
            "sd.threshold": float(rng.choice([0.0, 0.0, 0.0, 0.1, 0.5, 1.0])),
            "verbose": False,
        }
        if not args["scale"]:
            args["D.lower"] = float(rng.choice([0.0, 1.0, 2.0]))
            args["D.higher"] = float(rng.choice([0.0, 1.0, 2.0]))
        cases.append(
            {
                "id": f"{prefix}_{k + 1:05d}",
                "kind": "search",
                "family": "SEARCH(), random matrices",
                "data": _data(x, rownames, colnames),
                "seeds": seeds,
                "universe": universe,
                "args": args,
                "assoc_args": assoc,
            }
        )
    return cases


# ---- SEARCH.randomise() and summarise.iterations() ----------------------------------------------
def fuzz_randomise_cases(n_cases: int, seed: int = 20142, prefix: str = "randomise") -> list[Json]:
    rng = np.random.default_rng(seed)
    cases: list[Json] = []
    for k in range(n_cases):
        x, rownames, colnames = _random_matrix(rng, min_samples=5, max_samples=12, p_missing=0.06)
        g = len(rownames)
        seeds = [
            rownames[i]
            for i in rng.choice(g, size=min(int(rng.choice([2, 3, 3, 4])), g), replace=False)
        ]
        universe = _universe(rng, rownames, exotic=False)
        max_iter = int(rng.choice([1, 1, 2, 3, 4]))
        n_boot = int(rng.choice([1, 2, 3, 5]))
        args: Json = {
            "scale": bool(rng.random() < 0.9),
            "D.lower": float(rng.choice([0.0, 0.0, 0.1, 0.3])),
            "D.higher": float(rng.choice([0.0, 0.1, 0.24, 0.24, 0.5])),
            "O.lower": float(rng.choice([0.0, 0.01, 0.01, 0.1, 0.3])),
            "keep.seeds": bool(rng.random() < 0.8),
            "max.iter": max_iter,
            "sd.threshold": float(rng.choice([0.0, 0.0, 0.0, 0.1, 0.5])),
            "verbose": False,
            "n.boot": n_boot,
            "sampling": str(rng.choice(["bootstrap", "bootstrap", "jacknife"])),
            "sampling.seed": int(rng.integers(1, 10000)),
        }
        if rng.random() < 0.1:
            args["n.boot.start"] = int(rng.integers(2, 5))
        odd = rng.random()
        if odd < 0.01:
            args["n.boot"] = 0  # 1:0 in R counts down: replicates 1 and 0 are run
        elif odd < 0.02:
            args["sampling"] = "permutation"  # not a sampling the package knows
        elif odd < 0.03:
            seeds = ["not_in_matrix"]
        assoc: Json = {
            "correlation.threshold": float(rng.choice([0.5, 0.5, 0.3, 0.7])),
            "correlation.method": str(rng.choice(["spearman", "spearman", "spearman", "pearson"])),
        }
        batches = int(rng.choice([1, 2, 2, 3]))
        n_boot_null = int(rng.choice([1, 2, 3]))
        null_runs = []
        candidates = [r for r in rownames if r not in seeds]
        for b in range(batches):
            if len(candidates) < len(seeds) or rng.random() < 0.05:
                continue  # a batch that was never run
            null_seeds = [
                candidates[i] for i in rng.choice(len(candidates), size=len(seeds), replace=False)
            ]
            null_runs.append(
                {
                    "dir": f"null/batch_{b + 1}_0_FALSE_0_0/",
                    "seeds": null_seeds,
                    "n.boot": n_boot_null,
                }
            )
        # most often the summary is asked for the first iteration, which every replicate has
        summarise_iter = (
            1 if rng.random() < 0.6 else (max_iter if rng.random() < 0.9 else max_iter + 1)
        )
        case: Json = {
            "id": f"{prefix}_{k + 1:05d}",
            "kind": "randomise",
            "family": "SEARCH.randomise(), random matrices",
            "data": _data(x, rownames, colnames),
            "seeds": seeds,
            "universe": universe,
            "args": args,
            "assoc_args": assoc,
            "null_runs": null_runs,
            "summarise": {
                "max.iter": summarise_iter,
                "n.boot": n_boot
                if "n.boot.start" not in args
                else n_boot + int(args["n.boot.start"]),
                "n.boot.NULL": n_boot_null,
                "batches.NULL": batches,
            },
        }
        if rng.random() < 0.25:
            case["sample_kind"] = "Rounding"
        cases.append(case)
    return cases


# ---- regression functions -----------------------------------------------------------------------
def fuzz_glm_cases(n_cases: int, seed: int = 20143, prefix: str = "regression") -> list[Json]:
    rng = np.random.default_rng(seed)
    cases: list[Json] = []
    for k in range(n_cases):
        # enough samples, and no small whole numbers: with few of either a model can fit the
        # seed exactly, and what both languages return is then rounding noise
        x, rownames, colnames = _random_matrix(
            rng, min_samples=12, max_samples=18, styles=("continuous", "blocks"), p_missing=0.0
        )
        g, n = x.shape
        function = str(rng.choice(["nodes.glm", "nodes.glm", "nodes.fastLm"]))
        family = "gaussian"
        if function == "nodes.glm":
            family = str(
                rng.choice(
                    [
                        "gaussian",
                        "gaussian",
                        "Gamma",
                        "inverse.gaussian",
                        "poisson",
                        "quasi",
                        "quasipoisson",
                    ]
                )
            )
        if rng.random() < 0.02:
            x = x - 9.0  # negative values: the families that need positive ones stop
        if rng.random() < 0.12:
            x[int(rng.integers(g)), int(rng.integers(n))] = np.nan
        seeds = [
            rownames[i]
            for i in rng.choice(g, size=min(int(rng.choice([1, 2, 3, 3, 4])), g), replace=False)
        ]
        max_iter = int(rng.choice([1, 1, 2, 3]))
        constant = {rownames[i] for i in range(g) if np.nanmax(x[i]) == np.nanmin(x[i])}
        # A constant seed is fitted exactly and R's answer for it is noise (R defect R11). With
        # nodes.fastLm a constant gene can enter a neighbourhood and become a seed (R12).
        noise = any(s in constant for s in seeds) or (
            function == "nodes.fastLm" and bool(constant) and max_iter > 1
        )
        assoc: Json = {
            "significance.threshold": float(rng.choice([0.01, 0.05, 0.05, 0.2, 0.5, 1e-5])),
            "glm.family": family,
        }
        if rng.random() < 0.3:
            assoc["adjust.p"] = str(
                rng.choice(["BH", "bonferroni", "holm", "BY", "hochberg", "hommel", "fdr", "none"])
            )
        case: Json = {
            "id": f"{prefix}_{k + 1:05d}",
            "kind": "search",
            "family": "regression, constant seed" if noise else "regression, random matrices",
            "tolerance": REGRESSION_TOLERANCE,
            "data": _data(x, rownames, colnames),
            "seeds": seeds,
            "universe": _universe(rng, rownames, exotic=False),
            "args": {
                "scale": bool(rng.random() < 0.9),
                "D.lower": float(rng.choice([0.0, 0.0, 0.1, 0.3])),
                "D.higher": float(rng.choice([0.0, 0.1, 0.24, 0.5])),
                "O.lower": float(rng.choice([0.0, 0.01, 0.1, 0.3])),
                "keep.Alist": bool(rng.random() < 0.7),
                "max.iter": max_iter,
                "sd.threshold": 0.0,
                "verbose": False,
                "nodes.association.function": function,
            },
            "assoc_args": assoc,
        }
        if rng.random() < 0.45:
            batch = [str(rng.choice(["a", "b", "c"])) for _ in range(n)]
            age = [float(np.round(rng.normal(60, 10), 1)) for _ in range(n)]
            columns: list[Json] = [
                {"name": "age", "type": "numeric", "values": age},
                {"name": "batch", "type": "character", "values": batch},
            ]
            if rng.random() < 0.15:
                columns[0]["values"][int(rng.integers(n))] = None
            case["covariates"] = {"columns": columns}
            if rng.random() < 0.8:  # without formula.string the covariates are ignored
                assoc["formula.string"] = str(
                    rng.choice(
                        [
                            "age",
                            "batch",
                            "age + batch",
                            "age + batch",
                            "x.node:batch",
                            "age*batch",
                            "batch + x.node:age",
                        ]
                    )
                )
        cases.append(case)
    return cases


# ---- the package on its own data ----------------------------------------------------------------
def _breast_names(name: str = "Breastdata1") -> tuple[list[str], list[str]]:
    import gzip

    path = REPO / "tests" / "data" / "breast" / name / "mRNA_abundance.txt.gz"
    with gzip.open(path, "rt") as handle:
        header = [h.strip('"') for h in handle.readline().rstrip("\n").split("\t")[1:]]
        rows = [line.split("\t", 1)[0].strip('"') for line in handle]
    return rows, header


def _breast_covariates(
    samples: list[str], columns: dict[str, str], name: str = "Breastdata1"
) -> Json:
    """Columns of the patient annotation, in the order of the samples of the matrix; R matches
    the two through ``make.names`` of the patient identifier."""
    path = REPO / "tests" / "data" / "breast" / name / "patient_annotation.txt"
    with path.open() as handle:
        table = {row["ID"].replace("-", "."): row for row in csv.DictReader(handle, delimiter="\t")}
    return {
        "columns": [
            {"name": new, "type": "numeric", "values": [float(table[s][old]) for s in samples]}
            for old, new in columns.items()
        ]
    }


def _seed_file(path: str) -> list[str]:
    with (REPO / path).open() as handle:
        return [row["gene"] for row in csv.DictReader(handle)]


def named_cases(with_gbm: bool = False) -> list[Json]:
    rows, samples = _breast_names()
    cluster = _seed_file("examples/breast/seeds_breast_cluster.csv")
    first10 = rows[:10]
    randomise_defaults = {"D.lower": 0.1, "D.higher": 0.24, "O.lower": 0.01, "sd.threshold": 0.0}
    spearman = {"correlation.threshold": 0.5, "correlation.method": "spearman"}
    cases: list[Json] = []

    def search(
        case_id: str,
        note: str,
        seeds: list[str] = cluster,
        universe: list[str] | None = None,
        assoc: Json | None = None,
        tolerance: float = 0.0,
        extra: Json | None = None,
        family: str | None = None,
        **args: Any,
    ) -> None:
        full = {
            "scale": True,
            **randomise_defaults,
            "keep.Alist": True,
            "max.iter": 10,
            "verbose": False,
        }
        full.update({k.replace("_", "."): v for k, v in args.items()})
        case: Json = {
            "id": case_id,
            "kind": "search",
            "family": family
            or ("regression, Breastdata1" if tolerance else "SEARCH(), Breastdata1"),
            "note": note,
            "data": {"dataset": "Breastdata1"},
            "seeds": seeds,
            "universe": universe or rows,
            "args": full,
            "assoc_args": dict(spearman if assoc is None else assoc),
        }
        if tolerance:
            case["tolerance"] = tolerance
        case.update(extra or {})
        cases.append(case)

    def glm(function: str = "nodes.glm", **assoc: Any) -> Json:
        out: Json = {"significance.threshold": 1e-05, "glm.family": "gaussian", "mc.cores": 8}
        out.update({k.replace("_", "."): v for k, v in assoc.items()})
        return out

    doc = {"D_lower": 0.0, "D_higher": 0.0, "O_lower": 0.1}
    search("breast_single_pass", "one pass with the thresholds of SEARCH.randomise", max_iter=1)
    search("breast_growth", "growth to convergence")
    search(
        "breast_doc_example",
        "the correlation example of the manual, through SEARCH()",
        seeds=first10,
        max_iter=3,
        **doc,
    )
    search(
        "breast_sd_prefilter",
        "every default of SEARCH.randomise, sd prefilter included (R1)",
        sd_threshold=0.1,
    )
    search("breast_seeds_not_immune", "keep.seeds = FALSE (R2, R3)", keep_Alist=False)
    search(
        "breast_doc_not_immune",
        "the manual's example with keep.seeds = FALSE (R2, R3)",
        seeds=first10,
        max_iter=3,
        keep_Alist=False,
        **doc,
    )
    search(
        "breast_pearson",
        "Pearson instead of Spearman",
        assoc={"correlation.threshold": 0.5, "correlation.method": "pearson"},
    )
    search(
        "breast_kendall",
        "Kendall instead of Spearman",
        assoc={"correlation.threshold": 0.4, "correlation.method": "kendall"},
        max_iter=2,
    )
    search(
        "breast_unscaled",
        "scale = FALSE: D is a count",
        scale=False,
        D_lower=1.0,
        D_higher=2.0,
        max_iter=3,
    )
    search("breast_universe_subset", "universe of the first 600 probesets", universe=rows[:600])
    search(
        "breast_universe_reversed",
        "the universe listed in reverse order, with the sd prefilter (R1)",
        universe=rows[::-1],
        sd_threshold=0.1,
    )
    search(
        "breast_pruned_members",
        "lower D.higher: seeds are added, then removed on D (R4)",
        D_higher=0.1,
    )
    search(
        "breast_omax_pruning",
        "higher O.lower: new seeds are removed on Omax (R2, R4)",
        D_higher=0.1,
        O_lower=0.3,
    )
    na_rows = [r for r in rows[24::50] if r not in cluster]
    na_cells = [[rows.index(r), 2] for r in na_rows]
    search(
        "breast_missing_values",
        "one missing value in 25 probesets (R7)",
        D_higher=0.1,
        extra={"na_cells": na_cells},
    )
    search(
        "breast_missing_pairwise",
        "the same with correlation.use = pairwise.complete.obs",
        D_higher=0.1,
        max_iter=2,
        assoc={**spearman, "correlation.use": "pairwise.complete.obs"},
        extra={"na_cells": na_cells},
    )
    search("breast_one_seed", "a single seed: D is divided by zero", seeds=cluster[:1], max_iter=2)
    search(
        "breast_glm_doc",
        "the GLM example of the manual, through SEARCH() (R5)",
        seeds=first10,
        max_iter=2,
        nodes_association_function="nodes.glm",
        assoc=glm(),
        tolerance=REGRESSION_TOLERANCE,
        **doc,
    )
    search(
        "breast_fastlm_doc",
        "the same with nodes.fastLm (R5)",
        seeds=first10,
        max_iter=2,
        nodes_association_function="nodes.fastLm",
        assoc=glm(),
        tolerance=REGRESSION_TOLERANCE,
        **doc,
    )
    search(
        "breast_fastlm",
        "nodes.fastLm, six seeds",
        max_iter=1,
        nodes_association_function="nodes.fastLm",
        assoc=glm(),
        tolerance=REGRESSION_TOLERANCE,
    )
    search(
        "breast_glm_bh",
        "nodes.glm with adjust.p = BH",
        max_iter=1,
        nodes_association_function="nodes.glm",
        assoc=glm(significance_threshold=0.01, adjust_p="BH"),
        tolerance=REGRESSION_TOLERANCE,
    )
    covariates = {"covariates": _breast_covariates(samples, {"t.rfs": "t.rfs", "e.rfs": "e.rfs"})}
    search(
        "breast_glm_covariates",
        "nodes.glm with two covariates through formula.string",
        max_iter=1,
        nodes_association_function="nodes.glm",
        assoc=glm(significance_threshold=0.001, formula_string="t.rfs + e.rfs"),
        tolerance=REGRESSION_TOLERANCE,
        extra=covariates,
    )
    search(
        "breast_glm_covariates_ignored",
        "covariates without formula.string: left out of the model (D3)",
        max_iter=1,
        nodes_association_function="nodes.glm",
        assoc=glm(significance_threshold=0.001),
        tolerance=REGRESSION_TOLERANCE,
        extra=covariates,
    )
    search(
        "breast_glm_interaction",
        "nodes.glm with an interaction between the gene and a covariate",
        max_iter=1,
        nodes_association_function="nodes.glm",
        assoc=glm(significance_threshold=0.001, formula_string="e.rfs + x.node:e.rfs"),
        tolerance=REGRESSION_TOLERANCE,
        extra=covariates,
    )
    search(
        "breast_glm_poisson",
        "poisson family at the documented threshold: no association passes (R6)",
        max_iter=1,
        nodes_association_function="nodes.glm",
        assoc=glm(glm_family="poisson"),
        tolerance=REGRESSION_TOLERANCE,
    )
    search(
        "breast_glm_poisson_loose",
        "poisson family, threshold 0.9",
        max_iter=1,
        nodes_association_function="nodes.glm",
        assoc=glm(glm_family="poisson", significance_threshold=0.9),
        tolerance=REGRESSION_TOLERANCE,
    )
    search(
        "breast_unknown_function",
        "an association function that does not exist",
        max_iter=1,
        nodes_association_function="nodes.spline",
    )
    search(
        "breast_unknown_correlation_type",
        "correlation.type other than standard",
        max_iter=1,
        assoc={**spearman, "correlation.type": "partial"},
    )
    search(
        "breast_unknown_family",
        "a GLM family the package does not list",
        max_iter=1,
        nodes_association_function="nodes.glm",
        assoc=glm(glm_family="binomial"),
        tolerance=REGRESSION_TOLERANCE,
    )
    search(
        "breast_seed_not_in_matrix",
        "a seed that is not a row of the matrix, given to SEARCH() directly",
        seeds=[*cluster, "not_a_probeset"],
        max_iter=1,
    )
    for fam in ("Gamma", "inverse.gaussian", "quasi", "quasipoisson"):
        search(
            f"breast_glm_{fam.replace('.', '_').lower()}",
            f"nodes.glm, {fam} family",
            max_iter=1,
            nodes_association_function="nodes.glm",
            assoc=glm(glm_family=fam),
            tolerance=REGRESSION_TOLERANCE,
        )

    def randomise(
        case_id: str,
        note: str,
        n_boot: int,
        seeds: list[str] = cluster,
        null: tuple[int, int] = (2, 5),
        sampling: str = "bootstrap",
        sampling_seed: int = 1,
        real_dir: str = "real/",
        sample_kind: str | None = None,
        assoc: Json | None = None,
        tolerance: float = 0.0,
        summarise_iter: int | None = None,
        **args: Any,
    ) -> None:
        full = {
            "scale": True,
            **randomise_defaults,
            "keep.seeds": True,
            "max.iter": 1,
            "verbose": False,
            "n.boot": n_boot,
            "sampling": sampling,
            "sampling.seed": sampling_seed,
        }
        full.update({k.replace("_", "."): v for k, v in args.items()})
        rng = np.random.default_rng(sum(map(ord, case_id)))
        others = [r for r in rows if r not in seeds]
        null_runs = [
            {
                "dir": f"null/batch_{b + 1}_0_FALSE_0_0/",
                "seeds": [
                    others[i] for i in rng.choice(len(others), size=len(seeds), replace=False)
                ],
                "n.boot": null[1],
            }
            for b in range(null[0])
        ]
        case: Json = {
            "id": case_id,
            "kind": "randomise",
            "family": "regression, Breastdata1" if tolerance else "SEARCH.randomise(), Breastdata1",
            "note": note,
            "data": {"dataset": "Breastdata1"},
            "seeds": seeds,
            "universe": rows,
            "args": full,
            "assoc_args": dict(spearman if assoc is None else assoc),
            "null_runs": null_runs,
            "real_dir": real_dir,
            "summarise": {
                "max.iter": summarise_iter or full["max.iter"],
                "n.boot": n_boot,
                "n.boot.NULL": null[1],
                "batches.NULL": null[0],
            },
        }
        if sample_kind:
            case["sample_kind"] = sample_kind
        if tolerance:
            case["tolerance"] = tolerance
        cases.append(case)

    randomise(
        "breast_bootstrap",
        "five replicates, two iterations, null runs and the summary",
        5,
        max_iter=2,
    )
    randomise(
        "breast_bootstrap_6iter",
        "eight replicates, up to six iterations: they stop at different ones (R9)",
        8,
        null=(1, 3),
        max_iter=6,
    )
    randomise(
        "breast_bootstrap_defaults",
        "every default of SEARCH.randomise, ten replicates (R1, R8)",
        10,
        sd_threshold=0.1,
    )
    randomise(
        "breast_jackknife",
        "sampling = jacknife",
        4,
        null=(1, 3),
        sampling="jacknife",
        sampling_seed=7,
    )
    randomise(
        "breast_bootstrap_rounding",
        "the bootstrap as R drew it before version 3.6.0",
        5,
        max_iter=2,
        sample_kind="Rounding",
    )
    randomise(
        "breast_manual_correlation",
        "the correlation example of the manual, as written there",
        2,
        seeds=first10,
        null=(2, 2),
        D_lower=0.0,
        D_higher=0.0,
        O_lower=0.1,
        max_iter=3,
    )
    randomise(
        "breast_blank_in_path",
        "an output directory with a blank in its name (R10)",
        3,
        real_dir="real run/",
        max_iter=2,
    )
    randomise(
        "breast_pearson_bootstrap",
        "Pearson on resampled columns",
        3,
        null=(2, 2),
        max_iter=2,
        assoc={"correlation.threshold": 0.5, "correlation.method": "pearson"},
    )
    randomise(
        "breast_manual_glm",
        "the GLM example of the manual, as written there",
        2,
        seeds=first10,
        null=(2, 2),
        D_lower=0.0,
        D_higher=0.0,
        O_lower=0.1,
        max_iter=2,
        nodes_association_function="nodes.glm",
        assoc={"significance.threshold": 1e-05, "glm.family": "gaussian", "mc.cores": 8},
        tolerance=REGRESSION_TOLERANCE,
    )

    # Gene names that read.table takes for numbers (R defect R13): summarise.iterations reads
    # them back changed. A small matrix with a planted block, so that the runs have neighbours.
    rng = np.random.default_rng(13)
    signal = rng.normal(0, 1, size=10)
    block = np.round(
        8
        + np.outer(np.array([1, 1, 1, -1, 1, 0, 0, 0]), signal)
        + rng.normal(0, 0.4, size=(8, 10)),
        2,
    )
    for case_id, names, note in (
        (
            "names_whole_numbers",
            ["101", "007", "20", "3", "45", "6", "88", "9"],
            "leading zeros are lost: 007 comes back as 7",
        ),
        (
            "names_decimal_numbers",
            ["101", "1.50", "100000", "3", "45", "6", "88", "9"],
            "1.50 comes back as 1.5 and 100000 as 1e+05",
        ),
        (
            "names_colliding_numbers",
            ["7", "07", "20", "3", "45", "6", "88", "9"],
            "7 and 07 become the same name: the summary stops",
        ),
        (
            "names_logical",
            ["T", "F", "TRUE", "FALSE", "T2", "g6", "g7", "g8"],
            "not all logical: left as they are",
        ),
        (
            "names_all_logical",
            ["T", "F", "TRUE", "FALSE"],
            "T comes back as TRUE: duplicated names, the summary stops",
        ),
        (
            "names_with_na",
            ["NA", "g2", "g3", "g4", "g5", "g6", "g7", "g8"],
            "a gene called NA: the summary stops",
        ),
    ):
        values = block[: len(names)]
        cases.append(
            {
                "id": case_id,
                "kind": "randomise",
                "family": "SEARCH.randomise(), gene names that look like numbers",
                "note": note,
                "data": _data(values, names, [f"s{j + 1}" for j in range(10)]),
                "seeds": names[:3],
                "universe": names,
                "args": {
                    "scale": True,
                    "D.lower": 0.0,
                    "D.higher": 0.24,
                    "O.lower": 0.01,
                    "keep.seeds": True,
                    "max.iter": 1,
                    "verbose": False,
                    "sd.threshold": 0.0,
                    "n.boot": 3,
                    "sampling": "bootstrap",
                    "sampling.seed": 3,
                },
                "assoc_args": dict(spearman),
                "null_runs": [
                    {
                        "dir": f"null/batch_{b + 1}_0_FALSE_0_0/",
                        "seeds": names[-3:] if b else names[3:6][: len(names) - 3] or names[-3:],
                        "n.boot": 2,
                    }
                    for b in range(2)
                ],
                "summarise": {"max.iter": 1, "n.boot": 3, "n.boot.NULL": 2, "batches.NULL": 2},
            }
        )

    if with_gbm:
        gbm = _seed_file("examples/gbm/seeds_buffa2010_setA.csv")
        for case_id, note, extra_args in (
            ("it_growth", "GBM IT cohort, ten hypoxia seeds, three iterations", {"max.iter": 3}),
            (
                "it_sd_prefilter",
                "GBM IT cohort, one pass with the sd prefilter",
                {"max.iter": 1, "sd.threshold": 0.1},
            ),
        ):
            cases.append(
                {
                    "id": case_id,
                    "kind": "search",
                    "family": "SEARCH(), GBM IT cohort",
                    "note": note,
                    "data": {"dataset": "IT"},
                    "seeds": gbm,
                    "universe": None,
                    "args": {
                        "scale": True,
                        **randomise_defaults,
                        "keep.Alist": True,
                        "verbose": False,
                        **extra_args,
                    },
                    "assoc_args": dict(spearman),
                }
            )
    return cases


FAMILIES = ("search", "randomise", "regression", "named")


def cases_of(family: str, n_fuzz: int, seed: int, with_gbm: bool = False) -> list[Json]:
    if family == "search":
        return fuzz_search_cases(n_fuzz, seed)
    if family == "randomise":
        return fuzz_randomise_cases(max(n_fuzz // 8, 1), seed + 1)
    if family == "regression":
        return fuzz_glm_cases(max(n_fuzz // 8, 1), seed + 2)
    if family == "named":
        return named_cases(with_gbm)
    raise ValueError(family)
