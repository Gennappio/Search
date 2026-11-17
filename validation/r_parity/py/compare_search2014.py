"""Parity harness, Python side: the port against the SEARCH 0.0.2 R package.

Reads the cases written by ``r/search2014_run.R``, runs the same cases through ``search`` and
compares every quantity the R package returns. Usage::

    python validation/r_parity/py/compare_search2014.py <repository root> <output directory>
"""

from __future__ import annotations

import json
import os
import sys
import time
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from common import (  # noqa: E402
    TOLERANCE,
    Row,
    compare_scalars,
    compare_sets,
    compare_vectors,
    print_report,
    read_lines,
    read_named,
    read_r_matrix,
    write_report,
)

from search.config import RunConfig  # noqa: E402
from search.core.connectivity import nodes_degree  # noqa: E402
from search.core.overlap import jaccard_binary  # noqa: E402
from search.data import ExpressionSet, read_matrix  # noqa: E402
from search.growth import GrowthResult, _membership_from_neighbourhoods, run_growth  # noqa: E402
from search.null import empirical_p, per_gene_counts  # noqa: E402
from search.presets import get_preset  # noqa: E402

ITERATIVE = 1e-8
"""Tolerance for the GLM families fitted by iteratively reweighted least squares: R and
statsmodels stop the iterations at different points."""
ROUNDING = 5.1e-6
"""write.data.files rounds D to five decimals before summarise.iterations reads it back."""

P_ADJUST = {None: "none", "BH": "bh", "bonferroni": "bonferroni"}
GLM_FAMILY = {
    "gaussian": "gaussian",
    "poisson": "poisson",
    "Gamma": "gamma",
    "inverse.gaussian": "inverse_gaussian",
}


def as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    return list(value) if isinstance(value, list) else [value]


def load_dataset(repo: Path, name: str) -> ExpressionSet:
    if name.startswith("Breast"):
        folder = repo / "tests" / "data" / "breast" / name
        return read_matrix(
            folder / "mRNA_abundance.txt.gz", sample_table=folder / "patient_annotation.txt"
        )
    gbm = Path(os.environ["SEARCH_GBM_DATA"])
    return read_matrix(gbm / "DATA" / "NORMALIZED" / f"{name}_GBM_RNA_Deseq2_vst_Combat.csv")


def config_for(
    case: dict[str, Any], method: str = "ols", with_covariates: bool = True
) -> tuple[RunConfig, dict[str, Any]]:
    """The search2014 preset plus what the R call sets differently from SEARCH.randomise."""
    args = case["assoc_args"]
    overrides: dict[str, Any] = {
        "growth.max_iter": case["max.iter"],
        "growth.D_lower": case["D.lower"],
        "growth.D_higher": case["D.higher"],
        "growth.O_lower": case["O.lower"],
        "growth.keep_seeds": case["keep.Alist"],
        "growth.sd_threshold": case["sd.threshold"],
        "connectivity.normalise": "fraction" if case["scale"] else "count",
        "resampling.scheme": "none",
        "null.scheme": "none",
    }
    if case["association"] == "nodes.correlation":
        overrides["association.method"] = args["correlation.method"]
        overrides["edge.rule"] = "abs_threshold"
        overrides["edge.threshold"] = args["correlation.threshold"]
    else:
        overrides["association.method"] = method
        overrides["association.glm_family"] = GLM_FAMILY[args.get("glm.family", "gaussian")]
        overrides["association.covariates"] = (
            as_list(case.get("covariate_names")) or as_list(case.get("covariates"))
            if with_covariates
            else []
        )
        overrides["edge.rule"] = "p_threshold"
        overrides["edge.alpha"] = args["significance.threshold"]
        overrides["edge.bonferroni"] = "none"
        overrides["edge.p_adjust"] = P_ADJUST[args.get("adjust.p")]
    preset = get_preset("search2014")
    cfg = preset.with_overrides(overrides)
    base = preset.model_dump()
    changed = {}
    for key, value in overrides.items():
        section, name = key.split(".")
        if base[section][name] != value:
            changed[key] = value
    return cfg, changed


# ---- pieces of a growth result ------------------------------------------------------------------
def overlap_table(growth: GrowthResult, seeds: list[str]) -> pd.Series:
    """O between the seeds that hold a neighbourhood, flattened to 'a|b' for a before b."""
    delta, _ = _membership_from_neighbourhoods(growth.final_universe, seeds, growth.neighbourhoods)
    o = jaccard_binary(delta)
    out = {}
    for a, sa in enumerate(seeds):
        for b in range(a, len(seeds)):
            out[f"{sa}|{seeds[b]}"] = o[a, b]
    return pd.Series(out, dtype=float)


def r_overlap_table(path: Path) -> tuple[pd.Series, list[str]]:
    matrix = read_r_matrix(path)
    seeds = [str(s) for s in matrix.index]
    values = matrix.to_numpy(dtype=float)
    out = {}
    for a, sa in enumerate(seeds):
        for b in range(a, len(seeds)):
            out[f"{sa}|{seeds[b]}"] = values[a, b]
    return pd.Series(out, dtype=float), seeds


def degree_table(growth: GrowthResult) -> pd.Series:
    seeds = list(growth.neighbourhoods)
    universe = growth.final_universe
    delta, _ = _membership_from_neighbourhoods(universe, seeds, growth.neighbourhoods)
    seed_idx = universe.get_indexer(pd.Index(seeds))
    return pd.Series(nodes_degree(delta, seed_idx), index=seeds)


def neighbourhood_table(growth: GrowthResult) -> pd.Series:
    out = {}
    universe = set(growth.final_universe)
    for seed, hood in growth.neighbourhoods.items():
        for gene, value in hood.items():
            if gene in universe:
                out[f"{seed}|{gene}"] = float(value)
    return pd.Series(out, dtype=float)


def r_iterations(rdir: Path) -> list[int]:
    return sorted(int(p.stem.split("_")[-1]) for p in rdir.glob("D_iter_*.tsv"))


def compare_growth(
    label: str,
    variant: str,
    rdir: Path,
    expr: ExpressionSet,
    seeds: list[str],
    cfg: RunConfig,
    *,
    sample_idx: np.ndarray | None = None,
    universe: list[str] | None = None,
    light: bool = False,
    prefix: str = "",
) -> tuple[list[Row], GrowthResult | None]:
    rows: list[Row] = []
    iterative = cfg.association.method == "glm" and cfg.association.glm_family != "gaussian"
    tolerance = ITERATIVE if iterative else TOLERANCE
    growth = run_growth(expr, seeds, cfg, universe=universe, sample_idx=sample_idx)
    py_iters = [it for it in growth.iterations if len(it.D)]
    if (rdir / "error.txt").exists():
        message = (rdir / "error.txt").read_text().strip()
        rows.append(
            Row(
                label,
                variant,
                prefix + "run",
                0,
                1,
                float("nan"),
                "r_error",
                f"R stopped with '{message}'; Python: {len(py_iters)} iterations, "
                f"{growth.stop_reason}, {len(growth.final_seeds)} seeds",
            )
        )
        return rows, growth
    iters = r_iterations(rdir)
    rows.append(compare_scalars(label, variant, prefix + "iterations", len(iters), len(py_iters)))
    all_names = {str(g) for g in growth.universe}
    last_names: set[str] = all_names
    for i in iters:
        r_d = read_named(rdir / f"D_iter_{i}.tsv", "D")
        last_names = set(r_d.index)
        if i > len(py_iters):
            continue
        it = py_iters[i - 1]
        py_d = pd.Series(it.D.to_numpy(dtype=float), index=[str(g) for g in it.D.index])
        rows.append(compare_vectors(label, variant, f"{prefix}D, iteration {i}", r_d, py_d))
        r_added = read_lines(rdir / f"added_iter_{i}.txt")
        row = compare_sets(label, variant, f"{prefix}seeds added, iteration {i}", r_added, it.added)
        if row.status == "equal" and r_added != list(it.added):
            row.detail += "; same genes, different order"
        rows.append(row)
    r_removed_d = set(read_lines(rdir / "removed_by_D.txt"))
    py_removed_d = {g for it in growth.iterations for g in it.removed_D}
    rows.append(
        compare_sets(label, variant, prefix + "seeds removed on D", r_removed_d, py_removed_d)
    )
    py_low_sd = {str(r["gene"]) for r in growth.removed if r["reason"] == "low_sd"}
    py_removed_o = {g for it in growth.iterations for g in it.removed_O}
    r_gone = all_names - last_names - r_removed_d
    if cfg.growth.sd_threshold > 0:
        # R does not return what the prefilter or Omax removed: both are what left the universe
        rows.append(
            compare_sets(
                label,
                variant,
                prefix + "genes removed by the sd prefilter or on Omax",
                r_gone,
                py_low_sd | py_removed_o,
            )
        )
    else:
        rows.append(
            compare_sets(label, variant, prefix + "seeds removed on Omax", r_gone, py_removed_o)
        )
    r_sizes = read_named(rdir / "final_cluster_sizes.tsv", "size")
    rows.append(
        compare_sets(
            label,
            variant,
            prefix + "seeds holding a neighbourhood at the end",
            r_sizes.index,
            growth.neighbourhoods.keys(),
        )
    )
    if light:
        return rows, growth
    path = rdir / "final_neighbourhoods.tsv"
    if path.exists():
        table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        r_hoods = pd.Series(
            pd.to_numeric(table["value"].replace("NA", "nan")).to_numpy(dtype=float),
            index=table["seed"] + "|" + table["gene"],
        )
        rows.append(
            compare_vectors(
                label,
                variant,
                prefix + "neighbourhood members and association values",
                r_hoods,
                neighbourhood_table(growth),
                tolerance,
            )
        )
    path = rdir / "O_final.tsv"
    if path.exists():
        r_o, r_seeds = r_overlap_table(path)
        kept = [s for s in r_seeds if s in growth.neighbourhoods]
        rows.append(
            compare_vectors(
                label,
                variant,
                prefix + "O of the final neighbourhoods",
                r_o,
                overlap_table(growth, kept) if len(kept) >= 1 else pd.Series(dtype=float),
            )
        )
    for path in sorted(rdir.glob("O_iter_*.tsv")):
        i = int(path.stem.split("_")[-1])
        partial = run_growth(
            expr,
            seeds,
            cfg.with_overrides({"growth.max_iter": i}),
            universe=universe,
            sample_idx=sample_idx,
        )
        r_o, r_seeds = r_overlap_table(path)
        kept = [s for s in r_seeds if s in partial.neighbourhoods]
        rows.append(
            compare_vectors(
                label,
                variant,
                f"{prefix}O matrix stored at iteration {i}",
                r_o,
                overlap_table(partial, kept) if kept else pd.Series(dtype=float),
            )
        )
        degree_path = rdir / f"degree_iter_{i}.tsv"
        if degree_path.exists():
            r_degree = read_named(degree_path, "degree")
            rows.append(
                compare_vectors(
                    label,
                    variant,
                    f"{prefix}nodes degree, iteration {i}",
                    r_degree,
                    degree_table(partial),
                )
            )
    return rows, growth


def compare_bootstrap(case: dict[str, Any], case_dir: Path, expr: ExpressionSet) -> list[Row]:
    rows: list[Row] = []
    label = case["id"]
    cfg, _ = config_for(case)
    seeds = [str(s) for s in case["seeds"]]
    draws = np.loadtxt(case_dir / "sample_draws_1based.tsv", dtype=int, ndmin=2) - 1
    n_boot = int(case["n.boot"])
    for variant in ("verbatim", "patched"):
        rdir = case_dir / f"r_{variant}"
        if not rdir.is_dir():
            continue
        if (rdir / "error.txt").exists():
            rows.append(
                Row(
                    label,
                    variant,
                    "run",
                    0,
                    1,
                    float("nan"),
                    "r_error",
                    (rdir / "error.txt").read_text().strip(),
                )
            )
            continue
        per_rep: list[GrowthResult] = []
        worst: dict[str, Row] = {}
        for i in range(1, n_boot + 1):
            rep_rows, growth = compare_growth(
                label,
                variant,
                rdir / f"rep_{i}",
                expr,
                seeds,
                cfg,
                sample_idx=draws[i - 1],
                light=True,
            )
            assert growth is not None
            per_rep.append(growth)
            for row in rep_rows:
                key = row.quantity.split(", iteration")[0]
                kept = worst.get(key)
                if kept is None:
                    worst[key] = Row(
                        label,
                        variant,
                        f"per replicate: {key}",
                        row.n,
                        row.n_mismatch,
                        row.max_abs_diff,
                        row.status,
                        "" if row.status == "equal" else f"replicate {i}: {row.detail}",
                    )
                else:
                    kept.n += row.n
                    kept.n_mismatch += row.n_mismatch
                    if not np.isnan(row.max_abs_diff):
                        kept.max_abs_diff = max(kept.max_abs_diff, row.max_abs_diff)
                    if row.status != "equal" and kept.status == "equal":
                        kept.status = row.status
                        kept.detail = f"replicate {i}: {row.detail}"
        for row in worst.values():
            row.detail = (f"{n_boot} replicates on the samples drawn by R. " + row.detail).strip()
            rows.append(row)

        if (rdir / "summarise_error.txt").exists():
            rows.append(
                Row(
                    label,
                    variant,
                    "summarise.iterations",
                    0,
                    1,
                    float("nan"),
                    "r_error",
                    "R stopped with '" + (rdir / "summarise_error.txt").read_text().strip() + "'",
                )
            )
        # summarise.iterations: mean over the replicates, per growth iteration
        for path in sorted(rdir.glob("summary_real_iter_*.tsv")):
            k = int(path.stem.split("_")[-1])
            real = read_r_matrix(path)
            r_mean = real["mean"]
            columns = []
            for growth in per_rep:
                done = [it for it in growth.iterations if len(it.D)]
                if len(done) >= k:
                    d = done[k - 1].D
                    columns.append(
                        pd.Series(d.to_numpy(dtype=float), index=[str(g) for g in d.index])
                    )
            stacked = pd.concat(columns, axis=1) if columns else pd.DataFrame()
            py_mean = stacked.mean(axis=1, skipna=True).dropna()
            n_reached = len(columns)
            row = compare_vectors(
                label,
                variant,
                f"bootstrap mean of D, iteration {k}",
                r_mean.dropna(),
                py_mean,
                ROUNDING,
            )
            row.detail = (
                f"{n_reached} of {n_boot} replicates reached iteration {k}; R reads D back from "
                "files rounded to five decimals. " + row.detail
            ).strip()
            rows.append(row)
            if k == int(case["max.iter"]):
                finals = []
                for growth in per_rep:
                    finals.append(
                        pd.Series(growth.final_D(), index=[str(g) for g in growth.universe])
                    )
                production = pd.concat(finals, axis=1).mean(axis=1, skipna=True).dropna()
                row = compare_vectors(
                    label,
                    variant,
                    f"bootstrap mean as the pipeline takes it (final D), against R iteration {k}",
                    r_mean.dropna(),
                    production,
                    ROUNDING,
                )
                last = [len([it for it in g.iterations if len(it.D)]) for g in per_rep]
                row.detail = (f"last iteration of each replicate: {last}. " + row.detail).strip()
                rows.append(row)
            null_path = rdir / f"summary_null_iter_{k}.tsv"
            if null_path.exists():
                null = read_r_matrix(null_path).reindex(real.index)
                observed = r_mean.to_numpy(dtype=float)
                chunk = null.to_numpy(dtype=float)
                count, n_null = per_gene_counts(observed, chunk)
                py_p = pd.Series(empirical_p(count, n_null), index=real.index)
                complete = ~np.isnan(chunk).any(axis=1) & ~np.isnan(observed)
                groups = {
                    "genes with a mean and present in every null replicate": complete,
                    "genes absent from some null replicate": ~complete & ~np.isnan(observed),
                    "genes without a mean (pruned in every replicate)": np.isnan(observed),
                }
                for name, mask in groups.items():
                    if not mask.any():
                        continue
                    row = compare_vectors(
                        label,
                        variant,
                        f"null p-value, iteration {k}: {name}",
                        real["P"][mask],
                        py_p[mask],
                    )
                    row.detail = (
                        f"R's means and null values through per_gene_counts; {chunk.shape[1]} "
                        "null replicates. " + row.detail
                    ).strip()
                    rows.append(row)
    return rows


def main() -> None:
    warnings.filterwarnings("ignore")  # statsmodels warns on every perfect self-regression
    repo = Path(sys.argv[1]).resolve()
    out_root = Path(sys.argv[2]).resolve()
    only = set(sys.argv[3:])
    rows: list[Row] = []
    mapping: dict[str, Any] = {}
    cache: dict[str, ExpressionSet] = {}
    for case_dir in sorted(p for p in out_root.iterdir() if (p / "case.json").exists()):
        case = json.loads((case_dir / "case.json").read_text())
        if only and case["id"] not in only:
            continue
        if not case["dataset"].startswith("Breast") and not os.environ.get("SEARCH_GBM_DATA"):
            continue
        t0 = time.perf_counter()
        if case["dataset"] not in cache:
            cache[case["dataset"]] = load_dataset(repo, case["dataset"])
        expr = cache[case["dataset"]]
        renamed = as_list(case.get("covariate_names"))
        if renamed:
            assert expr.sample_table is not None
            table = expr.sample_table.rename(
                columns=dict(zip(as_list(case["covariates"]), renamed, strict=True))
            )
            expr = ExpressionSet(expr.values, expr.genes, expr.samples, expr.gene_table, table)
        seeds = [str(s) for s in case["seeds"]]
        na_rows = as_list(case.get("na_rows"))
        if na_rows:
            values = expr.values.copy()
            values[expr.gene_positions(na_rows), 2] = np.nan  # the third sample, as in R
            expr = ExpressionSet(
                values, expr.genes, expr.samples, expr.gene_table, expr.sample_table
            )
        universe: list[str] | None = as_list(case.get("universe")) or None
        if universe:
            # SEARCH() adds the seeds to the universe itself; search asks for them explicitly
            universe = universe + [s for s in seeds if s not in set(universe)]
        if case["kind"] == "bootstrap":
            rows.extend(compare_bootstrap(case, case_dir, expr))
            mapping[case["id"]] = config_for(case)[1]
        else:
            # which configurations of search are set against this R case
            runs: list[tuple[str, str, bool]] = [("", "ols", True)]
            if case["association"] != "nodes.correlation":
                gaussian = case["assoc_args"].get("glm.family", "gaussian") == "gaussian"
                if not gaussian:
                    runs = [("", "glm", True)]
                elif int(case["max.iter"]) == 1:
                    runs.append((" [association.method=glm]", "glm", True))
                if as_list(case.get("covariates")) and "formula.string" not in case["assoc_args"]:
                    runs = [
                        (" [covariates declared]", "ols", True),
                        (" [no covariates]", "ols", False),
                    ]
            for suffix, method, with_covariates in runs:
                cfg, changed = config_for(case, method, with_covariates)
                label = case["id"] + suffix
                mapping[label] = changed
                for variant in ("verbatim", "patched"):
                    rdir = case_dir / f"r_{variant}"
                    if rdir.is_dir():
                        got, _ = compare_growth(
                            label,
                            variant,
                            rdir,
                            expr,
                            seeds,
                            cfg,
                            universe=universe,
                            light=bool(case["light"]),
                        )
                        rows.extend(got)
        print(f"{case['id']:<28} compared in {time.perf_counter() - t0:.1f} s", flush=True)
    frame = write_report(rows, out_root / "report_search2014.csv")
    (out_root / "config_overrides.json").write_text(json.dumps(mapping, indent=1))
    print_report(frame)
    counts = frame.groupby(["variant", "status"]).size()
    print("\n", counts.to_string())


if __name__ == "__main__":
    main()
