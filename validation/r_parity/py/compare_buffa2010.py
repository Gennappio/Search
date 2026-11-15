"""Parity harness, Python side: the port (``buffa2010`` preset) against the lab's GBM R script.

The R side (``run_gbm_script.sh``) runs the script unmodified and exports its tables and its
random draws. Here the unmodified pipeline of ``search`` is run twice on the same matrices:

1. **same draws**: the three functions that produce random numbers (``resample_indices`` for
   the bootstrap, ``derive_rng`` and ``random_seed_set`` for the null model) are replaced by
   lookups in the draws of the R run. Every replicate then works on the same resample, or the
   same random seed set, as in R, and every table is compared number by number.
2. **own draws**: the pipeline as shipped, with its own random generator and the same number of
   replicates. Agreement is then statistical: it measures what a different random stream does.

Usage::

    SEARCH_GBM_DATA=<GBM project folder> python validation/r_parity/py/compare_buffa2010.py \
        <repository root> <output directory> [n_jobs for the own-draws run]
"""

from __future__ import annotations

import os
import sys
import warnings
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats as sps

sys.path.insert(0, str(Path(__file__).parent))

from common import (  # noqa: E402
    Row,
    compare_scalars,
    compare_sets,
    compare_vectors,
    print_report,
    read_lines,
    write_report,
)

import search.bootstrap as search_bootstrap  # noqa: E402
import search.null as search_null  # noqa: E402
import search.run as search_run  # noqa: E402
from search.checkpoint import ChunkStore, chunk_ranges  # noqa: E402
from search.config import RunConfig  # noqa: E402
from search.core.corepass import CoreResult, run_core  # noqa: E402
from search.data import ExpressionSet  # noqa: E402
from search.presets import get_preset  # noqa: E402
from search.run import CohortManifest, CohortRun, PipelineResult, run_pipeline  # noqa: E402

TRAINING = ("IT", "FR")
CORE_COLUMNS = [
    "C_raw",
    "C_rank",
    "n_edges",
    "n_positive_edges",
    "n_negative_edges",
    "positive_consistent",
    "weighted_seed_rho",
    "is_seed",
]
STAGE_COLUMNS = ["C_rank_boot_mean", "C_rank_boot_var", "C_rank_boot_sd", "p_MC_C"]
# column of the R script -> column of search
COMBINED_COLUMNS = {
    **{
        f"{c}_{d}": f"{c}_{d}"
        for d in TRAINING
        for c in (
            "C_raw",
            "C_rank",
            "C_rank_boot_mean",
            "C_rank_boot_var",
            "C_rank_boot_sd",
            "p_MC_C",
            "positive_consistent",
            "weighted_seed_rho",
        )
    },
    "var_IT_safe": "var_IT_safe",
    "var_FR_safe": "var_FR_safe",
    "is_seed": "is_seed",
    "pC_hat_ITxFR_obs": "pC_obs",
    "pC_hat_ITxFR_boot": "pC_boot",
    "C_hat_IVW_ITFR": "C_IVW",
    "positive_consistent_ITFR": "positive_consistent_all",
    "MC_significant_ITFR": "MC_significant_all",
    "Buffa_product_pass": "product_pass",
    "Buffa_exact_candidate": "candidate",
}


# ---- reading the R export ---------------------------------------------------------------------
def read_r_table(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for column in table.columns:
        values = table[column]
        if values.isin(["TRUE", "FALSE", "NA"]).all():
            table[column] = values.map({"TRUE": 1.0, "FALSE": 0.0, "NA": np.nan})
        else:
            converted = pd.to_numeric(values.replace("NA", "nan"), errors="coerce")
            if converted.notna().sum() == (values != "NA").sum():
                table[column] = converted
    return table


def numeric(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.astype(float)
    return pd.to_numeric(series, errors="coerce").astype(float)


def flatten(frame: pd.DataFrame) -> pd.Series:
    """genes x columns -> one value per 'gene|column'."""
    stacked = frame.stack(future_stack=True)
    stacked.index = [f"{g}|{c}" for g, c in stacked.index]
    return stacked.astype(float)


def read_gene_matrix(path: Path) -> pd.DataFrame:
    return read_r_table(path).set_index("gene").astype(float)


# ---- the pipeline --------------------------------------------------------------------------------
def make_manifest(repo: Path, gbm: Path) -> CohortManifest:
    data_dir = gbm / "DATA" / "NORMALIZED"
    spec: dict[str, Any] = {
        "name": "gbm_parity",
        "seeds": str(repo / "examples" / "gbm" / "seeds_buffa2010_setA.csv"),
        "gene_filter": {"strip_ensembl_version": True, "common_genes_across_training": True},
        "cohorts": [
            {
                "name": "IT",
                "role": "training",
                "format": "csv",
                "matrix": str(data_dir / "IT_GBM_RNA_Deseq2_vst_Combat.csv"),
            },
            {
                "name": "FR",
                "role": "training",
                "format": "csv",
                "matrix": str(data_dir / "FR_GBM_RNA_Deseq2_vst_Combat.csv"),
            },
            {
                "name": "TCGA",
                "role": "scoring",
                "format": "csv",
                "matrix": str(data_dir / "TCGA_GBM_RNA_Deseq2_vst.csv"),
            },
        ],
        "config": {"preset": "buffa2010"},
    }
    return CohortManifest.model_validate({**spec, "base_dir": repo})


class RDraws:
    def __init__(self, export: Path, cohorts: tuple[str, ...] = TRAINING) -> None:
        self.boot = {
            c: np.loadtxt(export / f"{c}_bootstrap_indices_1based.tsv", dtype=int, ndmin=2) - 1
            for c in cohorts
        }
        self.sets = {
            c: [line.split("\t") for line in read_lines(export / f"{c}_random_seed_sets.tsv")]
            for c in cohorts
        }
        self.eligible = {c: set(read_lines(export / f"{c}_eligible_genes.txt")) for c in cohorts}
        self.eligible_python: dict[str, set[str]] = {}


@contextmanager
def r_draws(draws: RDraws) -> Iterator[None]:
    """Hand the draws of the R run to the pipeline; nothing else is replaced."""
    current: dict[str, str] = {}
    original = (
        search_run.run_cohort,
        search_bootstrap.resample_indices,
        search_null.derive_rng,
        search_null.random_seed_set,
    )

    def run_cohort(expr: Any, seeds: Any, cfg: Any, name: str, run_dir: Any, **kw: Any) -> Any:
        current["cohort"] = name
        return original[0](expr, seeds, cfg, name, run_dir, **kw)

    def resample_indices(scheme: str, n_samples: int, seed_base: int, rep: int) -> Any:
        idx = draws.boot[current["cohort"]][rep - 1]
        assert scheme == "bootstrap" and len(idx) == n_samples
        return np.asarray(idx, dtype=np.intp)

    def derive_rng(seed_base: int, rep: int) -> Any:
        return rep  # the replicate number travels in place of the generator

    def random_seed_set(eligible: Any, k: int, rng: Any) -> list[str]:
        draws.eligible_python[current["cohort"]] = {str(g) for g in eligible}
        chosen = draws.sets[current["cohort"]][int(rng) - 1]
        assert len(chosen) == k
        return list(chosen)

    search_run.run_cohort = run_cohort  # type: ignore[assignment]
    search_bootstrap.resample_indices = resample_indices  # type: ignore[assignment]
    search_null.derive_rng = derive_rng  # type: ignore[assignment]
    search_null.random_seed_set = random_seed_set  # type: ignore[assignment]
    try:
        yield
    finally:
        (
            search_run.run_cohort,
            search_bootstrap.resample_indices,
            search_null.derive_rng,
            search_null.random_seed_set,
        ) = original


def python_chunks(
    run_dir: Path, cfg: RunConfig, cohort: str, stage: str, n_rep: int
) -> pd.DataFrame:
    size = cfg.resampling.chunk_size if stage == "bootstrap" else cfg.null.chunk_size
    store = ChunkStore(run_dir / "checkpoints" / cohort, stage, cfg.checkpoint_dtype)
    frame = store.read_all(len(chunk_ranges(n_rep, size)))
    frame.columns = [str(int(c.split("_")[1])) for c in frame.columns]
    return frame


# ---- comparisons -----------------------------------------------------------------------------
def compare_core(
    ds: str, variant: str, export: Path, core: CoreResult, r_name: str | None = None
) -> list[Row]:
    """The deterministic pass: correlations, edges, memberships, overlaps, connectivity."""
    rows: list[Row] = []
    r_name = r_name or ds
    genes = [str(g) for g in core.genes]
    r_table = read_r_table(export / f"{r_name}_gene_connectivity.tsv").set_index("gene")
    py_table = core.gene_table().set_index("gene")
    rows.append(compare_sets(ds, variant, "genes of the training matrix", r_table.index, genes))
    for name, py_matrix in (
        ("rho", core.stat),
        ("delta", core.delta.astype(float)),
        ("gamma", core.gamma),
    ):
        r_matrix = read_gene_matrix(export / f"{r_name}_core_{name}.tsv")
        py_frame = pd.DataFrame(py_matrix, index=genes, columns=core.seeds)
        rows.append(
            compare_vectors(
                ds, variant, f"core: {name} (genes x seeds)", flatten(r_matrix), flatten(py_frame)
            )
        )
    sizes = read_r_table(export / f"{r_name}_cluster_sizes.tsv").set_index("seed")["cluster_size"]
    rows.append(
        compare_vectors(
            ds,
            variant,
            "core: cluster size of every seed",
            sizes.astype(float),
            pd.Series(core.cluster_size.astype(float), index=core.seeds),
        )
    )
    for name, py_matrix in (("S", core.S), ("F", core.F)):
        r_matrix = read_gene_matrix(export / f"{r_name}_{name}.tsv")
        py_frame = pd.DataFrame(py_matrix, index=core.seeds, columns=core.seeds)
        rows.append(
            compare_vectors(
                ds, variant, f"core: {name} between seeds", flatten(r_matrix), flatten(py_frame)
            )
        )
    for column in CORE_COLUMNS:
        rows.append(
            compare_vectors(
                ds, variant, f"core: {column}", numeric(r_table[column]), numeric(py_table[column])
            )
        )
    return rows


def equal_cluster_sizes(
    expr: ExpressionSet, seeds: list[str], cfg: RunConfig, boot: np.ndarray
) -> list[int]:
    """Bootstrap replicates (1-based) in which two seeds have clusters of the same size.

    There two different combinations of seeds can give a gene the same connectivity. The two
    sums are equal in exact arithmetic and may differ in the last bit in floating point, so
    whether the genes are tied in the ranking depends on the order of the additions, which is
    not the same in R and in NumPy.
    """
    hits = []
    for b in range(len(boot)):
        sizes = run_core(expr, seeds, cfg, sample_idx=boot[b]).cluster_size.tolist()
        if len(set(sizes)) < len(sizes):
            hits.append(b + 1)
    return hits


def compare_cohort(
    ds: str,
    variant: str,
    export: Path,
    run: CohortRun,
    run_dir: Path,
    cfg: RunConfig,
    draws: RDraws,
    b_boot: int,
    n_mc: int,
    expr: ExpressionSet,
) -> list[Row]:
    """One dataset: the core pass, every bootstrap replicate, every random seed set."""
    rows = compare_core(ds, variant, export, run.core)
    rows.append(
        compare_sets(
            ds,
            variant,
            "genes eligible as random seeds",
            draws.eligible[ds],
            draws.eligible_python.get(ds, set()),
        )
    )
    r_table = read_r_table(export / f"{ds}_gene_connectivity.tsv").set_index("gene")
    py_table = run.gene_table().set_index("gene")
    for column in STAGE_COLUMNS:
        stage = "bootstrap summary" if "boot" in column else "null model"
        rows.append(
            compare_vectors(
                ds,
                variant,
                f"{stage}: {column}",
                numeric(r_table[column]),
                numeric(py_table[column]),
            )
        )
    summary_rows = [r for r in rows if r.quantity.startswith("bootstrap summary")]
    r_boot = read_gene_matrix(export / f"{ds}_bootstrap_C_rank.tsv")
    r_boot.columns = [c.replace("boot_", "") for c in r_boot.columns]
    py_boot = python_chunks(run_dir, cfg, ds, "bootstrap", b_boot)
    tied = [str(b) for b in equal_cluster_sizes(expr, run.core.seeds, cfg, draws.boot[ds])]
    plain = [c for c in r_boot.columns if c not in tied]
    row = compare_vectors(
        ds,
        variant,
        "bootstrap: C_rank of every gene, replicates where the seed clusters differ in size",
        flatten(r_boot[plain]),
        flatten(py_boot[plain]),
    )
    row.detail = (f"{len(plain)} of {b_boot} replicates. " + row.detail).strip()
    rows.append(row)
    if tied:
        row = compare_vectors(
            ds,
            variant,
            "bootstrap: C_rank of every gene, "
            "replicates where two seed clusters have the same size",
            flatten(r_boot[tied]),
            flatten(py_boot[tied]),
        )
        explanation = (
            "connectivities equal in exact arithmetic differ in the last bit, and R and NumPy "
            "split different ties. "
            if row.status == "differs"
            else "no tie is split differently by R and NumPy. "
        )
        row.detail = (f"replicates {', '.join(tied)}: " + explanation + row.detail).strip()
        rows.append(row)
        for summary in summary_rows:
            if summary.status == "differs":
                summary.detail = (
                    f"follows from the tied replicates {', '.join(tied)}. " + summary.detail
                )
    r_mc = read_gene_matrix(export / f"{ds}_mc_C_raw.tsv")
    r_mc.columns = [c.replace("mc_", "") for c in r_mc.columns]
    row = compare_vectors(
        ds,
        variant,
        "null model: C_raw of every gene for every random seed set",
        flatten(r_mc),
        flatten(python_chunks(run_dir, cfg, ds, "null", n_mc)),
    )
    row.detail = (f"{n_mc} random seed sets. " + row.detail).strip()
    rows.append(row)
    assert run.null is not None
    r_summary = read_r_table(export / f"{ds}_mc_random_seed_summary.tsv").set_index("mc")
    py_summary = run.null.random_seed_summary.set_index("rep")
    for column in ("mean_S", "median_S"):
        rows.append(
            compare_vectors(
                ds,
                variant,
                f"null model: {column} of every random seed set",
                pd.Series(
                    r_summary[column].to_numpy(dtype=float), index=r_summary.index.astype(str)
                ),
                pd.Series(
                    py_summary[column].to_numpy(dtype=float), index=py_summary.index.astype(str)
                ),
            )
        )
    return rows


def compare_same_draws(export: Path, result: PipelineResult, draws: RDraws) -> list[Row]:
    rows: list[Row] = []
    variant = "same draws"
    info = dict(line.split("\t", 1) for line in read_lines(export / "run_info.tsv"))
    b_boot, n_mc = int(info["b_boot"]), int(info["n_mc"])
    for ds in TRAINING:
        rows.extend(
            compare_cohort(
                ds,
                variant,
                export,
                result.cohorts[ds],
                result.run_dir,
                result.config,
                draws,
                b_boot,
                n_mc,
                result.inputs.cohorts[ds],
            )
        )

    assert result.combined is not None
    r_full = read_r_table(export / "combined_full_ranking.tsv")
    py_full = result.combined.table
    r_indexed, py_indexed = r_full.set_index("gene"), py_full.set_index("gene")
    tied_upstream = any("two seed clusters have the same size" in r.quantity for r in rows)
    for r_column, py_column in COMBINED_COLUMNS.items():
        row = compare_vectors(
            "IT x FR",
            variant,
            f"combination: {r_column}",
            numeric(r_indexed[r_column]),
            numeric(py_indexed[py_column]),
        )
        if row.status == "differs" and tied_upstream:
            row.detail = "follows from the tied bootstrap replicates. " + row.detail
        rows.append(row)
    same_order = list(r_full["gene"]) == list(py_full["gene"])
    first = next(
        (
            i
            for i, (a, b) in enumerate(zip(r_full["gene"], py_full["gene"], strict=False))
            if a != b
        ),
        None,
    )
    rows.append(
        Row(
            "IT x FR",
            variant,
            "combination: row order of the full ranking",
            len(r_full),
            0 if same_order else 1,
            0.0,
            "equal" if same_order else "differs",
            ""
            if same_order
            else (
                ("follows from the tied bootstrap replicates. " if tied_upstream else "")
                + f"first different row: {first + 1 if first is not None else '-'}"
            ),
        )
    )
    for r_name, py_frame, label in (
        ("combined_signature.tsv", result.combined.signature, "signature"),
        ("combined_signature_top51.tsv", result.combined.top, "top 51 of the signature"),
    ):
        r_genes = list(read_r_table(export / r_name)["gene"])
        py_genes = [str(g) for g in py_frame["gene"]]
        row = compare_sets(
            "IT x FR", variant, f"selection: genes of the {label}", r_genes, py_genes
        )
        if row.status == "equal":
            row.detail += "; same order" if r_genes == py_genes else "; different order"
            if r_genes != py_genes:
                row.status, row.n_mismatch = "differs", 1
        rows.append(row)
    for r_name, key in (
        ("scores_all_signature.tsv", "signature_all"),
        ("scores_top51.tsv", "signature_top"),
    ):
        path = export / r_name
        if not path.exists() or key not in result.scores:
            rows.append(
                compare_scalars(
                    "scores",
                    variant,
                    f"score table {key} written",
                    path.exists(),
                    key in result.scores,
                )
            )
            continue
        r_scores, py_scores = read_r_table(path), result.scores[key]
        r_scores.index = r_scores["dataset"] + "|" + r_scores["sample"].astype(str)
        py_scores = py_scores.set_index(
            py_scores["dataset"] + "|" + py_scores["sample"].astype(str)
        )
        for column in (
            "E_median_abs_expression",
            "HS_rank01",
            "n_signature_genes_present",
            "n_signature_genes_total",
        ):
            row = compare_vectors(
                "scores",
                variant,
                f"{key}: {column} (IT, FR and TCGA samples)",
                numeric(r_scores[column]),
                numeric(py_scores[column]),
            )
            # the reports are kept in the repository: no sample identifier may end up in them
            row.detail = f"{row.n_mismatch} samples differ" if row.n_mismatch else ""
            rows.append(row)
    return rows


def compare_own_draws(export: Path, result: PipelineResult) -> list[Row]:
    """What changes when only the random stream changes: information, not pass or fail."""
    rows: list[Row] = []
    variant = "own draws"
    cutoff = result.config.combine.mc_p_cutoff
    for ds in TRAINING:
        r_table = read_r_table(export / f"{ds}_gene_connectivity.tsv").set_index("gene")
        py_table = result.cohorts[ds].gene_table().set_index("gene").reindex(r_table.index)
        for column in ("C_raw", "C_rank"):
            rows.append(
                compare_vectors(
                    ds,
                    variant,
                    f"core: {column} (no random numbers involved)",
                    numeric(r_table[column]),
                    numeric(py_table[column]),
                )
            )
        for column in ("C_rank_boot_mean", "p_MC_C"):
            a, b = numeric(r_table[column]), numeric(py_table[column])
            ok = a.notna() & b.notna()
            rho = float(sps.spearmanr(a[ok], b[ok]).statistic)
            diff = float(np.abs(a[ok] - b[ok]).max())
            detail = f"Spearman {rho:.4f}, largest absolute difference {diff:.4f}"
            if column == "p_MC_C":
                same_side = int(((a[ok] < cutoff) == (b[ok] < cutoff)).sum())
                detail += (
                    f"; p < {cutoff}: R {int((a[ok] < cutoff).sum())}, Python "
                    f"{int((b[ok] < cutoff).sum())}, same side of the cutoff {same_side} "
                    f"of {int(ok.sum())}"
                )
            rows.append(Row(ds, variant, column, int(ok.sum()), 0, diff, "info", detail))
    assert result.combined is not None
    r_genes = set(read_r_table(export / "combined_signature.tsv")["gene"])
    py_genes = {str(g) for g in result.combined.signature["gene"]}
    rows.append(
        Row(
            "IT x FR",
            variant,
            "selection: genes of the signature",
            len(r_genes | py_genes),
            0,
            0.0,
            "info",
            f"R {len(r_genes)}, Python {len(py_genes)}, shared {len(r_genes & py_genes)}",
        )
    )
    return rows


def main() -> None:
    warnings.filterwarnings("ignore")
    repo = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve()
    n_jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    export = out / "r_export"
    gbm = Path(os.environ["SEARCH_GBM_DATA"]).resolve()
    info = dict(line.split("\t", 1) for line in read_lines(export / "run_info.tsv"))
    cfg: RunConfig = get_preset("buffa2010").with_overrides(
        {
            "resampling.n_rep": int(info["b_boot"]),
            "null.n_rep": int(info["n_mc"]),
            "resampling.chunk_size": 10,
            "null.chunk_size": 10,
        }
    )
    manifest = make_manifest(repo, gbm)
    draws = RDraws(export)
    with r_draws(draws):
        same = run_pipeline(manifest, cfg, out / "python_same_draws", n_jobs=1)
    rows = compare_same_draws(export, same, draws)
    own = run_pipeline(manifest, cfg, out / "python_own_draws", n_jobs=n_jobs)
    rows.extend(compare_own_draws(export, own))
    for label, run in (("same draws", same), ("own draws", own)):
        timings = {
            name: {k: round(v, 1) for k, v in c.timings.items()} for name, c in run.cohorts.items()
        }
        rows.append(Row("timing", label, "seconds per stage", 0, 0, 0.0, "info", str(timings)))
    rows.append(
        Row(
            "timing",
            "R script",
            "whole run",
            0,
            0,
            0.0,
            "info",
            (out / "r_script_seconds.txt").read_text().strip(),
        )
    )
    frame = write_report(rows, out / "report_buffa2010.csv")
    print_report(frame)
    print("\n", frame.groupby(["variant", "status"]).size().to_string())


if __name__ == "__main__":
    main()
