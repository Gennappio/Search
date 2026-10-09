"""The cohort manifest and the pipeline that runs every stage for every cohort."""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
import yaml
from pydantic import BaseModel, ConfigDict, Field

from search.bootstrap import BootstrapSummary, run_bootstrap, statistic_of
from search.combine import CombinedResult, combine_cohorts
from search.config import RunConfig
from search.core.corepass import CoreResult, resolve_seeds, run_core
from search.data import (
    ExpressionSet,
    SeedSpec,
    read_gene_list,
    read_matrix,
    read_seeds,
    read_table,
)
from search.errors import ConfigError, DataError
from search.export.manifest import expression_hash, file_sha256, now_iso, versions, write_manifest
from search.export.tables import (
    cohort_gene_table,
    write_cohort_tables,
    write_combined_tables,
    write_growth_tables,
    write_scores,
)
from search.growth import GrowthResult, growth_statistic, run_growth
from search.logging_utils import close_file_handlers, get_logger
from search.null import NullResult, run_internal_null
from search.parallel import default_n_jobs
from search.presets import get_preset, matching_preset
from search.qc import QC_INFORMATIONAL, QCWarning, format_warnings, qc_expression
from search.signature import score_samples

log = logging.getLogger("search")


# ---- manifest ----------------------------------------------------------------------------
class _Spec(BaseModel):
    model_config = ConfigDict(extra="forbid")


class CohortSpec(_Spec):
    name: str
    role: Literal["training", "scoring"] = "training"
    matrix: str
    format: Literal["csv", "tsv", "parquet", "h5ad"] | None = None
    orientation: Literal["genes_x_samples", "samples_x_genes"] = "genes_x_samples"
    sample_table: str | None = None


class GeneFilterSpec(_Spec):
    strip_ensembl_version: bool = True
    common_genes_across_training: bool = True
    annotation_table: str | None = None
    annotation_query: str | None = None


class ConfigSpec(_Spec):
    preset: str | None = None
    file: str | None = None
    overrides: dict[str, Any] = Field(default_factory=dict)


class CohortManifest(_Spec):
    name: str = "search_run"
    seeds: str
    universe: str | None = None
    gene_filter: GeneFilterSpec = Field(default_factory=GeneFilterSpec)
    cohorts: list[CohortSpec]
    config: ConfigSpec = Field(default_factory=ConfigSpec)
    base_dir: Path = Field(default=Path("."), exclude=True)
    source: Path | None = Field(default=None, exclude=True)

    @classmethod
    def load(cls, path: Path | str) -> CohortManifest:
        p = Path(path)
        try:
            data = yaml.safe_load(p.read_text()) or {}
        except (OSError, yaml.YAMLError) as exc:
            raise ConfigError(f"cannot read manifest {p}: {exc}") from exc
        try:
            manifest = cls.model_validate({**data, "base_dir": p.parent, "source": p})
        except ValueError as exc:
            raise ConfigError(f"invalid manifest {p}: {exc}") from exc
        names = [c.name for c in manifest.cohorts]
        if len(set(names)) != len(names):
            raise ConfigError("cohort names must be unique")
        if not any(c.role == "training" for c in manifest.cohorts):
            raise ConfigError("at least one cohort must have role 'training'")
        return manifest

    def resolve(self, path: str) -> Path:
        p = Path(path).expanduser()
        return p if p.is_absolute() else (self.base_dir / p).resolve()

    def build_config(
        self, preset: str | None = None, overrides: dict[str, Any] | None = None
    ) -> RunConfig:
        """Preset or config file from the manifest, manifest overrides, then CLI overrides."""
        chosen = preset or self.config.preset
        if chosen:
            cfg = get_preset(chosen)
        elif self.config.file:
            cfg = RunConfig.from_yaml(self.resolve(self.config.file).read_text())
        else:
            cfg = RunConfig()
        if self.config.overrides:
            cfg = cfg.with_overrides(self.config.overrides)
        if overrides:
            cfg = cfg.with_overrides(overrides)
        return cfg


# ---- inputs ------------------------------------------------------------------------------
@dataclass
class LoadedInputs:
    cohorts: dict[str, ExpressionSet]
    roles: dict[str, str]
    seeds: SeedSpec
    universe: list[str] | None
    symbols: dict[str, str]
    qc: dict[str, list[QCWarning]]
    sources: dict[str, dict[str, Any]]


def load_inputs(manifest: CohortManifest) -> LoadedInputs:
    seeds = read_seeds(manifest.resolve(manifest.seeds))
    universe: list[str] | None = None
    if manifest.universe:
        universe = read_gene_list(manifest.resolve(manifest.universe))
    symbols = dict(seeds.symbols)
    gf = manifest.gene_filter
    if gf.annotation_table:
        annotation = read_table(manifest.resolve(gf.annotation_table))
        symbol_column = next(
            (c for c in ("symbol", "symbol_consensus", "gene_symbol") if c in annotation.columns),
            None,
        )
        if symbol_column:
            found = annotation[symbol_column].dropna()
            symbols = {**{str(k): str(v) for k, v in found.items()}, **symbols}
        if gf.annotation_query:
            try:
                selected = annotation.query(gf.annotation_query)
            except Exception as exc:  # pandas raises many kinds
                raise ConfigError(f"annotation_query failed: {exc}") from exc
            allowed = [str(g) for g in selected.index]
            universe = allowed if universe is None else [g for g in universe if g in set(allowed)]
    cohorts: dict[str, ExpressionSet] = {}
    roles: dict[str, str] = {}
    sources: dict[str, dict[str, Any]] = {}
    for spec in manifest.cohorts:
        matrix_path = manifest.resolve(spec.matrix)
        table_path = manifest.resolve(spec.sample_table) if spec.sample_table else None
        expr = read_matrix(
            matrix_path,
            format=spec.format,
            orientation=spec.orientation,
            strip_version=gf.strip_ensembl_version,
            sample_table=table_path,
        )
        cohorts[spec.name] = expr
        roles[spec.name] = spec.role
        sources[spec.name] = {
            "matrix": str(matrix_path),
            "matrix_sha256": file_sha256(matrix_path),
            "sample_table": str(table_path) if table_path else None,
            "n_genes_loaded": expr.n_genes,
            "n_samples": expr.n_samples,
            "role": spec.role,
        }
    training = [name for name, role in roles.items() if role == "training"]
    if gf.common_genes_across_training and len(training) > 1:
        common = ExpressionSet.common_genes(*(cohorts[n] for n in training))
        if len(common) == 0:
            raise DataError("training cohorts share no genes")
        for name in training:
            cohorts[name] = cohorts[name].subset_genes(common)
            sources[name]["n_genes_common"] = len(common)
    qc = {name: qc_expression(expr, name) for name, expr in cohorts.items()}
    return LoadedInputs(cohorts, roles, seeds, universe, symbols, qc, sources)


# ---- one cohort ----------------------------------------------------------------------
@dataclass
class CohortRun:
    name: str
    role: str
    core: CoreResult
    bootstrap: BootstrapSummary | None
    null: NullResult | None
    input_hash: str
    timings: dict[str, float] = field(default_factory=dict)
    growth: GrowthResult | None = None

    def gene_table(self) -> pd.DataFrame:
        table = cohort_gene_table(self.name, self.core, self.bootstrap, self.null)
        if self.growth is not None:
            final = pd.DataFrame(
                {"gene": np.asarray(self.growth.universe), "D_final": self.growth.final_D()}
            )
            table = table.merge(final, on="gene", how="left")
        return table


def run_cohort(
    expr: ExpressionSet,
    seeds: list[str],
    cfg: RunConfig,
    name: str,
    run_dir: Path,
    *,
    universe: list[str] | None = None,
    resume: bool = True,
    n_jobs: int = 1,
) -> CohortRun:
    timings: dict[str, float] = {}
    work, seed_list = resolve_seeds(expr, seeds, universe)
    input_hash = expression_hash(work, seed_list, None)
    ckpt_dir = run_dir / "checkpoints" / name

    t0 = time.perf_counter()
    growth: GrowthResult | None = None
    if cfg.growth.enabled:
        log.info(
            "%s: growth loop on %d genes x %d samples, %d seeds",
            name,
            work.n_genes,
            work.n_samples,
            len(seed_list),
        )
        growth = run_growth(work, seed_list, cfg)
        log.info(
            "%s: growth stopped (%s) with %d seeds after %d iterations",
            name,
            growth.stop_reason,
            len(growth.final_seeds),
            growth.n_iterations,
        )
        write_growth_tables(run_dir / "results", name, growth)
        if len(growth.final_seeds) < 2:
            raise DataError(
                f"{name}: growth ended with fewer than two seeds ({growth.stop_reason}); "
                "the growth tables under results/ show what was pruned"
            )
        core = run_core(work.subset_genes(growth.final_universe), growth.final_seeds, cfg)
        observed = growth_statistic(growth, cfg.null.statistic)
    else:
        log.info(
            "%s: core pass on %d genes x %d samples, %d seeds",
            name,
            work.n_genes,
            work.n_samples,
            len(seed_list),
        )
        core = run_core(work, seed_list, cfg)
        observed = statistic_of(core, cfg.null.statistic)
    timings["core"] = time.perf_counter() - t0

    bootstrap: BootstrapSummary | None = None
    if cfg.resampling.scheme != "none" and cfg.resampling.n_rep > 0:
        t0 = time.perf_counter()
        log.info("%s: %s with %d replicates", name, cfg.resampling.scheme, cfg.resampling.n_rep)
        bootstrap = run_bootstrap(
            work,
            seed_list,
            cfg,
            universe=None,
            ckpt_dir=ckpt_dir,
            input_hash=input_hash,
            resume=resume,
            n_jobs=n_jobs,
        )
        timings["bootstrap"] = time.perf_counter() - t0

    null: NullResult | None = None
    if cfg.null.scheme == "internal_random_seeds" and cfg.null.n_rep > 0:
        t0 = time.perf_counter()
        log.info(
            "%s: null model with %d random seed sets (%s)", name, cfg.null.n_rep, cfg.null.pooling
        )
        null = run_internal_null(
            work,
            seed_list,
            cfg,
            universe=None,
            observed=observed,
            ckpt_dir=ckpt_dir,
            input_hash=input_hash,
            resume=resume,
            n_jobs=n_jobs,
        )
        timings["null"] = time.perf_counter() - t0
    elif cfg.null.scheme == "external_runs":
        log.info("%s: null scheme external_runs, p-values come later from 'search summarise'", name)

    run = CohortRun(name, "training", core, bootstrap, null, input_hash, timings, growth)
    write_cohort_tables(run_dir / "results", name, core, bootstrap, null, run.gene_table())
    return run


def write_symbols(results_dir: Path, inputs: LoadedInputs) -> Path | None:
    """gene -> symbol for the genes of the loaded cohorts, used by the plots and the tables."""
    if not inputs.symbols:
        return None
    present: set[str] = set()
    for expr in inputs.cohorts.values():
        present.update(str(g) for g in expr.genes)
    rows = [(gene, symbol) for gene, symbol in inputs.symbols.items() if gene in present]
    if not rows:
        return None
    results_dir.mkdir(parents=True, exist_ok=True)
    path = results_dir / "gene_symbols.csv"
    pd.DataFrame(rows, columns=["gene", "symbol"]).to_csv(path, index=False)
    return path


# ---- the pipeline ----------------------------------------------------------------------
@dataclass
class PipelineResult:
    run_dir: Path
    config: RunConfig
    inputs: LoadedInputs
    cohorts: dict[str, CohortRun]
    combined: CombinedResult | None
    scores: dict[str, pd.DataFrame]
    manifest: dict[str, Any]


def new_run_dir(out_dir: Path, name: str) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = Path(out_dir) / f"{name}_{stamp}"
    counter = 1
    while run_dir.exists():
        counter += 1
        run_dir = Path(out_dir) / f"{name}_{stamp}_{counter}"
    run_dir.mkdir(parents=True)
    return run_dir


def resolve_run_dir(path: Path | str) -> Path:
    """A run directory, or a ``latest`` directory holding run_path.txt."""
    p = Path(path)
    pointer = p / "run_path.txt"
    if pointer.exists():
        return Path(pointer.read_text().strip())
    return p


def run_pipeline(
    manifest: CohortManifest,
    cfg: RunConfig,
    out_dir: Path | str,
    *,
    resume: bool = False,
    run_dir: Path | None = None,
    n_jobs: int | None = None,
) -> PipelineResult:
    out = Path(out_dir)
    if run_dir is None:
        run_dir = new_run_dir(out, manifest.name)
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    latest = out / "latest"
    latest.mkdir(parents=True, exist_ok=True)
    (latest / "run_path.txt").write_text(str(run_dir.resolve()))
    logger = get_logger(run_dir / "log.txt")
    jobs = n_jobs if n_jobs is not None else (cfg.n_jobs if cfg.n_jobs > 1 else default_n_jobs())
    preset_name = matching_preset(cfg)
    started = now_iso()
    manifest_dict: dict[str, Any] = {
        "name": manifest.name,
        "status": "running",
        "started": started,
        "run_dir": str(run_dir.resolve()),
        "manifest_path": str(manifest.source.resolve()) if manifest.source else None,
        "versions": versions(),
        "config": cfg.model_dump(mode="json"),
        "config_hash": cfg.hash(),
        "preset_matched": preset_name,
        "preset_validated": preset_name is not None,
        "n_jobs": jobs,
        "resume": resume,
    }
    (run_dir / "config.yaml").write_text(cfg.to_yaml())
    write_manifest(run_dir, manifest_dict)
    logger.info(
        "run %s in %s (preset %s, validated %s)",
        manifest.name,
        run_dir,
        preset_name,
        preset_name is not None,
    )
    try:
        if preset_name is None:
            if cfg.strict:
                raise ConfigError("strict mode: the configuration matches no validated preset")
            logger.warning("configuration matches no validated preset; results are not validated")
        inputs = load_inputs(manifest)
        manifest_dict["inputs"] = inputs.sources
        manifest_dict["seeds"] = inputs.seeds.seeds
        manifest_dict["qc"] = {name: [w.__dict__ for w in ws] for name, ws in inputs.qc.items()}
        blocking: dict[str, list[str]] = {}
        for name, warnings in inputs.qc.items():
            if warnings:
                logger.warning("%s %s", name, format_warnings(warnings))
            else:
                logger.info("%s QC: no warnings", name)
            codes = sorted({w.code for w in warnings if w.code not in QC_INFORMATIONAL})
            if codes:
                blocking[name] = codes
        manifest_dict["qc_blocking"] = blocking
        if cfg.null.scheme == "external_runs":
            manifest_dict["notes"] = [
                "null scheme external_runs: p-values come from 'search summarise --null-runs' "
                "against runs with random seed sets"
            ]
        write_symbols(run_dir / "results", inputs)
        write_manifest(run_dir, manifest_dict)
        if cfg.strict and blocking:
            detail = "; ".join(f"{name}: {', '.join(codes)}" for name, codes in blocking.items())
            raise ConfigError(f"strict mode: QC warnings on the data ({detail})")

        cohort_runs: dict[str, CohortRun] = {}
        for name, role in inputs.roles.items():
            if role != "training":
                continue
            cohort_runs[name] = run_cohort(
                inputs.cohorts[name],
                inputs.seeds.seeds,
                cfg,
                name,
                run_dir,
                universe=inputs.universe,
                resume=resume,
                n_jobs=jobs,
            )
            manifest_dict.setdefault("timings", {})[name] = cohort_runs[name].timings
            write_manifest(run_dir, manifest_dict)

        combined: CombinedResult | None = None
        scores: dict[str, pd.DataFrame] = {}
        if cfg.combine.enabled and cohort_runs:
            tables = {name: run.gene_table() for name, run in cohort_runs.items()}
            combined = combine_cohorts(tables, cfg.combine, inputs.symbols)
            write_combined_tables(run_dir / "results", combined)
            manifest_dict["combine"] = {
                "cohorts": combined.cohorts,
                "n_candidates": int(len(combined.signature)),
                "n_top": int(len(combined.top)),
                "notes": combined.notes,
            }
            logger.info(
                "combined %s: %d candidate genes",
                "+".join(combined.cohorts),
                len(combined.signature),
            )
            if cfg.score.enabled:
                for label, frame in (
                    ("signature_all", combined.signature),
                    ("signature_top", combined.top),
                ):
                    genes = [str(g) for g in frame["gene"]]
                    if len(genes) < cfg.score.min_genes:
                        logger.warning(
                            "score %s skipped: %d genes, fewer than %d",
                            label,
                            len(genes),
                            cfg.score.min_genes,
                        )
                        continue
                    parts = [
                        score_samples(expr, genes, cfg.score, name)
                        for name, expr in inputs.cohorts.items()
                    ]
                    scores[label] = pd.concat(parts, ignore_index=True)
                    write_scores(run_dir / "results", scores[label], f"scores_{label}.csv")
        manifest_dict["outputs"] = sorted(
            str(p.relative_to(run_dir)) for p in (run_dir / "results").glob("*") if p.is_file()
        )
        manifest_dict["status"] = "finished"
        manifest_dict["finished"] = now_iso()
        write_manifest(run_dir, manifest_dict)
        logger.info("finished: %s", run_dir)
        return PipelineResult(run_dir, cfg, inputs, cohort_runs, combined, scores, manifest_dict)
    except BaseException as exc:
        manifest_dict["status"] = "failed"
        manifest_dict["error"] = {"type": type(exc).__name__, "message": str(exc)}
        manifest_dict["finished"] = now_iso()
        write_manifest(run_dir, manifest_dict)
        raise
    finally:
        close_file_handlers()
