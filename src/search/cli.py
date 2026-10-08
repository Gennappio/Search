"""Command line interface."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Any

import typer

from search import __version__
from search.errors import SearchError

app = typer.Typer(no_args_is_help=True, help="Seed-based gene co-expression clustering.")
config_app = typer.Typer(no_args_is_help=True, help="Inspect and validate configurations.")
app.add_typer(config_app, name="config")


def _print_version(value: bool) -> None:
    if value:
        typer.echo(__version__)
        raise typer.Exit()


def _fail(exc: Exception) -> None:
    typer.secho(f"error: {exc}", fg=typer.colors.RED, err=True)
    raise typer.Exit(code=1)


def _overrides(items: list[str] | None) -> dict[str, Any]:
    from search.config import parse_override

    return dict(parse_override(item) for item in items or [])


@app.callback()
def main(
    version: bool = typer.Option(
        False, "--version", callback=_print_version, is_eager=True, help="Print the version."
    ),
) -> None:
    """Search: seed-based gene co-expression clustering with named presets."""


@app.command()
def version() -> None:
    """Print the installed version."""
    typer.echo(__version__)


@app.command()
def run(
    manifest: Annotated[Path, typer.Option("--manifest", "-m", help="Cohort manifest (YAML).")],
    preset: Annotated[str | None, typer.Option(help="Preset name, overrides the manifest.")] = None,
    set_: Annotated[
        list[str] | None, typer.Option("--set", help="Override key.path=value (YAML value).")
    ] = None,
    out: Annotated[Path, typer.Option(help="Directory holding the run folders.")] = Path("runs"),
    n_jobs: Annotated[int | None, typer.Option(help="Worker processes for replicates.")] = None,
    resume: Annotated[bool, typer.Option(help="Reuse checkpoints of the latest run.")] = False,
    strict: Annotated[
        bool, typer.Option(help="Refuse configurations outside a validated preset.")
    ] = False,
    name: Annotated[str | None, typer.Option(help="Run name, overrides the manifest.")] = None,
) -> None:
    """Run every stage of the pipeline on the cohorts of a manifest."""
    from search.run import CohortManifest, resolve_run_dir, run_pipeline

    try:
        spec = CohortManifest.load(manifest)
        if name:
            spec = spec.model_copy(update={"name": name})
        cfg = spec.build_config(preset, _overrides(set_))
        if strict:
            cfg = cfg.with_overrides({"strict": True})
        run_dir = None
        if resume and (out / "latest" / "run_path.txt").exists():
            run_dir = resolve_run_dir(out / "latest")
        result = run_pipeline(spec, cfg, out, resume=resume, run_dir=run_dir, n_jobs=n_jobs)
    except (SearchError, NotImplementedError) as exc:
        _fail(exc)
        return
    typer.echo(f"run directory: {result.run_dir}")
    if result.combined is not None:
        typer.echo(f"candidate genes: {len(result.combined.signature)}")


@app.command()
def resume(
    run_dir: Annotated[
        Path, typer.Argument(help="Run directory, or the 'latest' directory.")
    ] = Path("runs/latest"),
    n_jobs: Annotated[int | None, typer.Option(help="Worker processes for replicates.")] = None,
) -> None:
    """Resume an interrupted run, reusing its checkpoints and configuration."""
    from search.config import RunConfig
    from search.export.manifest import read_manifest
    from search.run import CohortManifest, resolve_run_dir, run_pipeline

    try:
        target = resolve_run_dir(run_dir)
        stored = read_manifest(target)
        manifest_path = stored.get("manifest_path")
        if not manifest_path:
            raise SearchError("the run manifest does not record the cohort manifest path")
        spec = CohortManifest.load(manifest_path)
        cfg = RunConfig.from_yaml((target / "config.yaml").read_text())
        result = run_pipeline(spec, cfg, target.parent, resume=True, run_dir=target, n_jobs=n_jobs)
    except (SearchError, NotImplementedError, OSError) as exc:
        _fail(exc)
        return
    typer.echo(f"run directory: {result.run_dir}")


@app.command()
def qc(manifest: Annotated[Path, typer.Option("--manifest", "-m")]) -> None:
    """Load the cohorts of a manifest and print the quality warnings."""
    from search.qc import format_warnings
    from search.run import CohortManifest, load_inputs

    try:
        inputs = load_inputs(CohortManifest.load(manifest))
    except SearchError as exc:
        _fail(exc)
        return
    for name, expr in inputs.cohorts.items():
        typer.echo(
            f"{name}: {expr.n_genes} genes x {expr.n_samples} samples ({inputs.roles[name]})"
        )
        typer.echo(format_warnings(inputs.qc[name]))
    missing = [
        s for s in inputs.seeds.seeds if any(s not in e.genes for e in inputs.cohorts.values())
    ]
    typer.echo(f"seeds: {len(inputs.seeds.seeds)}, missing in some cohort: {missing or 'none'}")


@app.command()
def combine(
    run_dir: Annotated[Path, typer.Argument(help="Run directory, or 'latest'.")] = Path(
        "runs/latest"
    ),
    mc_p_cutoff: Annotated[float | None, typer.Option()] = None,
    product_cutoff: Annotated[float | None, typer.Option()] = None,
    top_n: Annotated[int | None, typer.Option()] = None,
    positive_consistent: Annotated[
        bool | None, typer.Option("--positive-consistent/--no-positive-consistent")
    ] = None,
    exclude_seeds: Annotated[bool | None, typer.Option("--exclude-seeds/--keep-seeds")] = None,
    prefix: Annotated[str, typer.Option(help="File name prefix of the new tables.")] = "recombined",
) -> None:
    """Recombine the cohort tables of a run with different selection criteria."""
    from search.combine import combine_cohorts
    from search.config import RunConfig
    from search.export.tables import read_cohort_tables, write_combined_tables
    from search.run import resolve_run_dir

    try:
        target = resolve_run_dir(run_dir)
        cfg = RunConfig.from_yaml((target / "config.yaml").read_text()).combine
        updates = {
            "mc_p_cutoff": mc_p_cutoff,
            "product_cutoff": product_cutoff,
            "top_n": top_n,
            "require_positive_consistent": positive_consistent,
            "exclude_seeds": exclude_seeds,
        }
        cfg = cfg.model_copy(update={k: v for k, v in updates.items() if v is not None})
        tables = read_cohort_tables(target / "results")
        if not tables:
            raise SearchError(f"no cohort tables under {target / 'results'}")
        symbols = None
        full = target / "results" / "combined_full_ranking.csv"
        if full.exists():
            import pandas as pd

            previous = pd.read_csv(full)
            if "symbol" in previous.columns:
                symbols = dict(
                    zip(previous["gene"].astype(str), previous["symbol"].astype(str), strict=True)
                )
        combined = combine_cohorts(tables, cfg, symbols)
        written = write_combined_tables(target / "results", combined, prefix)
    except (SearchError, OSError) as exc:
        _fail(exc)
        return
    typer.echo(
        f"candidate genes: {len(combined.signature)}; written: {', '.join(p.name for p in written)}"
    )


@app.command()
def score(
    manifest: Annotated[Path, typer.Option("--manifest", "-m")],
    signature: Annotated[Path, typer.Option(help="Table with a 'gene' column.")],
    out: Annotated[Path, typer.Option(help="Output CSV.")] = Path("scores.csv"),
    method: Annotated[str, typer.Option()] = "median_abs_rank01",
    min_genes: Annotated[int, typer.Option()] = 3,
) -> None:
    """Score every sample of every cohort with a signature."""
    import pandas as pd

    from search.config import ScoreConfig
    from search.data import read_gene_list
    from search.run import CohortManifest, load_inputs
    from search.signature import score_samples

    try:
        inputs = load_inputs(CohortManifest.load(manifest))
        genes = read_gene_list(signature)
        cfg = ScoreConfig(enabled=True, method=method, min_genes=min_genes)
        parts = [score_samples(expr, genes, cfg, name) for name, expr in inputs.cohorts.items()]
        table = pd.concat(parts, ignore_index=True)
        table.to_csv(out, index=False)
    except (SearchError, ImportError) as exc:
        _fail(exc)
        return
    typer.echo(f"scores for {len(table)} samples written to {out}")


@config_app.command("show")
def config_show(
    preset: Annotated[str | None, typer.Option(help="Preset name.")] = None,
    file: Annotated[Path | None, typer.Option(help="Configuration YAML.")] = None,
    set_: Annotated[list[str] | None, typer.Option("--set")] = None,
) -> None:
    """Print a configuration as YAML, with overrides applied."""
    from search.config import RunConfig
    from search.presets import get_preset, matching_preset

    try:
        cfg = (
            get_preset(preset)
            if preset
            else (RunConfig.from_yaml(file.read_text()) if file else RunConfig())
        )
        cfg = cfg.with_overrides(_overrides(set_))
    except (SearchError, OSError) as exc:
        _fail(exc)
        return
    typer.echo(cfg.to_yaml())
    typer.echo(f"# matches preset: {matching_preset(cfg)}")


@config_app.command("validate")
def config_validate(file: Annotated[Path, typer.Argument(help="Configuration YAML.")]) -> None:
    """Check that a configuration file is valid and say which preset it matches."""
    from search.config import RunConfig
    from search.presets import matching_preset

    try:
        cfg = RunConfig.from_yaml(file.read_text())
    except (SearchError, OSError) as exc:
        _fail(exc)
        return
    typer.echo(f"valid; matches preset: {matching_preset(cfg)}; hash {cfg.hash()}")


if __name__ == "__main__":  # pragma: no cover
    app()
