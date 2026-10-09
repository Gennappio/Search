"""Fixes after the usage audit: external null per cohort, strict QC, growth collapse, exports,
clean errors on missing files, symbols in plots, score method validation, breast example."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from typer.testing import CliRunner

from search.cli import app
from search.config import RunConfig
from search.core.corepass import run_core
from search.data import ExpressionSet, read_gene_list, read_matrix, read_seeds
from search.errors import ConfigError, DataError
from search.export.manifest import read_manifest
from search.export.tables import write_cohort_tables
from search.null import collect_null_chunks, external_null
from search.plots import load_symbols, node_labels
from search.presets import BUFFA2010, SEARCH2014
from search.run import CohortManifest, load_inputs, run_pipeline
from tests.synthetic_data import make_planted
from tests.unit.test_pipeline import write_manifest_files

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"


def test_external_null_is_restricted_to_the_same_cohort(tmp_path: Path) -> None:
    path = write_manifest_files(tmp_path / "proj", n_boot=4, n_null=0)
    manifest = CohortManifest.load(path)
    real = run_pipeline(manifest, manifest.build_config(), tmp_path / "runs")
    null_cfg = manifest.build_config(overrides={"resampling.seed_base": 777})
    null_run = run_pipeline(
        manifest.model_copy(update={"name": "nullrun"}), null_cfg, tmp_path / "runs"
    )
    table = real.cohorts["A"].gene_table()
    genes = pd.Index(table["gene"].astype(str))
    observed = table["C_rank_boot_mean"].to_numpy(dtype=float)
    mixed = external_null(observed, genes, [null_run.run_dir])
    only_a = external_null(observed, genes, [null_run.run_dir], cohort="A")
    assert mixed.n_rep == 8 and only_a.n_rep == 4
    assert (only_a.n_null <= 4).all() and (mixed.n_null <= 8).all()
    with pytest.raises(ValueError, match="no checkpoints for cohort"):
        collect_null_chunks([null_run.run_dir], cohort="Z")
    res = CliRunner().invoke(
        app, ["summarise", str(real.run_dir), "--null-runs", str(null_run.run_dir)]
    )
    assert res.exit_code == 0, res.output
    assert "against 4 null replicates" in res.output
    assert "A_external_null.csv" in res.output and "B_external_null.csv" in res.output


def _cohort_files(root: Path, sample_names: list[str] | None = None) -> Path:
    """One training cohort written to disk, optionally with chosen sample names."""
    root.mkdir(parents=True, exist_ok=True)
    data = make_planted(n_genes=40, n_samples=20, module_size=10, n_seeds=3, strength=2.0, seed=3)
    expr = data.expr
    if sample_names is not None:
        expr = ExpressionSet(expr.values, expr.genes, pd.Index(sample_names))
    expr.to_frame().to_csv(root / "A.csv")
    (root / "seeds.csv").write_text(
        "gene,symbol\n" + "\n".join(f"{g},SYM_{g}" for g in data.seeds) + "\n"
    )
    manifest = {
        "name": "strict",
        "seeds": "seeds.csv",
        "cohorts": [{"name": "A", "role": "training", "matrix": "A.csv"}],
        "config": {"preset": "buffa2010", "overrides": {"resampling.n_rep": 0, "null.n_rep": 0}},
    }
    (root / "cohorts.yaml").write_text(yaml.safe_dump(manifest))
    return root / "cohorts.yaml"


def test_strict_blocks_on_data_warnings_but_not_on_a_missing_sample_table(tmp_path: Path) -> None:
    names = [f"s{i}" for i in range(19)] + ["GbExp47_ctrl"]
    path = _cohort_files(tmp_path / "bad", names)
    manifest = CohortManifest.load(path)
    cfg = manifest.build_config(overrides={"strict": True})
    with pytest.raises(ConfigError, match="control_like_names"):
        run_pipeline(manifest, cfg, tmp_path / "runs_bad")
    run_dir = Path((tmp_path / "runs_bad" / "latest" / "run_path.txt").read_text())
    stored = read_manifest(run_dir)
    assert stored["status"] == "failed" and stored["error"]["type"] == "ConfigError"
    assert stored["qc_blocking"] == {"A": ["control_like_names"]}
    clean = CohortManifest.load(_cohort_files(tmp_path / "ok"))
    result = run_pipeline(
        clean, clean.build_config(overrides={"strict": True}), tmp_path / "runs_ok"
    )
    assert result.manifest["status"] == "finished" and result.manifest["qc_blocking"] == {}
    assert any(w.code == "sample_table_missing" for w in result.inputs.qc["A"])


def test_strict_on_unvalidated_config_marks_the_run_failed(tmp_path: Path) -> None:
    manifest = CohortManifest.load(_cohort_files(tmp_path / "proj"))
    cfg = manifest.build_config(overrides={"strict": True, "membership.mode": "binary"})
    with pytest.raises(ConfigError, match="no validated preset"):
        run_pipeline(manifest, cfg, tmp_path / "runs")
    run_dir = Path((tmp_path / "runs" / "latest" / "run_path.txt").read_text())
    assert read_manifest(run_dir)["status"] == "failed"


def test_growth_collapse_leaves_tables_and_an_error_in_the_manifest(tmp_path: Path) -> None:
    path = write_manifest_files(tmp_path / "proj", n_boot=0, n_null=0)
    manifest = CohortManifest.load(path)
    # thresholds nobody can satisfy: every seed is pruned at the first iteration
    cfg = manifest.build_config(
        preset="search2014",
        overrides={"growth.D_lower": 1.5, "growth.D_higher": 1.5, "growth.sd_threshold": 0.0},
    )
    with pytest.raises(DataError, match="fewer than two seeds"):
        run_pipeline(manifest, cfg, tmp_path / "runs")
    run_dir = Path((tmp_path / "runs" / "latest" / "run_path.txt").read_text())
    stored = read_manifest(run_dir)
    assert stored["status"] == "failed" and stored["error"]["type"] == "DataError"
    events = pd.read_csv(run_dir / "results" / "A_growth_events.csv")
    assert (events["event"] == "removed_D").sum() == 4
    assert (run_dir / "results" / "A_growth_iterations.csv").exists()


def test_seed_and_coefficient_exports(tmp_path: Path) -> None:
    rng = np.random.default_rng(5)
    n = 30
    table = pd.DataFrame(
        {"age": rng.normal(size=n), "batch": rng.choice(["x", "y"], n)},
        index=[f"s{j}" for j in range(n)],
    )
    expr = ExpressionSet(
        rng.normal(size=(12, n)),
        pd.Index([f"g{i}" for i in range(12)]),
        pd.Index(table.index),
        sample_table=table,
    )
    cfg = RunConfig().with_overrides(
        {"association.method": "ols", "association.covariates": ["age", "batch"]}
    )
    core = run_core(expr, ["g3", "g7"], cfg)
    assert core.coef_table is not None and list(core.coef_table["seed"]) == ["g3", "g7"]
    assert {"intercept", "age", "batch_y"} <= set(core.coef_table.columns)
    written = write_cohort_tables(tmp_path, "T", core, None, None)
    names = {p.name for p in written}
    assert {"T_seeds.csv", "T_covariate_coefficients.csv"} <= names
    seeds = pd.read_csv(tmp_path / "T_seeds.csv")
    assert list(seeds.columns) == ["dataset", "seed", "cluster_size", "nodes_degree"]
    assert list(seeds["seed"]) == ["g3", "g7"] and seeds["nodes_degree"].between(0, 1).all()
    coef = pd.read_csv(tmp_path / "T_covariate_coefficients.csv")
    assert list(coef["seed"]) == ["g3", "g7"] and coef.columns[0] == "dataset"


def test_missing_files_give_one_line_errors(tmp_path: Path) -> None:
    with pytest.raises(DataError, match="not found"):
        read_seeds(tmp_path / "nope.csv")
    with pytest.raises(DataError, match="not found"):
        read_gene_list(tmp_path / "nope.txt")
    with pytest.raises(DataError, match="not found"):
        read_matrix(tmp_path / "nope.txt.gz")
    path = write_manifest_files(tmp_path / "proj", n_boot=0, n_null=0)
    (tmp_path / "proj" / "seeds.csv").unlink()
    runner = CliRunner()
    for args in (
        ["run", "--manifest", str(path), "--out", str(tmp_path / "runs")],
        ["qc", "--manifest", str(path)],
    ):
        res = runner.invoke(app, args)
        assert res.exit_code == 1, res.output
        assert "not found" in res.output and "Traceback" not in res.output


def test_symbols_are_written_and_used_for_labels(tmp_path: Path) -> None:
    path = write_manifest_files(tmp_path / "proj", n_boot=0, n_null=0)
    manifest = CohortManifest.load(path)
    result = run_pipeline(manifest, manifest.build_config(), tmp_path / "runs")
    symbols = load_symbols(result.run_dir / "results")
    seeds = result.inputs.seeds.seeds
    assert all(symbols[s] == f"SYM_{s}" for s in seeds)
    labels = node_labels(seeds + ["unknown"], symbols)
    assert labels[seeds[0]] == f"SYM_{seeds[0]}" and labels["unknown"] == "unknown"
    assert load_symbols(tmp_path) == {}


def test_score_method_is_validated() -> None:
    assert BUFFA2010.score.method == "median_abs_rank01"
    ok = BUFFA2010.with_overrides({"score.method": "pysigscore_median"})
    assert ok.score.method == "pysigscore_median"
    with pytest.raises(ConfigError, match="score.method"):
        BUFFA2010.with_overrides({"score.method": "median"})
    with pytest.raises(ConfigError):
        BUFFA2010.with_overrides({"score.method": "pysigscore_"})


def test_breast_example_manifest_loads() -> None:
    manifest = CohortManifest.load(EXAMPLES / "breast" / "cohorts.yaml")
    assert manifest.build_config().preset == "search2014"
    inputs = load_inputs(manifest)
    assert inputs.roles == {"Breast1": "training", "Breast2": "scoring"}
    assert all(s in inputs.cohorts["Breast1"].genes for s in inputs.seeds.seeds)
    assert len(inputs.seeds.seeds) == 6
    assert inputs.qc["Breast1"] == [] and inputs.qc["Breast2"] == []
    assert SEARCH2014.growth.enabled
