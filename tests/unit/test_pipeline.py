"""End-to-end run on synthetic cohorts, through the API and the command line."""

from pathlib import Path

import pandas as pd
import yaml
from typer.testing import CliRunner

from search.cli import app
from search.data import write_matrix
from search.export.manifest import read_manifest
from search.run import CohortManifest, load_inputs, run_pipeline
from tests.synthetic_data import make_planted


def write_manifest_files(root: Path, n_boot: int = 4, n_null: int = 6) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    cohorts = []
    for name, seed, role in (("A", 31, "training"), ("B", 32, "training"), ("C", 33, "scoring")):
        data = make_planted(
            n_genes=60, n_samples=40, module_size=15, n_seeds=4, strength=2.0, seed=seed
        )
        matrix = root / f"{name}.csv"
        write_matrix(data.expr, matrix)
        cohorts.append({"name": name, "role": role, "matrix": matrix.name})
        seeds = data.seeds
    (root / "seeds.csv").write_text(
        "gene,symbol\n" + "\n".join(f"{g},SYM_{g}" for g in seeds) + "\n"
    )
    manifest = {
        "name": "synthetic",
        "seeds": "seeds.csv",
        "cohorts": cohorts,
        "config": {
            "preset": "buffa2010",
            "overrides": {
                "resampling.n_rep": n_boot,
                "resampling.chunk_size": 2,
                "null.n_rep": n_null,
                "null.chunk_size": 3,
                "combine.mc_p_cutoff": 0.6,
                "combine.product_cutoff": 0.2,
            },
        },
    }
    path = root / "cohorts.yaml"
    path.write_text(yaml.safe_dump(manifest))
    return path


def test_manifest_loading_and_inputs(tmp_path: Path) -> None:
    path = write_manifest_files(tmp_path / "proj")
    manifest = CohortManifest.load(path)
    assert [c.name for c in manifest.cohorts] == ["A", "B", "C"]
    cfg = manifest.build_config()
    assert cfg.preset == "buffa2010" and cfg.resampling.n_rep == 4
    cfg2 = manifest.build_config(overrides={"null.pooling": "per_gene"})
    assert cfg2.null.pooling == "per_gene"
    inputs = load_inputs(manifest)
    assert set(inputs.cohorts) == {"A", "B", "C"} and inputs.roles["C"] == "scoring"
    assert len(inputs.seeds.seeds) == 4 and inputs.symbols[inputs.seeds.seeds[0]].startswith("SYM_")
    assert all("sample_table_missing" in {w.code for w in ws} for ws in inputs.qc.values())


def test_pipeline_end_to_end_and_resume(tmp_path: Path) -> None:
    path = write_manifest_files(tmp_path / "proj")
    manifest = CohortManifest.load(path)
    cfg = manifest.build_config()
    result = run_pipeline(manifest, cfg, tmp_path / "runs")
    run_dir = result.run_dir
    assert (run_dir / "config.yaml").exists() and (run_dir / "log.txt").exists()
    stored = read_manifest(run_dir)
    assert stored["status"] == "finished"
    # the manifest overrides the selection cutoffs, so the run matches no validated preset
    assert stored["preset_matched"] is None and stored["preset_validated"] is False
    assert "matches no validated preset" in (run_dir / "log.txt").read_text()
    assert (
        set(stored["inputs"]) == {"A", "B", "C"} and stored["inputs"]["A"]["n_genes_common"] == 60
    )
    results = run_dir / "results"
    for name in ("A", "B"):
        table = pd.read_csv(results / f"{name}_gene_connectivity.csv")
        assert {"C_rank_boot_mean", "p_MC_C", "positive_consistent"} <= set(table.columns)
        assert len(table) == 60
    assert (results / "combined_full_ranking.csv").exists()
    assert result.combined is not None and len(result.combined.signature) >= 1
    assert "symbol" in result.combined.table.columns
    assert "signature_all" in result.scores
    scores = pd.read_csv(results / "scores_signature_all.csv")
    assert set(scores["dataset"]) == {"A", "B", "C"} and len(scores) == 120
    assert (tmp_path / "runs" / "latest" / "run_path.txt").read_text() == str(run_dir.resolve())
    # resume into the same directory: checkpoints are reused, results identical
    (run_dir / "checkpoints" / "A" / "null_chunk_0002.parquet").unlink()
    again = run_pipeline(manifest, cfg, tmp_path / "runs", resume=True, run_dir=run_dir)
    assert again.run_dir == run_dir
    log = (run_dir / "log.txt").read_text()
    assert "null chunk 1/2 already done" in log and "bootstrap chunk 1/2 already done" in log
    pd.testing.assert_frame_equal(
        again.cohorts["A"].gene_table(), result.cohorts["A"].gene_table(), check_exact=False
    )


def test_cli_run_combine_score_and_config(tmp_path: Path) -> None:
    path = write_manifest_files(tmp_path / "proj", n_boot=2, n_null=3)
    runner = CliRunner()
    out = tmp_path / "runs"
    res = runner.invoke(
        app, ["run", "--manifest", str(path), "--out", str(out), "--set", "null.chunk_size=3"]
    )
    assert res.exit_code == 0, res.output
    assert "run directory" in res.output and "candidate genes" in res.output
    res = runner.invoke(
        app, ["combine", str(out / "latest"), "--product-cutoff", "0.1", "--prefix", "loose"]
    )
    assert res.exit_code == 0, res.output
    run_dir = Path((out / "latest" / "run_path.txt").read_text())
    assert (run_dir / "results" / "loose_signature.csv").exists()
    res = runner.invoke(app, ["qc", "--manifest", str(path)])
    assert res.exit_code == 0 and "sample_table_missing" in res.output
    signature = run_dir / "results" / "combined_signature.csv"
    res = runner.invoke(
        app,
        [
            "score",
            "--manifest",
            str(path),
            "--signature",
            str(signature),
            "--out",
            str(tmp_path / "s.csv"),
        ],
    )
    assert res.exit_code == 0, res.output
    assert len(pd.read_csv(tmp_path / "s.csv")) == 120
    res = runner.invoke(
        app, ["config", "show", "--preset", "search2014", "--set", "growth.max_iter=3"]
    )
    assert (
        res.exit_code == 0 and "max_iter: 3" in res.output and "matches preset: None" in res.output
    )
    (tmp_path / "cfg.yaml").write_text(res.output.split("# matches")[0])
    res = runner.invoke(app, ["config", "validate", str(tmp_path / "cfg.yaml")])
    assert res.exit_code == 0 and "valid" in res.output
    res = runner.invoke(app, ["resume", str(out / "latest")])
    assert res.exit_code == 0, res.output
    res = runner.invoke(
        app,
        [
            "run",
            "--manifest",
            str(path),
            "--out",
            str(out),
            "--strict",
            "--set",
            "membership.mode=binary",
        ],
    )
    assert res.exit_code == 1 and "strict" in res.output


def test_pipeline_with_growth_preset(tmp_path: Path) -> None:
    path = write_manifest_files(tmp_path / "proj", n_boot=2, n_null=2)
    manifest = CohortManifest.load(path)
    cfg = manifest.build_config(
        preset="search2014",
        overrides={
            "growth.max_iter": 3,
            "growth.D_higher": 0.5,
            "growth.D_lower": 0.1,
            "growth.sd_threshold": 0.0,
            "resampling.n_rep": 2,
            "resampling.chunk_size": 2,
            "null.scheme": "internal_random_seeds",
            "null.n_rep": 2,
            "null.chunk_size": 2,
            "combine.enabled": True,
            "combine.mc_p_cutoff": 0.9,
            "combine.product_cutoff": 0.0,
        },
    )
    result = run_pipeline(manifest, cfg, tmp_path / "runs")
    run = result.cohorts["A"]
    assert run.growth is not None and run.growth.n_iterations >= 1
    results = result.run_dir / "results"
    assert (results / "A_growth_iterations.csv").exists() and (
        results / "A_growth_events.csv"
    ).exists()
    table = pd.read_csv(results / "A_gene_connectivity.csv")
    assert "D_final" in table.columns and "C_raw_boot_mean" in table.columns
