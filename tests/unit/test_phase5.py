"""Comparison of runs, variant grids, plots and the h5ad round trip."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from search.cli import app
from search.compare import compare_runs, summarise_variants, variant_label, variants_grid
from search.data import ExpressionSet, read_matrix, write_matrix
from search.export.tables import read_cohort_tables
from search.run import CohortManifest, run_pipeline
from tests.unit.test_pipeline import write_manifest_files


def test_grid_and_labels() -> None:
    combos = variants_grid(
        {"membership.mode": ["binary", "cluster_size"], "null.pooling": ["pooled", "per_gene"]}
    )
    assert len(combos) == 4 and combos[0] == {"membership.mode": "binary", "null.pooling": "pooled"}
    assert variant_label(combos[-1]) == "mode-cluster_size__pooling-per_gene"
    assert variant_label({}) == "baseline"


def test_compare_runs_and_variants(tmp_path: Path) -> None:
    path = write_manifest_files(tmp_path / "proj", n_boot=2, n_null=3)
    manifest = CohortManifest.load(path)
    base = run_pipeline(manifest, manifest.build_config(), tmp_path / "runs")
    other_cfg = manifest.build_config(overrides={"membership.mode": "binary"})
    other = run_pipeline(
        manifest.model_copy(update={"name": "binary"}), other_cfg, tmp_path / "runs"
    )
    same = compare_runs(base.run_dir, base.run_dir)
    assert same.signature_jaccard == pytest.approx(1.0) and not same.only_a and not same.only_b
    assert (same.cohort_metrics["spearman"].dropna() > 0.999).all()
    diff = compare_runs(base.run_dir, other.run_dir, top_n=10)
    assert {"cohort", "column", "pearson", "spearman", "top10_overlap"} <= set(
        diff.cohort_metrics.columns
    )
    assert "signature Jaccard" in diff.summary()
    table = summarise_variants(
        {"baseline": base.run_dir, "binary": other.run_dir}, "baseline", top_n=10
    )
    assert list(table["variant"]) == ["baseline", "binary"]
    assert table.loc[0, "spearman_vs_baseline"] == pytest.approx(1.0)
    assert table["top10"].str.len().gt(0).all()


def test_saved_results_distinguish_padded_gene_ids(tmp_path: Path) -> None:
    for name, candidate in [("a", "007"), ("b", "7")]:
        results = tmp_path / name / "results"
        results.mkdir(parents=True)
        (results / "001_gene_connectivity.csv").write_text(
            "dataset,gene,C_raw\n001,007,0.5\n001,7,NA\n"
        )
        (results / "combined_full_ranking.csv").write_text(
            f"gene,candidate\n007,{candidate == '007'}\n7,{candidate == '7'}\n"
        )
    table = read_cohort_tables(tmp_path / "a" / "results")["001"]
    assert table["gene"].tolist() == ["007", "7"]
    np.testing.assert_equal(table["C_raw"].to_numpy(), [0.5, np.nan])
    comparison = compare_runs(tmp_path / "a", tmp_path / "b")
    assert comparison.only_a == ["007"] and comparison.only_b == ["7"]
    assert comparison.signature_jaccard == 0.0


def test_cli_compare_variants_and_plot(tmp_path: Path) -> None:
    pytest.importorskip("matplotlib")
    path = write_manifest_files(tmp_path / "proj", n_boot=2, n_null=3)
    runner = CliRunner()
    out = tmp_path / "variants"
    res = runner.invoke(
        app,
        [
            "variants",
            "--manifest",
            str(path),
            "--grid",
            "membership.mode=cluster_size,binary",
            "--out",
            str(out),
            "--top-n",
            "10",
        ],
    )
    assert res.exit_code == 0, res.output
    summary = pd.read_csv(out / "variants_summary.csv")
    assert list(summary["variant"]) == ["mode-cluster_size", "mode-binary"]
    runs = [Path(p) for p in summary["run_dir"]]
    res = runner.invoke(
        app, ["compare", str(runs[0]), str(runs[1]), "--out", str(tmp_path / "cmp.csv")]
    )
    assert res.exit_code == 0 and (tmp_path / "cmp.csv").exists()
    res = runner.invoke(app, ["plot", str(runs[0])])
    assert res.exit_code == 0, res.output
    plots = sorted((runs[0] / "plots").glob("*.png"))
    names = {p.name for p in plots}
    assert "rank_scatter.png" in names and "product_histogram.png" in names
    assert any(n.startswith("A_seed_network") for n in names)


def test_h5ad_round_trip(tmp_path: Path) -> None:
    pytest.importorskip("anndata")
    values = np.arange(12, dtype=float).reshape(3, 4)
    table = pd.DataFrame({"batch": ["a", "a", "b", "b"]}, index=["s1", "s2", "s3", "s4"])
    expr = ExpressionSet(
        values, pd.Index(["g1", "g2", "g3"]), pd.Index(table.index), sample_table=table
    )
    path = tmp_path / "m.h5ad"
    write_matrix(expr, path)
    back = read_matrix(path)
    np.testing.assert_allclose(back.values, values)
    assert list(back.genes) == ["g1", "g2", "g3"] and list(back.samples) == list(table.index)
    assert back.sample_table is not None and list(back.sample_table["batch"]) == [
        "a",
        "a",
        "b",
        "b",
    ]
