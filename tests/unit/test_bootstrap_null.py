from pathlib import Path

import numpy as np
import pandas as pd

from search.bootstrap import run_bootstrap
from search.core.corepass import run_core
from search.null import empirical_p, per_gene_counts, pooled_counts, run_internal_null
from search.presets import BUFFA2010
from tests.synthetic_data import make_planted


def cfg(n_boot: int, n_null: int, chunk: int = 4, **extra):  # type: ignore[no-untyped-def]
    overrides = {
        "resampling.n_rep": n_boot,
        "resampling.chunk_size": chunk,
        "null.n_rep": n_null,
        "null.chunk_size": chunk,
        **extra,
    }
    return BUFFA2010.with_overrides(overrides)


def test_bootstrap_summary_and_resume(tmp_path: Path) -> None:
    data = make_planted(n_genes=40, n_samples=24, module_size=10, n_seeds=3, seed=21)
    c = cfg(6, 0)
    summary = run_bootstrap(
        data.expr, data.seeds, c, universe=None, ckpt_dir=tmp_path, input_hash="x"
    )
    assert summary.n_rep == 6 and summary.statistic == "C_rank"
    table = summary.table()
    assert list(table.columns) == [
        "gene",
        "C_rank_boot_mean",
        "C_rank_boot_var",
        "C_rank_boot_sd",
        "boot_n_ok",
    ]
    assert table["boot_n_ok"].eq(6).all()
    assert table["C_rank_boot_mean"].between(0, 1).all()
    # resume: delete one chunk, rerun, the other chunk is untouched and the result identical
    first = tmp_path / "bootstrap_chunk_0001.parquet"
    second = tmp_path / "bootstrap_chunk_0002.parquet"
    mtime_second = second.stat().st_mtime_ns
    first.unlink()
    again = run_bootstrap(
        data.expr, data.seeds, c, universe=None, ckpt_dir=tmp_path, input_hash="x"
    )
    assert second.stat().st_mtime_ns == mtime_second
    np.testing.assert_allclose(again.mean, summary.mean)
    # a different input hash invalidates everything
    changed = run_bootstrap(
        data.expr, data.seeds, c, universe=None, ckpt_dir=tmp_path, input_hash="y"
    )
    assert second.stat().st_mtime_ns != mtime_second
    np.testing.assert_allclose(changed.mean, summary.mean)
    indices = pd.read_parquet(tmp_path / "bootstrap_indices_0001.parquet")
    assert indices.shape == (24, 4)


def test_bootstrap_matches_direct_replicates(tmp_path: Path) -> None:
    data = make_planted(n_genes=30, n_samples=16, module_size=8, n_seeds=2, seed=22)
    c = cfg(3, 0, chunk=10)
    summary = run_bootstrap(
        data.expr, data.seeds, c, universe=None, ckpt_dir=tmp_path, input_hash="x"
    )
    chunk = pd.read_parquet(tmp_path / "bootstrap_chunk_0001.parquet").set_index("gene")
    idx = pd.read_parquet(tmp_path / "bootstrap_indices_0001.parquet")
    for rep in (1, 2, 3):
        direct = run_core(data.expr, data.seeds, c, sample_idx=idx[f"rep_{rep:05d}"].to_numpy())
        assert direct.C_rank is not None
        np.testing.assert_allclose(chunk[f"rep_{rep:05d}"].to_numpy(), direct.C_rank)
    np.testing.assert_allclose(summary.mean, chunk.to_numpy().mean(axis=1))


def test_null_counting_helpers() -> None:
    observed = np.array([0.5, 0.1, np.nan])
    chunk = np.array([[0.4, 0.6, np.nan], [0.1, 0.0, 0.2], [0.9, 0.9, 0.9]])
    counts, n = pooled_counts(observed, chunk)
    assert n == 8
    np.testing.assert_array_equal(counts, [4, 7, 0])  # values >= 0.5: 0.6, 0.9 x3 ; >= 0.1: 7 of 8
    counts_g, n_g = per_gene_counts(observed, chunk)
    np.testing.assert_array_equal(counts_g, [1, 2, 0])
    np.testing.assert_array_equal(n_g, [2, 3, 3])
    p = empirical_p(np.array([0, 3, 0]), np.array([9, 9, 0]))
    np.testing.assert_allclose(p[:2], [0.1, 0.4])
    assert np.isnan(p[2])


def test_internal_null_pooled_and_per_gene(tmp_path: Path) -> None:
    data = make_planted(n_genes=40, n_samples=40, module_size=10, n_seeds=3, strength=2.0, seed=23)
    observed = run_core(data.expr, data.seeds, BUFFA2010).C_raw
    pooled = run_internal_null(
        data.expr,
        data.seeds,
        cfg(0, 8),
        universe=None,
        observed=observed,
        ckpt_dir=tmp_path / "a",
        input_hash="x",
    )
    assert pooled.pooling == "pooled" and (pooled.n_null == 8 * 40).all()
    assert pooled.p.min() > 0 and pooled.p.max() <= 1
    assert len(pooled.random_seed_summary) == 8 and "seeds" in pooled.random_seed_summary.columns
    # module genes should look less random than background genes
    table = pooled.table().set_index("gene")
    module = [g for g in data.module if g not in data.seeds]
    assert table.loc[module, "p_MC_C"].median() < table.loc[data.background, "p_MC_C"].median()
    per_gene = run_internal_null(
        data.expr,
        data.seeds,
        cfg(0, 8, **{"null.pooling": "per_gene"}),
        universe=None,
        observed=observed,
        ckpt_dir=tmp_path / "b",
        input_hash="x",
    )
    assert (per_gene.n_null == 8).all()
    # the same chunks serve both poolings: identical null values on disk
    a = pd.read_parquet(tmp_path / "a" / "null_chunk_0001.parquet")
    b = pd.read_parquet(tmp_path / "b" / "null_chunk_0001.parquet")
    pd.testing.assert_frame_equal(a, b)
