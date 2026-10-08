"""The recursive growth loop on planted structure."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from search.bootstrap import growth_evaluator, run_bootstrap
from search.core.corepass import run_core
from search.growth import run_growth
from search.null import external_null
from search.presets import SEARCH2014
from tests.synthetic_data import make_planted


def growth_cfg(**overrides):  # type: ignore[no-untyped-def]
    base = {
        "growth.enabled": True,
        "growth.max_iter": 4,
        "growth.D_lower": 0.1,
        "growth.D_higher": 0.5,
        "growth.O_lower": 0.05,
        "growth.keep_seeds": True,
        "growth.sd_threshold": 0.0,
    }
    base.update(overrides)
    return SEARCH2014.with_overrides(base)


def test_growth_recruits_the_planted_module() -> None:
    data = make_planted(n_genes=200, n_samples=60, module_size=30, n_seeds=5, strength=2.0, seed=41)
    result = run_growth(data.expr, data.seeds, growth_cfg())
    final = set(result.final_seeds)
    module = set(data.module)
    assert len(final & module) >= 25, "most module genes should be recruited as seeds"
    assert len(final - module) <= 3, "background genes should not be recruited"
    assert result.stop_reason in {"converged", "max_iter"}
    assert result.iterations[0].added and not result.iterations[0].removed_O
    table = result.iteration_table()
    assert list(table["iteration"]) == list(range(1, result.n_iterations + 1))
    events = result.event_table()
    assert (events["event"] == "added").sum() >= 20


def test_foreign_seed_is_pruned_by_D_then_by_omax() -> None:
    data = make_planted(
        n_genes=200,
        n_samples=60,
        module_size=30,
        n_seeds=4,
        strength=2.0,
        seed=42,
        second_module=30,
    )
    foreign = str(data.expr.genes[30])  # first gene of the unrelated module
    seeds = data.seeds + [foreign]
    # D pruning: the foreign seed is in no other seed's neighbourhood -> D = 0 < D_lower
    by_d = run_growth(data.expr, seeds, growth_cfg(**{"growth.max_iter": 1}))
    assert foreign in by_d.iterations[0].removed_D
    assert foreign not in by_d.final_seeds
    assert any(r["gene"] == foreign and r["reason"] == "D_lower" for r in by_d.removed)
    # O pruning: with D_lower = 0 the foreign seed survives iteration 1; at iteration 2 seeds are
    # not protected and its neighbourhood overlaps with nobody -> Omax < O_lower
    by_o = run_growth(
        data.expr,
        seeds,
        growth_cfg(**{"growth.D_lower": 0.0, "growth.keep_seeds": False, "growth.max_iter": 2}),
    )
    assert by_o.n_iterations >= 2
    assert foreign in by_o.iterations[1].removed_O
    assert foreign not in by_o.final_seeds and foreign not in by_o.final_universe


def test_one_iteration_without_change_equals_single_pass() -> None:
    data = make_planted(n_genes=120, n_samples=50, module_size=20, n_seeds=4, strength=2.0, seed=43)
    cfg = growth_cfg(**{"growth.D_lower": 0.0, "growth.D_higher": 2.0, "growth.max_iter": 3})
    result = run_growth(data.expr, data.seeds, cfg)
    assert result.n_iterations == 1 and result.stop_reason == "converged"
    core = run_core(data.expr, data.seeds, cfg)
    np.testing.assert_allclose(result.final_D(), core.C_raw)
    assert result.final_seeds == data.seeds


def test_growth_stops_when_seeds_vanish() -> None:
    data = make_planted(n_genes=80, n_samples=30, module_size=10, n_seeds=3, strength=0.0, seed=44)
    result = run_growth(data.expr, data.seeds, growth_cfg(**{"growth.D_lower": 0.5}))
    assert result.stop_reason in {"fewer_than_two_seeds", "converged"}
    assert len(result.final_seeds) < 3 or result.iterations[0].removed_D == []


def test_sd_prefilter_is_keyed_by_gene_name() -> None:
    data = make_planted(n_genes=60, n_samples=30, module_size=10, n_seeds=3, strength=2.0, seed=45)
    data.expr.values[50] = 1.0  # constant gene
    result = run_growth(
        data.expr, data.seeds, growth_cfg(**{"growth.sd_threshold": 0.01, "growth.max_iter": 1})
    )
    assert "g0050" not in result.final_universe
    assert any(r["gene"] == "g0050" and r["reason"] == "low_sd" for r in result.removed)
    shuffled = data.expr.subset_genes(list(reversed(list(data.expr.genes))))
    again = run_growth(
        shuffled, data.seeds, growth_cfg(**{"growth.sd_threshold": 0.01, "growth.max_iter": 1})
    )
    assert "g0050" not in again.final_universe
    pd.testing.assert_series_equal(
        result.iterations[0].D.sort_index(), again.iterations[0].D.sort_index()
    )


def test_bootstrap_over_growth_and_external_null(tmp_path: Path) -> None:
    data = make_planted(n_genes=120, n_samples=50, module_size=20, n_seeds=4, strength=2.0, seed=46)
    cfg = growth_cfg(
        **{"resampling.n_rep": 4, "resampling.chunk_size": 2, "resampling.statistic": "C_raw"}
    )
    summary = run_bootstrap(
        data.expr, data.seeds, cfg, universe=None, ckpt_dir=tmp_path / "real", input_hash="x"
    )
    assert summary.n_rep == 4 and len(summary.genes) == 120
    assert summary.table()["C_raw_boot_mean"].notna().sum() > 100
    direct = growth_evaluator(data.expr, data.seeds, cfg, None)
    assert direct.shape == (120,)
    # random-seed runs as external null
    rng = np.random.default_rng(0)
    for k in range(2):
        random_seeds = [
            str(g) for g in rng.choice(np.asarray(data.expr.genes[40:]), 4, replace=False)
        ]
        run_bootstrap(
            data.expr,
            random_seeds,
            cfg,
            universe=None,
            ckpt_dir=tmp_path / f"null{k}",
            input_hash="x",
        )
    observed = summary.mean
    result = external_null(
        observed, summary.genes, [tmp_path / "null0", tmp_path / "null1"], pooling="per_gene"
    )
    assert result.n_rep == 8 and (result.n_null <= 8).all()
    table = pd.DataFrame({"gene": summary.genes, "p": result.p}).set_index("gene")
    module = [g for g in data.module if g not in data.seeds]
    assert table.loc[module, "p"].median() < table.loc[data.background, "p"].median()
    with pytest.raises(ValueError):
        external_null(observed, summary.genes, [tmp_path / "nowhere"])
