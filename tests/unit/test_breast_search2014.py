"""The search2014 preset on the breast matrix bundled with the original package."""

from pathlib import Path

import pandas as pd

from search.data import read_matrix
from search.growth import run_growth
from search.presets import SEARCH2014


def test_search2014_grows_a_network_on_breastdata1(breast_dir: Path) -> None:
    expr = read_matrix(breast_dir / "Breastdata1" / "mRNA_abundance.txt.gz")
    seeds = [str(g) for g in expr.genes[:10]]  # the seeds of the original vignette
    cfg = SEARCH2014.with_overrides(
        {"growth.max_iter": 3, "growth.D_lower": 0.0, "growth.D_higher": 0.0, "growth.O_lower": 0.1}
    )
    result = run_growth(expr, seeds, cfg)
    assert result.n_iterations >= 1
    assert result.stop_reason in {"converged", "max_iter", "fewer_than_two_seeds"}
    table = result.iteration_table()
    assert len(table) == result.n_iterations
    events = result.event_table()
    assert set(events["event"]) <= {"added", "removed_O", "removed_D"}
    assert len(result.final_seeds) >= 2
    d = result.iterations[0].D
    # the preset's sd threshold removes flat genes before the loop
    assert d.between(0, 1).all() and len(d) == len(result.final_universe) <= expr.n_genes
    assert any(r["reason"] == "low_sd" for r in result.removed)
    assert isinstance(result.connectivity_table(), pd.DataFrame)
