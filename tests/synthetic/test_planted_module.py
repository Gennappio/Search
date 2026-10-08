"""A planted co-expression module must be recovered from a few of its genes used as seeds."""

import numpy as np
import pytest

from search.core.corepass import run_core
from search.presets import BUFFA2010, SEARCH2014
from tests.synthetic_data import auc, make_planted


@pytest.mark.parametrize("preset", [BUFFA2010, SEARCH2014], ids=["buffa2010", "search2014"])
def test_module_genes_rank_above_background(preset) -> None:
    data = make_planted(n_genes=300, n_samples=60, module_size=50, n_seeds=5, strength=1.5, seed=11)
    res = run_core(data.expr, data.seeds, preset)
    table = res.gene_table().set_index("gene")
    score = table["C_raw"]
    module = score.loc[[g for g in data.module if g not in data.seeds]].to_numpy()
    background = score.loc[data.background].to_numpy()
    assert auc(module, background) > 0.95
    assert np.isfinite(res.S).all() and np.nanmean(res.S[np.triu_indices(5, 1)]) > 0.3


def test_unrelated_second_module_does_not_overlap() -> None:
    data = make_planted(
        n_genes=300,
        n_samples=60,
        module_size=50,
        n_seeds=4,
        strength=1.5,
        seed=12,
        second_module=40,
    )
    foreign = data.expr.genes[50]  # first gene of the second module
    res = run_core(data.expr, data.seeds + [str(foreign)], BUFFA2010)
    o = res.O
    among_true = o[:4, :4][np.triu_indices(4, 1)]
    with_foreign = o[4, :4]
    assert among_true.mean() > 0.3 and with_foreign.max() < 0.1
