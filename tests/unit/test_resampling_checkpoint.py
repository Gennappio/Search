from pathlib import Path

import numpy as np
import pandas as pd

from search.checkpoint import ChunkMeta, ChunkStore, chunk_ranges
from search.resampling import (
    bootstrap_indices,
    derive_rng,
    eligible_null_genes,
    jackknife_indices,
    random_seed_set,
    resample_indices,
)
from tests.synthetic_data import make_planted


def test_resampling_is_reproducible_and_well_formed() -> None:
    a = resample_indices("bootstrap", 20, 100, 3)
    b = resample_indices("bootstrap", 20, 100, 3)
    c = resample_indices("bootstrap", 20, 100, 4)
    np.testing.assert_array_equal(a, b)
    assert not np.array_equal(a, c)
    assert len(a) == 20 and a.min() >= 0 and a.max() < 20
    j = jackknife_indices(20, derive_rng(5, 1))
    assert len(j) == 19 and len(set(j.tolist())) == 19
    assert len(bootstrap_indices(7, derive_rng(0, 1))) == 7


def test_eligible_genes_and_random_seed_sets() -> None:
    data = make_planted(n_genes=30, n_samples=10, module_size=5, n_seeds=2, seed=1)
    data.expr.values[7] = 1.0  # constant gene
    eligible = eligible_null_genes(data.expr, data.seeds)
    assert "g0007" not in eligible and not set(data.seeds) & set(eligible)
    assert len(eligible) == 30 - 1 - 2
    picked = random_seed_set(eligible, 4, derive_rng(9, 1))
    assert len(set(picked)) == 4 and set(picked) <= set(eligible)
    assert picked == random_seed_set(eligible, 4, derive_rng(9, 1))


def test_chunk_ranges() -> None:
    assert chunk_ranges(250, 100) == [(1, 100), (101, 200), (201, 250)]
    assert chunk_ranges(0, 100) == []
    assert chunk_ranges(3, 0) == [(1, 1), (2, 2), (3, 3)]


def test_chunk_store_round_trip_and_validity(tmp_path: Path) -> None:
    store = ChunkStore(tmp_path / "ck", "bootstrap", "float32")
    genes = pd.Index(["a", "b", "c"])
    meta = ChunkMeta("bootstrap", 1, 2, 10, "cfg", "inp", {})
    store.write(
        1,
        genes,
        {1: np.array([0.1, 0.2, 0.3]), 2: np.array([1.0, np.nan, 3.0])},
        meta,
        {"indices": pd.DataFrame({"rep_00001": [0, 1]})},
    )
    assert store.is_valid(1, meta)
    assert not store.is_valid(1, ChunkMeta("bootstrap", 1, 2, 10, "other", "inp", {}))
    assert not store.is_valid(2, meta)
    frame = store.read(1)
    assert list(frame.columns) == ["rep_00001", "rep_00002"] and list(frame.index) == [
        "a",
        "b",
        "c",
    ]
    np.testing.assert_allclose(frame["rep_00001"].to_numpy(), [0.1, 0.2, 0.3], rtol=1e-6)
    assert np.isnan(frame.loc["b", "rep_00002"])
    assert list(store.read_aux(1, "indices")["rep_00001"]) == [0, 1]
    store.delete(1)
    assert not store.is_valid(1, meta)
