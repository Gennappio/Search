"""Reasonableness on the real GBM cohorts: with the ten seeds of Buffa 2010, the genes of the
published 51-gene metagene must be enriched at the top of the connectivity ranking."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats as sps

from search.core.corepass import run_core
from search.data import read_matrix, read_seeds, read_table
from search.presets import BUFFA2010

EXAMPLES = Path(__file__).resolve().parents[2] / "examples" / "gbm"

pytestmark = pytest.mark.realdata


@pytest.fixture(scope="module")
def metagene_symbols() -> set[str]:
    table = pd.read_csv(EXAMPLES / "buffa2010_metagene_51.csv")
    symbols = set(table["symbol"].astype(str)) | set(table["alias"].dropna().astype(str))
    return symbols


@pytest.mark.parametrize("cohort", ["IT", "FR"])
def test_metagene_enriched_in_top_connectivity(
    gbm_dir: Path, cohort: str, metagene_symbols: set[str]
) -> None:
    suffix = "_Deseq2_vst_Combat.csv"
    expr = read_matrix(
        gbm_dir / "DATA" / "NORMALIZED" / f"{cohort}_GBM_RNA{suffix}", strip_version=True
    )
    annotation = read_table(
        gbm_dir / "DATA" / "master_annotation_ITFR_Ensembl104_TCGA_GDC_consensus.csv"
    )
    symbols = annotation["symbol_consensus"].dropna().astype(str).to_dict()
    seeds = read_seeds(EXAMPLES / "seeds_buffa2010_setA.csv").seeds
    result = run_core(expr, seeds, BUFFA2010)
    table = result.gene_table()
    table["symbol"] = table["gene"].map(symbols)
    non_seed = table[~table["is_seed"]].sort_values("C_raw", ascending=False)
    reference = non_seed["symbol"].isin(metagene_symbols)
    top = 200
    hits = int(reference.head(top).sum())
    total_ref = int(reference.sum())
    universe = len(non_seed)
    p = sps.hypergeom.sf(hits - 1, universe, total_ref, top)
    print(
        f"\n{cohort}: {hits}/{total_ref} metagene genes in the top {top} of {universe} (p={p:.2e})"
    )
    assert total_ref >= 30, "the annotation did not map the metagene symbols"
    assert hits >= 8 and p < 1e-4
    # the seeds are strongly interconnected: every seed shares neighbours with the others
    assert result.S is not None and np.nanmean(result.S[np.triu_indices(10, 1)]) > 0.1
