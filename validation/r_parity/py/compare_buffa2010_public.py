"""Parity harness, Python side: the port (``buffa2010`` preset) against the functions of the
lab's GBM R script on public data (Breastdata1), so that the comparison can be repeated by
anyone who has the repository. See ``r/gbm_functions_public.R`` for the R side. Usage::

    python validation/r_parity/py/compare_buffa2010_public.py <repository root> <output dir>
"""

from __future__ import annotations

import shutil
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from common import Row, print_report, read_lines, write_report  # noqa: E402
from compare_buffa2010 import RDraws, compare_cohort, compare_core, r_draws  # noqa: E402

import search.run as search_run  # noqa: E402
from search.core.corepass import run_core  # noqa: E402
from search.data import ExpressionSet, read_matrix  # noqa: E402
from search.presets import get_preset  # noqa: E402

DATASET = "Breast1"


def main() -> None:
    warnings.filterwarnings("ignore")
    repo = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve()
    export = out / "r_export"
    info = dict(line.split("\t", 1) for line in read_lines(export / "run_info.tsv"))
    b_boot, n_mc = int(info["b_boot"]), int(info["n_mc"])
    seeds = read_lines(export / "seeds.txt")
    expr = read_matrix(repo / "tests" / "data" / "breast" / "Breastdata1" / "mRNA_abundance.txt.gz")
    cfg = get_preset("buffa2010").with_overrides(
        {
            "resampling.n_rep": b_boot,
            "null.n_rep": n_mc,
            "resampling.chunk_size": 5,
            "null.chunk_size": 5,
        }
    )
    run_dir = out / "python_same_draws"
    shutil.rmtree(run_dir, ignore_errors=True)
    draws = RDraws(export, (DATASET,))
    with r_draws(draws):
        # through the module, so that the wrapper installed by r_draws sees the cohort name
        run = search_run.run_cohort(expr, seeds, cfg, DATASET, run_dir, resume=False, n_jobs=1)
    rows: list[Row] = compare_cohort(
        DATASET, "same draws", export, run, run_dir, cfg, draws, b_boot, n_mc, expr
    )

    # missing values: the script correlates over the samples that are present, search gives
    # no association to a gene with a missing value
    na_rows = read_lines(export / "na_rows.txt")
    values = expr.values.copy()
    values[expr.gene_positions(na_rows), 2] = np.nan
    with_na = ExpressionSet(values, expr.genes, expr.samples)
    core = run_core(with_na, seeds, cfg)
    label = f"{DATASET}, one missing value in {len(na_rows)} probesets"
    rows.extend(compare_core(label, "by design", export, core, r_name="Breast1NA"))
    frame = write_report(rows, out / "report_buffa2010_public.csv")
    print_report(frame)
    print("\n", frame.groupby(["variant", "status"]).size().to_string())


if __name__ == "__main__":
    main()
