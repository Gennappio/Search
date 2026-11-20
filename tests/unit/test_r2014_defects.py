"""The verbatim port of SEARCH 0.0.2: the rules of R it rests on, and the defects of the package
it keeps, each on an example small enough to work out by hand.

These tests say what the port is expected to do and why. That it does what R does is the job
of ``test_r2014_golden.py``, which compares it with answers recorded from R.
"""

from __future__ import annotations

import random
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from search.r2014 import (
    NA_REAL,
    SEARCH,
    NamedVector,
    NotReproducible,
    RError,
    RList,
    RMatrix,
    RRandom,
    SEARCH_randomise,
    getOmax,
    nodes_correlation,
    nodes_glm,
    platform,
    rformat,
    rio,
    rsem,
    rstats,
    summarise_iterations,
    x87,
)
from search.r2014.summarise import _row_means

GENES = ["g1", "g2", "g3", "g4", "g5"]
SAMPLES = [f"s{i + 1}" for i in range(6)]
VALUES = [
    [5.0, 5.0, 5.0, 5.0, 5.0, 5.01],  # nearly flat: its standard deviation is 0.004
    [9.3, 8.95, 7.3, 6.73, 7.38, 8.04],
    [5.67, 7.78, 6.75, 7.27, 7.46, 7.68],
    [8.41, 9.04, 7.87, 9.37, 7.33, 8.35],
    [8.9, 8.09, 7.26, 7.08, 7.54, 8.22],
]


def matrix(values: list[list[float]] = VALUES) -> RMatrix:
    return RMatrix(np.array(values), list(GENES), list(SAMPLES))


# ---- the rules of R ----------------------------------------------------------------------------
def test_assigning_to_a_name_that_is_not_there_appends_it() -> None:
    vec = NamedVector(np.array([1.0, 2.0]), ["a", "b"])
    out = rsem.set_by_names(
        vec, [None, None, "q", "q", "a"], np.array([10.0, 11.0, 12.0, 13.0, 14.0])
    )
    # every missing name is an element of its own; a new name given twice is one element
    assert out.names == ["a", "b", None, None, "q"]
    assert out.values.tolist() == [14.0, 2.0, 10.0, 11.0, 13.0]


def test_reading_a_name_that_is_not_there_gives_a_missing_value_and_name() -> None:
    out = rsem.get_by_names(NamedVector(np.array([1.0, 2.0]), ["a", "b"]), ["a", "zz", None])
    assert out.names == ["a", None, None]
    assert out.values[0] == 1.0 and rsem.is_na_real(out.values).tolist() == [False, True, True]


def test_a_missing_condition_selects_a_missing_element() -> None:
    vec = NamedVector(np.array([0.9, np.nan, 0.2]), ["a", "b", "c"])
    out = rsem.subset_logical(vec, rsem.compare(np.abs(vec.values), ">", 0.5))
    assert out.names == ["a", None] and bool(rsem.is_na_real(out.values)[1])


def test_more_names_than_elements_is_the_error_r_gives() -> None:
    with pytest.raises(
        RError, match=r"'names' attribute \[3\] must be the same length as the vector \[2\]"
    ):
        rsem.set_names(NamedVector(np.array([1.0, 2.0])), ["a", "b", "c"])
    short = rsem.set_names(NamedVector(np.array([1.0, 2.0, 3.0])), ["a"])
    assert short is not None and short.names == ["a", None, None]


def test_numbers_are_written_as_write_table_writes_them() -> None:
    written = [rformat.format_real(v) for v in (0.0004, 0.00047, 1e5, 123456.0, 1 / 3, 1e-05)]
    assert written == ["4e-04", "0.00047", "1e+05", "123456", "0.333333333333333", "1e-05"]
    # as.character() tells NaN from NA; write.table writes both as NA
    assert rformat.format_real(NA_REAL) == "NA" and rformat.format_real(float("nan")) == "NaN"
    table = rio.table_of_vector(NamedVector(np.array([NA_REAL, float("nan")]), ["a", "b"]))
    assert table == "\tx\na\tNA\nb\tNA\n"


def test_a_number_of_17_digits_is_read_as_r_reads_it_not_as_the_nearest_double() -> None:
    text = "9.118641810491681"
    assert float(text) == 9.118641810491681
    assert rformat.r_strtod(text) == 9.11864181049168  # R 4.4.3 on arm64: one unit below
    with platform.use(platform.R_3_1_X86):
        assert rformat.r_strtod(text) == 9.118641810491681  # 80-bit digits: the nearest


def test_sample_draws_as_r_before_and_after_3_6() -> None:
    # set.seed(1); sample(1:10, 5)
    assert RRandom(1, "Rounding").sample(10, 5).tolist() == [3, 4, 5, 7, 2]
    assert RRandom(1, "Rejection").sample(10, 5).tolist() == [9, 4, 7, 1, 2]
    assert (
        platform.R_3_1_X86.sample_kind == "Rounding"
        and platform.R_4_ARM64.sample_kind == "Rejection"
    )


def test_round_changed_in_r_4() -> None:
    value = np.array([0.003125])
    assert rstats.r_round(value, 5)[0] == 0.00313
    with platform.use(platform.R_3_1_X86):
        assert rstats.r_round(value, 5)[0] == 0.00312


def test_the_mean_of_ten_tenths_depends_on_the_processor() -> None:
    tenths = np.full((1, 10), 0.1)
    assert _row_means(tenths)[0] == 0.09999999999999999  # doubles, rounded at every addition
    with platform.use(platform.R_3_1_X86):
        assert _row_means(tenths)[0] == 0.1  # 80 bits, rounded once


def test_extended_arithmetic_rounds_once_to_64_bits() -> None:
    def exact(v: x87.Ext) -> Fraction:
        return Fraction(v[0]) * Fraction(2) ** v[1]

    def rounded(value: Fraction) -> Fraction:
        if value == 0:
            return value
        sign = -1 if value < 0 else 1
        magnitude = abs(value)
        scale = 64 - (magnitude.numerator.bit_length() - magnitude.denominator.bit_length())
        for s in (scale - 1, scale, scale + 1):
            scaled = magnitude * Fraction(2) ** s
            whole = scaled.numerator // scaled.denominator
            if whole.bit_length() == 64:
                rest = scaled - whole
                if rest > Fraction(1, 2) or (rest == Fraction(1, 2) and whole & 1):
                    whole += 1
                return sign * Fraction(whole) / Fraction(2) ** s
        raise AssertionError

    rng = random.Random(1)
    for _ in range(2000):
        a = x87.mul(x87.from_double(rng.uniform(-1e3, 1e3)), x87.from_double(rng.random() + 0.5))
        b = x87.from_double(rng.choice([rng.uniform(-10, 10), float(rng.randint(1, 100)), 0.1]))
        assert exact(x87.add(a, b)) == rounded(exact(a) + exact(b))
        assert exact(x87.mul(a, b)) == rounded(exact(a) * exact(b))
        assert exact(x87.div(a, b)) == rounded(exact(a) / exact(b))
    assert x87.to_double(x87.sqrt(x87.from_double(2.0))) == 2.0**0.5


# ---- the defects of the package ----------------------------------------------------------------
def test_r1_the_sd_prefilter_removes_the_gene_that_carries_the_name_not_the_flat_one() -> None:
    # The standard deviations are computed in matrix order (g1 ... g5) and named in the order of
    # the universe with the seeds first (g3, g4, g5, g1, g2): the small one, g1's, is filed
    # under "g3". The seed g3 is removed and the flat g1 stays.
    result = SEARCH(matrix(), ["g3", "g4", "g5"], GENES, sd_threshold=0.1, max_iter=1)
    assert result.D["1"].names == ["g4", "g5", "g1", "g2"]


def test_r2_omax_compares_a_list_only_with_those_after_it() -> None:
    same = [NamedVector(np.array([1.0, 1.0]), ["a", "b"]) for _ in range(3)]
    lists = RList(same, ["s1", "s2", "s3"])
    omax = getOmax(lists, lists)
    # three identical neighbourhoods overlap completely, yet the last one gets 0
    assert omax is not None and omax.values.tolist() == [1.0, 1.0, 0.0]


def test_r7_a_missing_correlation_puts_an_entry_named_na_in_the_neighbourhood() -> None:
    values = [row[:] for row in VALUES]
    values[1][2] = float("nan")
    lists = nodes_correlation(["g3"], matrix(values), GENES).CL
    assert lists.items[0].names is not None and None in lists.items[0].names
    # and the run stops as soon as the matrices are built
    with pytest.raises(RError, match="subscript out of bounds"):
        SEARCH(matrix(values), ["g3", "g4"], GENES, D_higher=-1.0, max_iter=2)


def test_r5_r6_regression_neighbourhoods_of_one_and_of_none() -> None:
    # only the seed itself passes: one row and one column of R's table, the name is dropped
    alone = nodes_glm(["g3"], matrix(), GENES, significance_threshold=1e-12).CL
    assert alone.names == ["g3"] and alone.items[0].names is None and len(alone.items[0]) == 1
    # nothing passes: the seed is not in the result at all, and SEARCH stops
    assert len(nodes_glm(["g3"], matrix(), GENES, significance_threshold=0.0).CL) == 0
    with pytest.raises(RError, match="subscript out of bounds"):
        SEARCH(
            matrix(),
            ["g3"],
            GENES,
            nodes_association_function="nodes.glm",
            significance_threshold=0.0,
        )


def test_r11_a_constant_seed_in_a_regression_is_flagged() -> None:
    values = [[8.7] * 6, *VALUES[1:]]
    with pytest.warns(NotReproducible, match="same value in every sample"):
        nodes_glm(["g1"], matrix(values), GENES)


def test_one_seed_divides_by_zero() -> None:
    d = SEARCH(matrix(), ["g3"], GENES, max_iter=1).D["1"]
    assert d.names == ["g3", "g1", "g2", "g4", "g5"]
    assert np.isnan(d.values[0]) and d.values[1:].tolist() == [0.0, 0.0, 0.0, 0.0]


def _connectivity_file(folder: Path, replicate: int, values: dict[str, float]) -> None:
    (folder / "output").mkdir(parents=True, exist_ok=True)
    table = rio.table_of_vector(NamedVector(np.array(list(values.values())), list(values)))
    rio.write_and_gzip(table, f"{folder}/output/connectivity__iteration_{replicate}_1.txt")


def test_r9_replicates_after_a_missing_one_are_dropped_from_the_summary(tmp_path: Path) -> None:
    # replicates 1 and 3 have a file, replicate 2 has none
    _connectivity_file(tmp_path / "real", 1, {"a": 0.5, "b": 0.2})
    _connectivity_file(tmp_path / "real", 3, {"a": 0.9, "b": 0.4})
    for batch in (1, 2):
        _connectivity_file(
            tmp_path / "null" / f"batch_{batch}_0_FALSE_0_0", 1, {"a": 0.1, "b": 0.2}
        )
    summary = summarise_iterations(
        f"{tmp_path}/real/", max_iter=1, n_boot=3, output_directory_NULL=f"{tmp_path}/null/",
        n_boot_NULL=1, batches_NULL=2,
    )  # fmt: skip
    real = summary.real["1"]
    assert real.colnames == ["1", "2", "mean", "P"]
    # two files, so two columns, numbered 1 and 2; column 2 is looked up as replicate 2,
    # which does not exist, and replicate 3 is never used: the mean is replicate 1 alone
    assert real.rownames == ["a", "b"]
    assert real.values[:, 0].tolist() == [0.5, 0.2] and np.isnan(real.values[:, 1]).all()
    assert real.values[:, 2].tolist() == [0.5, 0.2]
    # p = (null values >= mean, + 1) / (null columns + 1): a 1/3, b 3/3
    assert real.values[:, 3].tolist() == [1 / 3, 1.0]


def test_summarise_stops_without_null_runs_as_r_does(tmp_path: Path) -> None:
    _connectivity_file(tmp_path / "real", 1, {"a": 0.5, "b": 0.2})
    with pytest.raises(RError, match="argument of length 0"):
        summarise_iterations(f"{tmp_path}/real/", max_iter=1, n_boot=1)  # batches.NULL = NULL
    with pytest.raises(RError, match="dimnames"):  # no replicate reached iteration 2
        summarise_iterations(
            f"{tmp_path}/real/", max_iter=2, n_boot=1, n_boot_NULL=1, batches_NULL=1
        )


def test_r10_a_blank_in_the_output_path_leaves_the_files_uncompressed(tmp_path: Path) -> None:
    def run(folder: str) -> list[str]:
        SEARCH_randomise(
            matrix(), ["g3", "g4"], GENES, sd_threshold=0.0, n_boot=1, sampling_seed=1,
            output_directory=f"{tmp_path}/{folder}",
        )  # fmt: skip
        return sorted(p.suffix for p in (tmp_path / folder / "output").iterdir())

    assert set(run("plain")) == {".gz"}
    assert set(run("with blank")) == {".txt"}  # gzip was handed two halves of the path


def test_n_boot_zero_runs_replicates_one_and_zero(tmp_path: Path) -> None:
    summaries = SEARCH_randomise(
        matrix(), ["g3", "g4"], GENES, sd_threshold=0.0, n_boot=0, sampling_seed=1,
        output_directory=str(tmp_path),
    )  # fmt: skip
    written = {p.name for p in (tmp_path / "output").iterdir() if p.name.startswith("connectivity")}
    assert len(summaries) == 2
    assert written == {"connectivity__iteration_1_1.txt.gz", "connectivity__iteration_0_1.txt.gz"}
