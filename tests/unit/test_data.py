from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from search.data import (
    ExpressionSet,
    read_gene_list,
    read_matrix,
    read_sample_table,
    read_seeds,
    strip_ensembl_version,
    write_matrix,
)
from search.errors import DataError


def small() -> ExpressionSet:
    values = np.arange(12, dtype=float).reshape(3, 4)
    table = pd.DataFrame({"batch": ["a", "a", "b", "b"]}, index=["s1", "s2", "s3", "s4"])
    return ExpressionSet(
        values, pd.Index(["g1", "g2", "g3"]), pd.Index(["s1", "s2", "s3", "s4"]), None, table
    )


def test_validation_errors() -> None:
    with pytest.raises(DataError, match="duplicated gene"):
        ExpressionSet(np.zeros((2, 2)), pd.Index(["a", "a"]), pd.Index(["x", "y"]))
    with pytest.raises(DataError, match="lacks"):
        ExpressionSet(
            np.zeros((1, 2)),
            pd.Index(["a"]),
            pd.Index(["x", "y"]),
            sample_table=pd.DataFrame(index=["x"]),
        )
    with pytest.raises(DataError):
        ExpressionSet(np.zeros((2, 2)), pd.Index(["a"]), pd.Index(["x", "y"]))


def test_sample_table_is_reordered_to_match_columns() -> None:
    table = pd.DataFrame({"batch": ["b", "a"]}, index=["s2", "s1"])
    expr = ExpressionSet(
        np.zeros((1, 2)), pd.Index(["g"]), pd.Index(["s1", "s2"]), sample_table=table
    )
    assert expr.sample_table is not None
    assert list(expr.sample_table["batch"]) == ["a", "b"]
    assert expr.sample_table_provided


def test_take_samples_keeps_table_aligned_with_repeats() -> None:
    expr = small()
    boot = expr.take_samples([3, 3, 0, 1])
    assert boot.n_samples == 4
    assert list(boot.samples) == ["s4", "s4", "s1", "s2"]
    assert boot.sample_table is not None
    assert list(boot.sample_table["batch"]) == ["b", "b", "a", "a"]
    np.testing.assert_array_equal(boot.values[:, 0], expr.values[:, 3])
    with pytest.raises(DataError):
        expr.take_samples([4])


def test_subset_and_positions_and_common() -> None:
    expr = small()
    sub = expr.subset_genes(["g3", "g1"])
    assert list(sub.genes) == ["g3", "g1"]
    np.testing.assert_array_equal(sub.values[0], expr.values[2])
    np.testing.assert_array_equal(expr.gene_positions(["g2"]), [1])
    with pytest.raises(DataError):
        expr.gene_positions(["nope"])
    other = ExpressionSet(np.zeros((2, 1)), pd.Index(["g3", "zz"]), pd.Index(["q"]))
    assert list(ExpressionSet.common_genes(expr, other)) == ["g3"]


def test_round_trips(tmp_path: Path) -> None:
    expr = small()
    for name in ["m.csv", "m.tsv", "m.parquet", "m.csv.gz"]:
        path = tmp_path / name
        write_matrix(expr, path)
        back = read_matrix(path)
        np.testing.assert_allclose(back.values, expr.values)
        assert list(back.genes) == list(expr.genes) and list(back.samples) == list(expr.samples)
    transposed = tmp_path / "t.csv"
    expr.to_frame().T.to_csv(transposed)
    back = read_matrix(transposed, orientation="samples_x_genes")
    np.testing.assert_allclose(back.values, expr.values)


def test_read_breast_fixture(breast_dir: Path) -> None:
    expr = read_matrix(
        breast_dir / "Breastdata1" / "mRNA_abundance.txt.gz",
        sample_table=breast_dir / "Breastdata1" / "patient_annotation.txt",
    )
    assert (expr.n_genes, expr.n_samples) == (1231, 77)
    assert expr.genes[0] == "1000_at"
    assert expr.sample_table_provided
    assert expr.sample_table is not None and "t.rfs" in expr.sample_table.columns
    table = read_sample_table(breast_dir / "Breastdata1" / "patient_annotation.txt")
    assert len(table) == 77


@pytest.mark.parametrize("filename", ["m.csv", "m.tsv.gz"])
def test_text_identifiers_survive_matrix_round_trip(tmp_path: Path, filename: str) -> None:
    expr = ExpressionSet(
        np.array([[1.5, np.nan], [2, 3], [4, 5], [6, 7]]),
        pd.Index(["007", "7", "1.50", "NA"]),
        pd.Index(["001", "1"]),
    )
    path = tmp_path / filename
    write_matrix(expr, path)
    metadata = tmp_path / "metadata.csv"
    metadata.write_text("sample,age\n1,54\n001,NA\n")
    back = read_matrix(path, sample_table=metadata)
    assert list(back.genes) == list(expr.genes)
    assert list(back.samples) == list(expr.samples)
    np.testing.assert_equal(back.values, expr.values)
    assert back.sample_table is not None
    np.testing.assert_equal(back.sample_table["age"].to_numpy(), [np.nan, 54])


def test_seed_and_signature_ids_stay_distinct(tmp_path: Path) -> None:
    path = tmp_path / "seeds.csv"
    path.write_text("symbol,gene\n001,007\n1,7\n,1.50\nNA,NA\n,\n")
    spec = read_seeds(path)
    assert spec.seeds == ["007", "7", "1.50", "NA"]
    assert spec.symbols == {"007": "001", "7": "1", "NA": "NA"}
    assert read_gene_list(path) == spec.seeds


@pytest.mark.parametrize("header", [",s1,s2", "s1,s2"])
def test_row_names_with_an_unnamed_or_omitted_header(tmp_path: Path, header: str) -> None:
    path = tmp_path / "matrix.csv"
    path.write_text(f"{header}\n007,1,2\n7,3,4\nNA,5,NA\n")
    expr = read_matrix(path)
    assert expr.genes.tolist() == ["007", "7", "NA"]
    assert expr.samples.tolist() == ["s1", "s2"]
    np.testing.assert_equal(expr.values, [[1, 2], [3, 4], [5, np.nan]])


def test_strip_version_and_seeds(tmp_path: Path) -> None:
    assert strip_ensembl_version(["ENSG1.15", "ENSG2", "x.1.2"]) == ["ENSG1", "ENSG2", "x.1"]
    path = tmp_path / "seeds.csv"
    path.write_text("gene,symbol\nENSG1,ADM\nENSG2,CA9\nENSG1,ADM\n")
    spec = read_seeds(path)
    assert spec.seeds == ["ENSG1", "ENSG2"] and spec.symbols["ENSG2"] == "CA9"
    bad = tmp_path / "bad.csv"
    bad.write_text("id\nx\n")
    with pytest.raises(DataError):
        read_seeds(bad)
