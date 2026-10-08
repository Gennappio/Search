"""Expression data container and readers.

The matrix is genes x samples. The sample table is validated against the sample names, so
covariates can never fall out of step with the columns, also after resampling.
"""

from __future__ import annotations

import gzip
import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from search.errors import DataError

Orientation = Literal["genes_x_samples", "samples_x_genes"]
TableFormat = Literal["csv", "tsv", "parquet", "h5ad"]

_VERSION_RE = re.compile(r"\.\d+$")


def strip_ensembl_version(ids: Iterable[Any]) -> list[str]:
    """ENSG00000000003.15 -> ENSG00000000003."""
    return [_VERSION_RE.sub("", str(g)) for g in ids]


@dataclass
class ExpressionSet:
    """A genes x samples matrix with aligned gene and sample tables."""

    values: NDArray[np.float64]
    genes: pd.Index
    samples: pd.Index
    gene_table: pd.DataFrame | None = None
    sample_table: pd.DataFrame | None = None
    sample_table_provided: bool | None = None

    def __post_init__(self) -> None:
        values = np.asarray(self.values, dtype=np.float64)
        if values.ndim != 2:
            raise DataError("expression values must be a 2-D array, genes x samples")
        self.values = np.ascontiguousarray(values)
        self.genes = pd.Index([str(g) for g in self.genes], name="gene")
        self.samples = pd.Index([str(s) for s in self.samples], name="sample")
        if len(self.genes) != values.shape[0]:
            raise DataError(f"{len(self.genes)} gene names for {values.shape[0]} rows")
        if len(self.samples) != values.shape[1]:
            raise DataError(f"{len(self.samples)} sample names for {values.shape[1]} columns")
        duplicated = self.genes[self.genes.duplicated()].unique()
        if len(duplicated):
            raise DataError(
                f"{len(duplicated)} duplicated gene identifiers, for example "
                f"{list(duplicated[:5])}; strip Ensembl versions or collapse duplicates upstream"
            )

        if self.gene_table is None:
            self.gene_table = pd.DataFrame(index=self.genes)
        else:
            table = self.gene_table.copy()
            table.index = pd.Index([str(g) for g in table.index])
            if table.index.has_duplicates:
                raise DataError("gene table index has duplicates")
            self.gene_table = table.reindex(self.genes)
        self.gene_table.index.name = "gene"

        if self.sample_table is None:
            self.sample_table = pd.DataFrame(index=self.samples)
            if self.sample_table_provided is None:
                self.sample_table_provided = False
        else:
            table = self.sample_table.copy()
            table.index = pd.Index([str(s) for s in table.index])
            if self.samples.has_duplicates:
                # a resampled matrix: the table must already be aligned row by row
                if not table.index.equals(self.samples):
                    raise DataError("a resampled matrix needs a sample table aligned row by row")
            else:
                if table.index.has_duplicates:
                    raise DataError("sample table index has duplicates")
                missing = self.samples.difference(table.index)
                if len(missing):
                    raise DataError(
                        f"sample table lacks {len(missing)} samples present in the matrix, "
                        f"for example {list(missing[:5])}"
                    )
                table = table.reindex(self.samples)
            self.sample_table = table
            if self.sample_table_provided is None:
                self.sample_table_provided = True
        self.sample_table.index.name = "sample"

    # ---- shape -------------------------------------------------------------------------
    @property
    def n_genes(self) -> int:
        return int(self.values.shape[0])

    @property
    def n_samples(self) -> int:
        return int(self.values.shape[1])

    # ---- selection ----------------------------------------------------------------------
    def gene_positions(self, genes: Iterable[str]) -> NDArray[np.intp]:
        wanted = [str(g) for g in genes]
        idx = self.genes.get_indexer(pd.Index(wanted))
        if (idx < 0).any():
            missing = [g for g, i in zip(wanted, idx, strict=True) if i < 0]
            raise DataError(f"{len(missing)} genes not in the matrix, for example {missing[:5]}")
        return np.asarray(idx, dtype=np.intp)

    def subset_genes(self, genes: Iterable[str]) -> ExpressionSet:
        idx = self.gene_positions(genes)
        assert self.gene_table is not None
        return ExpressionSet(
            self.values[idx],
            self.genes[idx],
            self.samples,
            self.gene_table.iloc[idx],
            self.sample_table,
            self.sample_table_provided,
        )

    def take_samples(self, idx: Iterable[int]) -> ExpressionSet:
        """Columns by position, repeats allowed (bootstrap); the sample table follows."""
        positions = np.asarray(list(idx), dtype=np.intp)
        if positions.size and (positions.min() < 0 or positions.max() >= self.n_samples):
            raise DataError("sample positions out of range")
        assert self.sample_table is not None
        return ExpressionSet(
            self.values[:, positions],
            self.genes,
            self.samples[positions],
            self.gene_table,
            self.sample_table.iloc[positions],
            self.sample_table_provided,
        )

    @staticmethod
    def common_genes(*sets: ExpressionSet) -> pd.Index:
        """Genes present in every set, in the order of the first."""
        if not sets:
            return pd.Index([], name="gene")
        common = sets[0].genes
        for other in sets[1:]:
            common = common[common.isin(other.genes)]
        return common

    # ---- conversion --------------------------------------------------------------------
    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.values, index=self.genes, columns=self.samples)

    @classmethod
    def from_frame(
        cls,
        frame: pd.DataFrame,
        sample_table: pd.DataFrame | None = None,
        gene_table: pd.DataFrame | None = None,
    ) -> ExpressionSet:
        return cls(
            frame.to_numpy(dtype=np.float64), frame.index, frame.columns, gene_table, sample_table
        )

    def to_anndata(self) -> Any:
        """AnnData is samples x genes: the matrix is transposed at this boundary."""
        import anndata

        assert self.sample_table is not None and self.gene_table is not None
        return anndata.AnnData(
            X=np.ascontiguousarray(self.values.T),
            obs=self.sample_table.copy(),
            var=self.gene_table.copy(),
        )

    @classmethod
    def from_anndata(cls, adata: Any) -> ExpressionSet:
        x = adata.X
        if hasattr(x, "toarray"):
            x = x.toarray()
        values = np.asarray(x, dtype=np.float64).T
        return cls(values, adata.var_names, adata.obs_names, adata.var.copy(), adata.obs.copy())


# ---- seeds ---------------------------------------------------------------------------------
@dataclass
class SeedSpec:
    """Seed genes with optional symbols and an optional universe restriction."""

    seeds: list[str]
    symbols: dict[str, str]
    universe: list[str] | None = None


# ---- files --------------------------------------------------------------------------------
def detect_format(path: Path | str) -> TableFormat:
    p = Path(path)
    suffixes = [s.lower() for s in p.suffixes]
    if suffixes and suffixes[-1] == ".gz":
        suffixes = suffixes[:-1]
    ext = suffixes[-1] if suffixes else ""
    if ext == ".parquet":
        return "parquet"
    if ext == ".h5ad":
        return "h5ad"
    if ext in (".tsv", ".tab"):
        return "tsv"
    if ext == ".csv":
        return "csv"
    # .txt or unknown: sniff the first line
    opener = gzip.open if p.suffix.lower() == ".gz" else open
    with opener(p, "rt") as handle:
        first = handle.readline()
    return "tsv" if "\t" in first else "csv"


def read_table(path: Path | str, format: TableFormat | None = None) -> pd.DataFrame:
    """A table whose first column holds row identifiers (csv, tsv, parquet)."""
    p = Path(path)
    if not p.exists():
        raise DataError(f"file not found: {p}")
    fmt = format or detect_format(p)
    if fmt == "parquet":
        frame = pd.read_parquet(p)
        if isinstance(frame.index, pd.RangeIndex) and frame.shape[1] > 0:
            frame = frame.set_index(frame.columns[0])
    elif fmt in ("csv", "tsv"):
        frame = pd.read_csv(p, sep="\t" if fmt == "tsv" else ",", index_col=0)
    else:
        raise DataError(f"read_table does not handle format {fmt}")
    frame.index = pd.Index([str(i) for i in frame.index])
    return frame


def read_sample_table(path: Path | str) -> pd.DataFrame:
    """One row per sample, first column = sample identifier."""
    table = read_table(path)
    if table.index.has_duplicates:
        raise DataError(f"duplicated sample identifiers in {path}")
    table.index.name = "sample"
    return table


def read_gene_list(path: Path | str, column: str = "gene") -> list[str]:
    """A gene list: a table with a ``gene`` column, or a one-column file without header."""
    p = Path(path)
    with gzip.open(p, "rt") if p.suffix.lower() == ".gz" else open(p) as handle:
        lines = [line.strip() for line in handle]
    lines = [line for line in lines if line]
    if not lines:
        return []
    sep = "\t" if "\t" in lines[0] else ","
    header = [h.strip().strip('"') for h in lines[0].split(sep)]
    if column in header:
        frame = pd.read_csv(p, sep=sep)
        return [str(g).strip() for g in frame[column].dropna()]
    return [line.split(sep)[0].strip().strip('"') for line in lines]


def read_seeds(path: Path | str) -> SeedSpec:
    """Seeds from a table with a ``gene`` column and an optional ``symbol`` column."""
    p = Path(path)
    sep = "\t" if detect_format(p) == "tsv" else ","
    frame = pd.read_csv(p, sep=sep)
    if "gene" not in frame.columns:
        raise DataError(f"seed file {p} needs a 'gene' column")
    seeds = [str(g).strip() for g in frame["gene"].dropna()]
    symbols: dict[str, str] = {}
    if "symbol" in frame.columns:
        for gene, symbol in zip(frame["gene"], frame["symbol"], strict=True):
            if pd.notna(gene) and pd.notna(symbol):
                symbols[str(gene).strip()] = str(symbol).strip()
    return SeedSpec(seeds=list(dict.fromkeys(seeds)), symbols=symbols)


def read_matrix(
    path: Path | str,
    *,
    format: TableFormat | None = None,
    orientation: Orientation = "genes_x_samples",
    strip_version: bool = False,
    sample_table: Path | str | pd.DataFrame | None = None,
    gene_table: pd.DataFrame | None = None,
) -> ExpressionSet:
    """Read an expression matrix from csv/tsv/parquet/h5ad."""
    p = Path(path)
    fmt = format or detect_format(p)
    if fmt == "h5ad":
        import anndata

        expr = ExpressionSet.from_anndata(anndata.read_h5ad(p))
        frame = expr.to_frame()
        if orientation == "samples_x_genes":
            frame = frame.T
    else:
        frame = read_table(p, fmt)
        if orientation == "samples_x_genes":
            frame = frame.T
    genes: Iterable[Any] = frame.index
    if strip_version:
        genes = strip_ensembl_version(genes)
    table: pd.DataFrame | None
    if isinstance(sample_table, (str, Path)):
        table = read_sample_table(sample_table)
    else:
        table = sample_table
    try:
        values = frame.to_numpy(dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise DataError(f"non-numeric values in {p}: {exc}") from exc
    return ExpressionSet(values, pd.Index(list(genes)), frame.columns, gene_table, table)


def write_matrix(expr: ExpressionSet, path: Path | str, format: TableFormat | None = None) -> None:
    p = Path(path)
    fmt = format or detect_format(p)
    frame = expr.to_frame()
    if fmt == "parquet":
        frame.reset_index().to_parquet(p, index=False)
    elif fmt == "csv":
        frame.to_csv(p)
    elif fmt == "tsv":
        frame.to_csv(p, sep="\t")
    elif fmt == "h5ad":
        expr.to_anndata().write_h5ad(p)
    else:
        raise DataError(f"unknown format {fmt}")
