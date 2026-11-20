"""The files of SEARCH 0.0.2: ``write.table`` as the package calls it, the ``gzip`` it runs, and
``read.table`` as ``summarise.iterations`` calls it.

Every call of the package is ``write.table(x, file, row.names = TRUE, col.names = NA, quote =
FALSE, sep = "\\t")``, and every read is ``as.matrix(read.table(file, row.names = 1, header =
TRUE, sep = "\\t"))``. Only those two forms are reproduced.
"""

from __future__ import annotations

import gzip
import math
import re
from pathlib import Path

import numpy as np

from search.r2014.rformat import format_real, parse_real
from search.r2014.rsem import Name, NamedVector, Names, RError, RMatrix, UnportedState

_SAFE_PATH = re.compile(r"[A-Za-z0-9_./+,=:@%-]*")


def _name(value: Name) -> str:
    return "NA" if value is None else value


def _automatic(n: int) -> list[str]:
    return [str(i + 1) for i in range(n)]


def _cell(value: float) -> str:
    """A number as ``write.table`` writes it. Its test for a missing value is ``ISNAN``
    (``isna`` in ``utils/src/io.c``), so a ``NaN`` is written ``NA`` like a missing value."""
    return "NA" if math.isnan(value) else format_real(value)


def table_of_vector(vec: NamedVector) -> str:
    """``write.table`` of a numeric vector: R makes a data frame with one column ``x`` of it.

    The names become the row names when they are all different; otherwise the rows are
    numbered. A single missing name is an error of ``data.frame``.
    """
    names = vec.names
    if names is None or len(names) == 0 or all(n == "" for n in names):
        rows = _automatic(len(vec))
    elif len(set(names)) != len(names):
        rows = _automatic(len(vec))
    elif any(n is None for n in names):
        raise RError("row names contain missing values", "data.frame(x)")
    else:
        rows = [_name(n) for n in names]
    lines = ["\tx"]
    lines += [f"{r}\t{_cell(float(v))}" for r, v in zip(rows, vec.values, strict=True)]
    return "\n".join(lines) + "\n"


def table_of_strings(values: Names) -> str:
    """``write.table`` of a character vector: one column ``x``, rows numbered."""
    lines = ["\tx"]
    lines += [f"{i + 1}\t{_name(v)}" for i, v in enumerate(values)]
    return "\n".join(lines) + "\n"


def table_of_matrix(m: RMatrix) -> str:
    """``write.table`` of a numeric matrix. Without column names R writes V1, V2, ...; without
    row names it numbers the rows."""
    rows = [_name(n) for n in m.rownames] if m.rownames is not None else _automatic(m.nrow)
    if m.colnames is not None:
        cols = [_name(n) for n in m.colnames]
    else:
        cols = [f"V{j + 1}" for j in range(m.ncol)]
    lines = ["\t".join(["", *cols])] if m.ncol > 0 else [""]
    for i in range(m.nrow):
        lines.append(
            "\t".join([rows[i], *(_cell(float(v)) for v in m.values[i])])
            if m.ncol > 0
            else rows[i] + "\t"
        )
    return "\n".join(lines) + "\n"


def write_and_gzip(text: str, path: str) -> None:
    """Write ``path`` and then do what ``system("gzip --force --best <path>")`` does.

    The package pastes the path into a shell command without quoting it. With a blank in the
    path the shell splits it, gzip finds none of the pieces, R ignores the failure, and the
    file stays uncompressed, where ``summarise.iterations`` does not look for it (R defect
    R10). The port reproduces this; it does not run a shell, and refuses a path with other
    characters the shell would interpret.
    """
    Path(path).write_text(text, encoding="utf-8")
    if _SAFE_PATH.fullmatch(path):
        with gzip.GzipFile(path + ".gz", "wb", compresslevel=9, mtime=0) as handle:
            handle.write(text.encode("utf-8"))
        Path(path).unlink()
    elif re.fullmatch(r"[A-Za-z0-9_./+,=:@% \t-]*", path):
        return  # gzip was given pieces of the path: nothing is compressed
    else:
        raise UnportedState(
            f"output path {path!r} holds characters the shell would interpret; the package "
            "passes it to system() unquoted"
        )


_LOGICAL = {"T": "TRUE", "TRUE": "TRUE", "F": "FALSE", "FALSE": "FALSE"}
_INTEGER = re.compile(r"[+-]?[0-9]+")


def _row_names_as_r_reads_them(tokens: list[str]) -> list[str]:
    """The first column of a file after ``type.convert`` and ``as.character``.

    ``read.table`` converts every column before it takes the row names from the first one
    (R defect R13). When every name looks like a number the names come back as numbers:
    ``007`` is read as ``7``, ``1.50`` as ``1.5``, and ``100000`` as ``1e+05`` as soon as one
    other name has a decimal point. Names that are all ``T``, ``F``, ``TRUE`` or ``FALSE``
    come back as ``TRUE`` and ``FALSE``. A gene called ``NA`` stops the function.
    """
    if any(t == "NA" for t in tokens):
        raise RError("missing values in 'row.names' are not allowed")
    if not tokens:
        return []
    if all(t in _LOGICAL for t in tokens):
        return [_LOGICAL[t] for t in tokens]
    if all(_INTEGER.fullmatch(t) and abs(int(t)) < 2**31 for t in tokens):
        return [str(int(t)) for t in tokens]
    try:
        values = [parse_real(t) for t in tokens]
    except ValueError:
        return tokens
    return [format_real(v) for v in values]


def read_connectivity(path: str) -> RMatrix | None:
    """``as.matrix(read.table(file, row.names = 1, header = TRUE, sep = "\\t"))`` for a file
    written by the package: a one-column matrix with the genes as row names. ``None`` when the
    file does not exist."""
    file = Path(path)
    if not file.exists():
        return None
    with gzip.open(file, "rt", encoding="utf-8") as handle:
        text = handle.read()
    if any(c in text for c in "'\"#"):
        raise UnportedState("a file with quotes or # in it: read.table treats them specially")
    lines = [line for line in text.split("\n") if line != ""]
    if not lines:
        raise RError("no lines available in input")
    names: list[str] = []
    values: list[float] = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != 2:
            raise UnportedState("a connectivity file whose lines do not have two fields")
        names.append(fields[0])
        values.append(parse_real(fields[1]))
    rownames = _row_names_as_r_reads_them(names)
    if len(set(rownames)) != len(rownames):
        raise RError("duplicate 'row.names' are not allowed")
    return RMatrix(
        np.array(values, dtype=np.float64).reshape(len(values), 1),
        list(rownames) if rownames else None,
        ["x"],
    )


def read_matrix(path: str | Path, sep: str = "\t") -> RMatrix:
    """A data matrix from a text file, read as the examples of the package read theirs,
    ``as.matrix(read.table(file, row.names = 1, header = TRUE, sep = sep))``: genes in rows
    with their names in the first column, samples in columns, the numbers as R reads them
    (:func:`search.r2014.rformat.r_strtod`). The file may be compressed with gzip. Quotes around
    the names are removed; R also rewrites sample names that are not valid R names, which the
    package never looks at, and this function leaves them as they are."""
    file = Path(path)
    opener = gzip.open if file.suffix == ".gz" else open
    with opener(file, "rt", encoding="utf-8") as handle:
        lines = [line.rstrip("\r\n") for line in handle if line.strip() != ""]
    header = [field.strip('"') for field in lines[0].split(sep)]
    rows: list[str] = []
    values: list[list[float]] = []
    for line in lines[1:]:
        fields = line.split(sep)
        rows.append(fields[0].strip('"'))
        values.append([parse_real(field.strip('"')) for field in fields[1:]])
    n_columns = len(values[0]) if values else 0
    columns = header[-n_columns:] if n_columns else []
    return RMatrix(
        np.array(values, dtype=np.float64).reshape(len(rows), n_columns), list(rows), list(columns)
    )
