"""Quality checks on the inputs. They warn, they never block, unless the run is strict."""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass

import numpy as np

from search.data import ExpressionSet
from search.stats import row_sds

_CONTROL_RE = re.compile(r"(ctrl|control|normal|healthy|benign)", re.IGNORECASE)
_REPLICATE_RE = re.compile(r"^(.+?)(?:\.|_rep|_r|-rep|-r)(\d{1,2})$", re.IGNORECASE)


@dataclass(frozen=True)
class QCWarning:
    code: str
    message: str
    count: int = 0

    def __str__(self) -> str:
        return f"[{self.code}] {self.message}"


def _examples(items: list[str], k: int = 5) -> str:
    shown = ", ".join(items[:k])
    return shown + (", ..." if len(items) > k else "")


def qc_expression(expr: ExpressionSet, name: str = "") -> list[QCWarning]:
    """Inspect a cohort and return warnings about anything that could bias the analysis."""
    tag = f"{name}: " if name else ""
    warnings: list[QCWarning] = []
    samples = [str(s) for s in expr.samples]

    if not expr.sample_table_provided:
        warnings.append(
            QCWarning(
                "sample_table_missing",
                f"{tag}no sample table: batches, patients, tissue types and clinical data "
                "are unknown; results cannot be checked against them",
            )
        )

    n_nan = int(np.isnan(expr.values).sum())
    if n_nan:
        warnings.append(
            QCWarning(
                "missing_values",
                f"{tag}{n_nan} missing values in the matrix; genes with missing values get "
                "no edges",
                n_nan,
            )
        )

    n_constant = int((row_sds(expr.values) == 0).sum())
    if n_constant:
        warnings.append(
            QCWarning(
                "constant_genes",
                f"{tag}{n_constant} genes are constant across samples and can never be "
                "associated with a seed",
                n_constant,
            )
        )

    dup = [s for s, c in Counter(samples).items() if c > 1]
    if dup:
        warnings.append(
            QCWarning(
                "duplicate_samples", f"{tag}duplicated sample names: {_examples(dup)}", len(dup)
            )
        )

    controls = [s for s in samples if _CONTROL_RE.search(s)]
    if controls:
        warnings.append(
            QCWarning(
                "control_like_names",
                f"{tag}{len(controls)} sample names look like controls or non-tumour tissue "
                f"({_examples(controls)}); if they are not tumours they distort co-expression",
                len(controls),
            )
        )

    roots: dict[str, list[str]] = {}
    for s in samples:
        match = _REPLICATE_RE.match(s)
        if match:
            roots.setdefault(match.group(1), []).append(s)
    groups = {root: members for root, members in roots.items() if len(members) > 1}
    if groups:
        example = next(iter(groups.values()))
        warnings.append(
            QCWarning(
                "replicate_like_names",
                f"{tag}{len(groups)} groups of sample names look like replicates of the same "
                f"patient ({_examples(example, 3)}); repeated patients break the independence "
                "assumed by correlations, p-values and bootstrap",
                len(groups),
            )
        )

    table = expr.sample_table
    if table is not None and expr.sample_table_provided:
        if "patient" in table.columns:
            counts = table["patient"].value_counts(dropna=True)
            repeated = counts[counts > 1]
            if len(repeated):
                warnings.append(
                    QCWarning(
                        "repeated_patients",
                        f"{tag}{len(repeated)} patients have more than one sample "
                        f"({_examples([str(p) for p in repeated.index])})",
                        len(repeated),
                    )
                )
        if "tissue" in table.columns:
            kinds = table["tissue"].astype(str).str.lower()
            non_tumour = int((~kinds.str.contains("tum")).sum())
            if non_tumour:
                warnings.append(
                    QCWarning(
                        "non_tumour_tissue",
                        f"{tag}{non_tumour} samples not labelled as tumour in the tissue column",
                        non_tumour,
                    )
                )
        if "batch" in table.columns:
            sizes = table["batch"].value_counts(dropna=True)
            if len(sizes) > 1 and (sizes.min() < 3 or sizes.min() / sizes.max() < 0.2):
                warnings.append(
                    QCWarning(
                        "unbalanced_batches",
                        f"{tag}tiny or unbalanced batches ({dict(sizes)}); check batch correction",
                        len(sizes),
                    )
                )
        for column in table.columns:
            n_missing = int(table[column].isna().sum())
            if n_missing:
                warnings.append(
                    QCWarning(
                        "missing_covariate_values",
                        f"{tag}sample table column '{column}' has {n_missing} missing values",
                        n_missing,
                    )
                )
    return warnings


def format_warnings(warnings: list[QCWarning]) -> str:
    if not warnings:
        return "QC: no warnings"
    return "QC warnings:\n" + "\n".join(f"  - {w}" for w in warnings)
