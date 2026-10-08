"""Named, validated configurations. Only presets are validated; free configurations run with a
warning and are marked as such in the manifest."""

from __future__ import annotations

from typing import Any

from search.config import (
    RUNTIME_FIELDS,
    AssociationConfig,
    CombineConfig,
    ConnectivityConfig,
    EdgeConfig,
    GrowthConfig,
    MembershipConfig,
    NullConfig,
    OverlapConfig,
    ResamplingConfig,
    RunConfig,
    ScoreConfig,
)
from search.errors import ConfigError

BUFFA2010 = RunConfig(
    preset="buffa2010",
    association=AssociationConfig(method="spearman"),
    edge=EdgeConfig(
        rule="p_threshold", alpha=0.05, bonferroni="n_genes", p_adjust="none", self_edge=True
    ),
    membership=MembershipConfig(mode="cluster_size"),
    connectivity=ConnectivityConfig(exclude_self=True, normalise="fraction", rank="fractional"),
    overlap=OverlapConfig(shared_neighbourhood="auto"),
    growth=GrowthConfig(enabled=False),
    resampling=ResamplingConfig(
        scheme="bootstrap", n_rep=1000, chunk_size=100, seed_base=1_000_000, statistic="C_rank"
    ),
    null=NullConfig(
        scheme="internal_random_seeds",
        n_rep=5000,
        chunk_size=100,
        seed_base=2_000_000,
        pooling="pooled",
        statistic="C_raw",
        eligible_min_sd=0.0,
    ),
    combine=CombineConfig(
        enabled=True,
        mc_p_cutoff=0.01,
        product_cutoff=0.5,
        require_positive_consistent=True,
        exclude_seeds=True,
        top_n=51,
        var_floor=1e-8,
    ),
    score=ScoreConfig(enabled=True, method="median_abs_rank01", min_genes=3),
)
"""Buffa et al. 2010 as implemented by the GBM script: single pass, Bonferroni edges,
cluster-size membership, fractional ranks, bootstrap of ranks, pooled Monte Carlo null,
cross-cohort rank product and hypoxia score."""

SEARCH2014 = RunConfig(
    preset="search2014",
    association=AssociationConfig(method="spearman"),
    edge=EdgeConfig(rule="abs_threshold", threshold=0.5, self_edge=True),
    membership=MembershipConfig(mode="binary"),
    connectivity=ConnectivityConfig(exclude_self=True, normalise="fraction", rank="none"),
    overlap=OverlapConfig(shared_neighbourhood="auto"),
    growth=GrowthConfig(
        enabled=True,
        max_iter=1,
        D_lower=0.1,
        D_higher=0.24,
        O_lower=0.01,
        keep_seeds=True,
        sd_threshold=0.1,
    ),
    resampling=ResamplingConfig(
        scheme="bootstrap", n_rep=10, chunk_size=100, seed_base=0, statistic="C_raw"
    ),
    null=NullConfig(scheme="external_runs", pooling="per_gene", statistic="C_raw"),
    combine=CombineConfig(enabled=False),
    score=ScoreConfig(enabled=False),
)
"""The SEARCH R package with the defaults of SEARCH.randomise: hard correlation threshold,
binary membership, recursive growth with D and O pruning, bootstrap of D, external null runs."""

PRESETS: dict[str, RunConfig] = {"buffa2010": BUFFA2010, "search2014": SEARCH2014}


def get_preset(name: str) -> RunConfig:
    try:
        return PRESETS[name].model_copy(deep=True)
    except KeyError:
        raise ConfigError(f"unknown preset {name!r}; available: {sorted(PRESETS)}") from None


def _comparable(cfg: RunConfig) -> dict[str, Any]:
    return cfg.model_dump(mode="json", exclude=RUNTIME_FIELDS)


def matching_preset(cfg: RunConfig) -> str | None:
    """Name of the preset this configuration reproduces, ignoring runtime fields, or None."""
    target = _comparable(cfg)
    for name, preset in PRESETS.items():
        if _comparable(preset) == target:
            return name
    return None


def is_validated(cfg: RunConfig) -> bool:
    return matching_preset(cfg) is not None
