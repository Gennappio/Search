"""Configuration models for every stage.

A RunConfig is a plain, validated, serialisable description of a run. Presets are concrete
RunConfig instances (see presets.py). Keys are addressed with dotted paths in overrides,
for example ``null.pooling=per_gene``.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Mapping
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator

from search.errors import ConfigError


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)


class AssociationConfig(_Model):
    """How the association between every seed and every gene is measured."""

    method: Literal["spearman", "pearson", "ols", "glm"] = "spearman"
    covariates: list[str] = Field(default_factory=list)
    glm_family: Literal["gaussian", "poisson", "gamma", "inverse_gaussian"] = "gaussian"
    formula: str | None = None


class EdgeConfig(_Model):
    """When an association becomes an edge of the seed network."""

    rule: Literal["abs_threshold", "p_threshold"] = "abs_threshold"
    threshold: float = 0.5
    alpha: float = 0.05
    bonferroni: Literal["n_genes", "none"] = "n_genes"
    p_adjust: Literal["none", "bonferroni", "bh"] = "none"
    self_edge: bool = True


class MembershipConfig(_Model):
    """How strongly a gene belongs to the cluster of a seed, given the edges."""

    mode: Literal["binary", "cluster_size", "gene_normalised", "sigmoid"] = "binary"
    sigmoid_k: float = 20.0
    sigmoid_x0: float = 0.5


class ConnectivityConfig(_Model):
    """How memberships are aggregated into one connectivity score per gene."""

    exclude_self: bool = True
    normalise: Literal["fraction", "count"] = "fraction"
    rank: Literal["none", "fractional"] = "none"


class OverlapConfig(_Model):
    """Seed-to-seed overlap measures. The binary Jaccard O is always computed."""

    shared_neighbourhood: Literal["auto", "always", "never"] = "auto"
    auto_max_seeds: int = 200


class GrowthConfig(_Model):
    """Recursive growth of the seed set (SEARCH 2014)."""

    enabled: bool = False
    max_iter: int = 1
    D_lower: float = 0.0
    D_higher: float = 0.0
    O_lower: float = 0.0
    keep_seeds: bool = True
    sd_threshold: float = 0.0


class ResamplingConfig(_Model):
    """Bootstrap or jackknife over samples."""

    scheme: Literal["none", "bootstrap", "jackknife"] = "none"
    n_rep: int = 0
    chunk_size: int = 100
    seed_base: int = 1_000_000
    statistic: Literal["C_rank", "C_raw"] = "C_rank"


class NullConfig(_Model):
    """Null model with random seed sets."""

    scheme: Literal["none", "internal_random_seeds", "external_runs"] = "none"
    n_rep: int = 0
    chunk_size: int = 100
    seed_base: int = 2_000_000
    pooling: Literal["pooled", "per_gene"] = "pooled"
    statistic: Literal["C_raw", "C_rank"] = "C_raw"
    eligible_min_sd: float = 0.0


class CombineConfig(_Model):
    """Combination of several training cohorts and signature selection."""

    enabled: bool = False
    mc_p_cutoff: float = 0.01
    product_cutoff: float = 0.5
    require_positive_consistent: bool = True
    exclude_seeds: bool = True
    top_n: int = 51
    var_floor: float = 1e-8


class ScoreConfig(_Model):
    """Per-sample signature score."""

    enabled: bool = False
    method: str = "median_abs_rank01"
    min_genes: int = 3

    @field_validator("method")
    @classmethod
    def _known_method(cls, value: str) -> str:
        if value == "median_abs_rank01" or (
            value.startswith("pysigscore_") and len(value) > len("pysigscore_")
        ):
            return value
        raise ValueError(
            "score.method must be 'median_abs_rank01' (built in) or 'pysigscore_<Scorer>' "
            "such as 'pysigscore_median' (needs pysigscore installed by hand)"
        )


# Fields that do not change the result and are ignored when matching a preset.
RUNTIME_FIELDS: dict[str, Any] = {
    "preset": True,
    "n_jobs": True,
    "checkpoint_dtype": True,
    "strict": True,
    "resampling": {"n_rep", "chunk_size", "seed_base"},
    "null": {"n_rep", "chunk_size", "seed_base"},
}

CORE_SECTIONS = ("association", "edge", "membership", "connectivity", "overlap")


class RunConfig(_Model):
    """Everything a run needs, apart from the data."""

    preset: str | None = None
    association: AssociationConfig = Field(default_factory=AssociationConfig)
    edge: EdgeConfig = Field(default_factory=EdgeConfig)
    membership: MembershipConfig = Field(default_factory=MembershipConfig)
    connectivity: ConnectivityConfig = Field(default_factory=ConnectivityConfig)
    overlap: OverlapConfig = Field(default_factory=OverlapConfig)
    growth: GrowthConfig = Field(default_factory=GrowthConfig)
    resampling: ResamplingConfig = Field(default_factory=ResamplingConfig)
    null: NullConfig = Field(default_factory=NullConfig)
    combine: CombineConfig = Field(default_factory=CombineConfig)
    score: ScoreConfig = Field(default_factory=ScoreConfig)
    n_jobs: int = 1
    checkpoint_dtype: Literal["float64", "float32"] = "float64"
    strict: bool = False

    def with_overrides(self, overrides: Mapping[str, Any]) -> RunConfig:
        """Copy with dotted-path overrides applied, e.g. {"null.pooling": "per_gene"}."""
        data = self.model_dump()
        for key, value in overrides.items():
            parts = key.split(".")
            node: Any = data
            for part in parts[:-1]:
                if not isinstance(node, dict) or part not in node:
                    raise ConfigError(f"unknown configuration path: {key}")
                node = node[part]
            if not isinstance(node, dict) or parts[-1] not in node:
                raise ConfigError(f"unknown configuration path: {key}")
            node[parts[-1]] = value
        try:
            return RunConfig.model_validate(data)
        except ValueError as exc:  # pydantic ValidationError is a ValueError
            raise ConfigError(str(exc)) from exc

    def section_dump(self, sections: Iterable[str]) -> dict[str, Any]:
        dump = self.model_dump(mode="json")
        return {name: dump[name] for name in sections}

    def hash(self, sections: Iterable[str] | None = None) -> str:
        """Short content hash of the result-relevant part of the configuration."""
        if sections is None:
            payload = self.model_dump(mode="json", exclude=RUNTIME_FIELDS)
        else:
            payload = self.section_dump(sections)
        text = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(text.encode()).hexdigest()[:16]

    def core_hash(self) -> str:
        return self.hash(CORE_SECTIONS)

    def to_yaml(self) -> str:
        return yaml.safe_dump(self.model_dump(mode="json"), sort_keys=False)

    @classmethod
    def from_yaml(cls, text: str) -> RunConfig:
        data = yaml.safe_load(text) or {}
        try:
            return cls.model_validate(data)
        except ValueError as exc:
            raise ConfigError(str(exc)) from exc


def parse_override(text: str) -> tuple[str, Any]:
    """Parse ``key.path=value`` from the command line; the value is read as YAML."""
    if "=" not in text:
        raise ConfigError(f"override must look like key.path=value, got {text!r}")
    key, raw = text.split("=", 1)
    return key.strip(), yaml.safe_load(raw)
