import pytest

from search.config import RunConfig, parse_override
from search.errors import ConfigError
from search.presets import BUFFA2010, SEARCH2014, get_preset, matching_preset


def test_overrides_and_errors() -> None:
    cfg = RunConfig().with_overrides({"null.pooling": "per_gene", "edge.threshold": 0.3})
    assert cfg.null.pooling == "per_gene" and cfg.edge.threshold == 0.3
    with pytest.raises(ConfigError):
        RunConfig().with_overrides({"null.nope": 1})
    with pytest.raises(ConfigError):
        RunConfig().with_overrides({"membership.mode": "weird"})


def test_parse_override_reads_yaml_values() -> None:
    assert parse_override("edge.alpha=0.01") == ("edge.alpha", 0.01)
    assert parse_override("growth.enabled=true") == ("growth.enabled", True)
    assert parse_override("association.covariates=[batch, age]") == (
        "association.covariates",
        ["batch", "age"],
    )
    with pytest.raises(ConfigError):
        parse_override("no-equals")


def test_hash_ignores_runtime_fields_and_tracks_results() -> None:
    base = BUFFA2010
    same = base.with_overrides({"n_jobs": 8, "resampling.n_rep": 3, "null.seed_base": 7})
    assert same.hash() == base.hash()
    assert same.core_hash() == base.core_hash()
    different = base.with_overrides({"membership.mode": "binary"})
    assert different.hash() != base.hash()
    assert different.core_hash() != base.core_hash()
    unrelated = base.with_overrides({"combine.top_n": 10})
    assert unrelated.core_hash() == base.core_hash() and unrelated.hash() != base.hash()


def test_yaml_round_trip() -> None:
    text = SEARCH2014.to_yaml()
    assert RunConfig.from_yaml(text) == SEARCH2014


def test_presets_match_themselves_and_not_each_other() -> None:
    assert matching_preset(BUFFA2010) == "buffa2010"
    assert matching_preset(SEARCH2014) == "search2014"
    assert matching_preset(BUFFA2010.with_overrides({"n_jobs": 4, "null.n_rep": 10})) == "buffa2010"
    assert matching_preset(BUFFA2010.with_overrides({"null.pooling": "per_gene"})) is None
    assert get_preset("buffa2010") == BUFFA2010
    with pytest.raises(ConfigError):
        get_preset("nope")
