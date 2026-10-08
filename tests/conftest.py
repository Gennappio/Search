"""Shared fixtures. Real patient data never lives in the repository."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).parent
DATA_DIR = TESTS_DIR / "data"


@pytest.fixture(scope="session")
def breast_dir() -> Path:
    """Two small breast cancer microarray matrices bundled with the original SEARCH package."""
    return DATA_DIR / "breast"


@pytest.fixture(scope="session")
def gbm_dir() -> Path:
    """GBM cohorts folder, only through SEARCH_GBM_DATA; tests using it are skipped otherwise."""
    value = os.environ.get("SEARCH_GBM_DATA")
    if not value:
        pytest.skip("SEARCH_GBM_DATA not set")
    path = Path(value)
    if not path.is_dir():
        pytest.skip(f"SEARCH_GBM_DATA does not point to a directory: {path}")
    return path
