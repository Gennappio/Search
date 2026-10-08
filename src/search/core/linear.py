"""Linear-model association (OLS with covariates, GLM families). Implemented in phase 4."""

from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from search.config import AssociationConfig
from search.core.association import AssociationResult


def associate_linear(
    values: NDArray[np.float64],
    seed_idx: NDArray[np.intp],
    cfg: AssociationConfig,
    sample_table: pd.DataFrame | None,
) -> AssociationResult:
    raise NotImplementedError("ols/glm association is implemented in phase 4")
