"""SEARCH 0.0.2 (Haider and Buffa, 2014) as it was written: a verbatim port, defects included.

The rest of Search computes the method as it was meant: its ``search2014`` preset has none of
the defects of the R package and differs from it wherever a defect comes into play. This
subpackage is the opposite: every function of the package transliterated line by line, with
R's own rules for names, missing values, rounding, random numbers and number formatting, so
that a result of the 2014 package can be obtained again without R, and checked.

==============================  =====================================================
R                               here
==============================  =====================================================
``SEARCH``                      :func:`search.r2014.core.SEARCH`
``SEARCH.randomise``            :func:`search.r2014.randomise.SEARCH_randomise`
``summarise.iterations``        :func:`search.r2014.summarise.summarise_iterations`
``nodes.correlation``           :func:`search.r2014.association.nodes_correlation`
``nodes.glm``, ``nodes.fastLm``  :func:`search.r2014.glm.nodes_glm`, ``nodes_fastLm``
``getO``, ``getOmax``           :func:`search.r2014.core.getO`, ``getOmax``
``get.matrix.format``           :func:`search.r2014.core.get_matrix_format`
``estimate.nodes.degree``       :func:`search.r2014.core.estimate_nodes_degree`
``write.data.files``            :func:`search.r2014.randomise.write_data_files`
``removeNode``, ``removeByName``  :mod:`search.r2014.rsem`
==============================  =====================================================

Arguments keep their R names with the dots written as underscores (``D.lower`` is
``D_lower``). Where R stops with an error the port raises :class:`~search.r2014.rsem.RError`
with R's message. The figures of the package are not ported.

How far the port reproduces R, case by case and with the evidence, is in ``docs/replication.md``.
"""

from search.r2014.association import AssociationData, nodes_correlation
from search.r2014.core import (
    SEARCH,
    SearchResult,
    estimate_nodes_degree,
    get_matrix_format,
    getO,
    getOmax,
)
from search.r2014.glm import Covariates, nodes_fastLm, nodes_glm
from search.r2014.randomise import SEARCH_randomise, write_data_files
from search.r2014.rio import read_matrix
from search.r2014.rng import RRandom
from search.r2014.rsem import (
    NA_REAL,
    NamedVector,
    NotReproducible,
    RError,
    RList,
    RMatrix,
    UnportedState,
)
from search.r2014.summarise import IterationSummary, summarise_iterations

__all__ = [
    "NA_REAL",
    "SEARCH",
    "AssociationData",
    "Covariates",
    "IterationSummary",
    "NamedVector",
    "NotReproducible",
    "RError",
    "RList",
    "RMatrix",
    "RRandom",
    "SEARCH_randomise",
    "SearchResult",
    "UnportedState",
    "estimate_nodes_degree",
    "getO",
    "getOmax",
    "get_matrix_format",
    "nodes_correlation",
    "nodes_fastLm",
    "nodes_glm",
    "read_matrix",
    "summarise_iterations",
    "write_data_files",
]
