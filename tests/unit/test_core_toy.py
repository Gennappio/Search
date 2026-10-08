"""Hand-computed toy: 3 seeds (s1, s2, s3) and 3 other genes (g1, g2, g3).

Clusters (delta columns):  s1 -> {s1, s2, g1, g2}   s2 -> {s2, s1, g1}   s3 -> {s3, g2, g3}
Every association is +0.8 except g2 with s3, which is -0.7.

Binary membership, self excluded, fraction:
    D = s1 0.5, s2 0.5, s3 0, g1 2/3, g2 2/3, g3 1/3
Cluster-size membership (sizes 4, 3, 3):
    C = s1 1/6, s2 1/8, s3 0, g1 7/36, g2 7/36, g3 1/9
Binary Jaccard O: (s1,s2) 3/4, (s1,s3) 1/6, (s2,s3) 0; Omax 0.75, 0.75, 1/6
Weighted S with cluster-size gamma, the two seeds' rows excluded:
    (s1,s2) 3/7, (s1,s3) 3/14, (s2,s3) 0
Sign concordance F: (s1,s2) +1, (s1,s3) -1, (s2,s3) undefined
"""

import numpy as np
import pytest

from search.config import ConnectivityConfig, MembershipConfig
from search.core.connectivity import compute_connectivity, edge_statistics, nodes_degree
from search.core.membership import compute_membership
from search.core.overlap import jaccard_binary, omax, shared_neighbourhood

SEED_IDX = np.array([0, 1, 2])


def toy() -> tuple[np.ndarray, np.ndarray]:
    delta = np.zeros((6, 3), dtype=bool)
    delta[[0, 1, 3, 4], 0] = True
    delta[[0, 1, 3], 1] = True
    delta[[2, 4, 5], 2] = True
    stat = delta * 0.8
    stat[4, 2] = -0.7
    return delta, stat


def test_binary_connectivity_is_search_d() -> None:
    delta, stat = toy()
    gamma, sizes = compute_membership(delta, stat, MembershipConfig(mode="binary"))
    np.testing.assert_array_equal(sizes, [4, 3, 3])
    d = compute_connectivity(gamma, SEED_IDX, ConnectivityConfig())
    np.testing.assert_allclose(d, [0.5, 0.5, 0.0, 2 / 3, 2 / 3, 1 / 3])
    counts = compute_connectivity(gamma, SEED_IDX, ConnectivityConfig(normalise="count"))
    np.testing.assert_allclose(counts, [1, 1, 0, 2, 2, 1])
    with_self = compute_connectivity(gamma, SEED_IDX, ConnectivityConfig(exclude_self=False))
    np.testing.assert_allclose(with_self, [2 / 3, 2 / 3, 1 / 3, 2 / 3, 2 / 3, 1 / 3])


def test_cluster_size_connectivity() -> None:
    delta, stat = toy()
    gamma, _ = compute_membership(delta, stat, MembershipConfig(mode="cluster_size"))
    c = compute_connectivity(gamma, SEED_IDX, ConnectivityConfig())
    np.testing.assert_allclose(c, [1 / 6, 1 / 8, 0.0, 7 / 36, 7 / 36, 1 / 9])


def test_gene_normalised_membership_is_degenerate_with_step_edges() -> None:
    delta, stat = toy()
    gamma, _ = compute_membership(delta, stat, MembershipConfig(mode="gene_normalised"))
    np.testing.assert_allclose(gamma.sum(axis=1), 1.0)
    c = compute_connectivity(gamma, SEED_IDX, ConnectivityConfig(exclude_self=False))
    np.testing.assert_allclose(c, 1 / 3)  # every connected gene gets the same score


def test_sigmoid_membership_is_bounded() -> None:
    delta, stat = toy()
    gamma, _ = compute_membership(
        delta, stat, MembershipConfig(mode="sigmoid", sigmoid_k=10, sigmoid_x0=0.5)
    )
    assert gamma[delta].min() > 0 and gamma.max() <= 1 and (gamma[~delta] == 0).all()


def test_overlaps() -> None:
    delta, stat = toy()
    o = jaccard_binary(delta)
    np.testing.assert_allclose(o[0, 1], 3 / 4)
    np.testing.assert_allclose(o[0, 2], 1 / 6)
    np.testing.assert_allclose(o[1, 2], 0.0)
    np.testing.assert_allclose(np.diag(o), 1.0)
    np.testing.assert_allclose(omax(o), [0.75, 0.75, 1 / 6])
    gamma, _ = compute_membership(delta, stat, MembershipConfig(mode="cluster_size"))
    s, f = shared_neighbourhood(gamma, stat, SEED_IDX)
    np.testing.assert_allclose(s[0, 1], 3 / 7)
    np.testing.assert_allclose(s[0, 2], 3 / 14)
    np.testing.assert_allclose(s[1, 2], 0.0)
    np.testing.assert_allclose(s, s.T)
    assert f[0, 1] == pytest.approx(1.0) and f[0, 2] == pytest.approx(-1.0) and np.isnan(f[1, 2])


def test_binary_s_equals_jaccard_without_seed_rows() -> None:
    """With binary membership S is the Jaccard index computed on the non-seed rows."""
    delta, stat = toy()
    gamma, _ = compute_membership(delta, stat, MembershipConfig(mode="binary"))
    s, _ = shared_neighbourhood(gamma, stat, SEED_IDX)
    # (s1, s2): rows s3, g1, g2, g3 -> clusters {g1, g2} and {g1}: 1/2
    np.testing.assert_allclose(s[0, 1], 0.5)


def test_edge_statistics_and_degree() -> None:
    delta, stat = toy()
    gamma, _ = compute_membership(delta, stat, MembershipConfig(mode="binary"))
    es = edge_statistics(delta, stat, gamma)
    np.testing.assert_array_equal(es["n_edges"], [2, 2, 1, 2, 2, 1])
    np.testing.assert_array_equal(es["n_negative_edges"], [0, 0, 0, 0, 1, 0])
    np.testing.assert_array_equal(es["positive_consistent"], [1, 1, 1, 1, 0, 1])
    np.testing.assert_allclose(es["weighted_seed_rho"], [0.8, 0.8, 0.8, 0.8, 0.05, 0.8])
    np.testing.assert_allclose(nodes_degree(delta, SEED_IDX), [2 / 3, 2 / 3, 1 / 3])
