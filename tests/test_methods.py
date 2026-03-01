"""Tests for prioritization methods."""
import numpy as np
import pandas as pd
from gwas_loop.methods.distance import DistanceMethod, TSSSDistanceMethod
from gwas_loop.methods.twas import TWASMethod
from gwas_loop.methods.coloc import ColocMethod
from gwas_loop.methods.network import NetworkMethod
from gwas_loop.methods.ensemble import RankEnsemble, BayesianEnsemble


def _make_locus(n=10):
    return pd.DataFrame({
        "gene_id": [f"GENE{i}" for i in range(n)],
        "distance_to_lead": np.random.randint(100, 500000, n),
        "tss_distance": np.random.randint(100, 500000, n),
        "twas_z": np.random.randn(n) * 3,
        "twas_p": np.random.uniform(0, 1, n),
        "pp_h4": np.random.uniform(0, 1, n),
    })


def test_distance_method():
    locus = _make_locus()
    m = DistanceMethod()
    scores = m.score(locus)
    assert len(scores) == 10
    assert scores.max() <= 1.0
    assert scores.min() >= 0.0
    # Closest gene should score highest
    closest = locus["distance_to_lead"].idxmin()
    closest_gene = locus.loc[closest, "gene_id"]
    assert scores[closest_gene] == scores.max()


def test_tss_distance():
    locus = _make_locus()
    m = TSSSDistanceMethod()
    scores = m.score(locus)
    assert len(scores) == 10


def test_twas_single_tissue():
    locus = _make_locus()
    m = TWASMethod()
    scores = m.score(locus)
    assert len(scores) == 10
    assert scores.max() <= 1.0


def test_twas_multi_tissue():
    rows = []
    for tissue in ["brain", "liver", "blood"]:
        for i in range(5):
            rows.append({"gene_id": f"G{i}", "twas_z": np.random.randn() * 3,
                        "twas_p": np.random.uniform(), "tissue": tissue})
    df = pd.DataFrame(rows)
    for agg in ["best", "mean", "stouffer"]:
        m = TWASMethod(aggregate=agg)
        scores = m.score(df)
        assert len(scores) == 5


def test_coloc():
    locus = _make_locus()
    m = ColocMethod()
    scores = m.score(locus)
    assert len(scores) == 10
    assert scores.max() <= 1.0


def test_network_no_edges():
    locus = _make_locus()
    m = NetworkMethod()
    scores = m.score(locus)
    assert (scores == 0.0).all()


def test_network_with_edges():
    locus = pd.DataFrame({
        "gene_id": ["A", "B", "C", "D"],
        "is_seed": [True, False, False, False],
    })
    edges = pd.DataFrame({"gene1": ["A", "A", "B"], "gene2": ["B", "C", "D"], "weight": [1, 1, 1]})
    locus.attrs["edges"] = edges
    m = NetworkMethod()
    scores = m.score(locus)
    assert scores["A"] > scores["D"]  # Seed should rank high


def test_rank_ensemble():
    locus = _make_locus()
    e = RankEnsemble()
    e.add_method(DistanceMethod(), weight=1.0)
    e.add_method(ColocMethod(), weight=2.0)
    scores = e.score(locus)
    assert len(scores) == 10


def test_bayesian_ensemble():
    genes = [f"G{i}" for i in range(5)]
    be = BayesianEnsemble(prior=0.1)
    be.add_scores("m1", pd.Series([0.8, 0.2, 0.6, 0.1, 0.5], index=genes))
    be.add_scores("m2", pd.Series([0.7, 0.3, 0.5, 0.2, 0.9], index=genes))
    scores = be.score()
    assert len(scores) == 5
    assert scores["G0"] > scores["G3"]  # Both methods agree G0 > G3
