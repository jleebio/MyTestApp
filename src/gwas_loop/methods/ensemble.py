"""Ensemble method combining multiple prioritization scores."""
import numpy as np
import pandas as pd
from .base import BaseMethod


class RankEnsemble(BaseMethod):
    """Combine methods using rank aggregation.

    Converts each method's scores to ranks, then averages.
    Robust to different score scales across methods.
    """
    name = "rank_ensemble"

    def __init__(self, methods: list[BaseMethod] | None = None, weights: dict[str, float] | None = None):
        self.methods = methods or []
        self.weights = weights or {}

    def add_method(self, method: BaseMethod, weight: float = 1.0):
        self.methods.append(method)
        self.weights[method.name] = weight

    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        """Weighted rank aggregation across methods."""
        if not self.methods:
            raise ValueError("No methods in ensemble")

        rank_df = pd.DataFrame()
        for method in self.methods:
            scores = method.score(locus_data)
            ranks = scores.rank(ascending=False, method="min")
            w = self.weights.get(method.name, 1.0)
            rank_df[method.name] = ranks * w

        total_weight = sum(self.weights.get(m.name, 1.0) for m in self.methods)
        avg_rank = rank_df.sum(axis=1) / total_weight
        # Invert so higher = better
        final_score = avg_rank.max() - avg_rank + 1
        final_score = final_score / final_score.max()
        return final_score.rename("ensemble_score")


class BayesianEnsemble(BaseMethod):
    """Bayesian score combination using calibrated probabilities."""
    name = "bayesian_ensemble"

    def __init__(self, prior: float = 0.05):
        self.prior = prior
        self.method_scores: dict[str, pd.Series] = {}

    def add_scores(self, method_name: str, scores: pd.Series):
        self.method_scores[method_name] = scores

    def score(self, locus_data: pd.DataFrame = None) -> pd.Series:
        """Combine via naive Bayes assumption (independent evidence)."""
        if not self.method_scores:
            raise ValueError("No scores added")

        genes = list(self.method_scores.values())[0].index
        log_odds = np.log(self.prior / (1 - self.prior))

        for name, scores in self.method_scores.items():
            s = scores.clip(lower=0.01, upper=0.99)
            log_odds = log_odds + np.log(s / (1 - s))

        prob = 1 / (1 + np.exp(-log_odds))
        return pd.Series(prob.values, index=genes, name="bayesian_score")
