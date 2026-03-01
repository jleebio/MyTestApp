"""INNOVATION 1: Disagreement-Guided Discovery.

Key insight: When methods disagree, the *pattern* of disagreement
is itself biologically informative. Coding loci show different
disagreement signatures than regulatory loci.

This is NOT just an ensemble — it explicitly models disagreement
as a predictive feature, which no existing method does.
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from dataclasses import dataclass


@dataclass
class DisagreementFeatures:
    """Features extracted from method disagreement patterns."""
    rank_variance: float          # How much methods disagree on rank
    direction_agreement: float    # Fraction of methods agreeing on top gene
    max_pairwise_divergence: float  # Worst-case disagreement
    entropy: float                # Information-theoretic uncertainty
    coding_regulatory_split: float  # Do coding methods agree but regulatory don't?


class DisagreementModel:
    """Model method disagreement as a predictive signal.

    Innovation: Instead of treating disagreement as noise to be averaged away
    (like standard ensembles do), we model disagreement explicitly.

    Hypothesis: The WAY methods disagree tells us about locus biology.
    - High distance + low coloc disagreement → coding variant, distance is right
    - High TWAS + low VEP disagreement → regulatory variant, TWAS is right
    - Universal disagreement → complex/polygenic locus, needs special handling
    """

    def __init__(self):
        self.model = GradientBoostingClassifier(
            n_estimators=200, max_depth=5, learning_rate=0.05, random_state=42
        )
        self.trained = False
        self.feature_names = [
            "rank_variance", "rank_range", "top1_agreement", "top3_overlap",
            "coding_vs_regulatory_discord", "pairwise_tau_mean", "pairwise_tau_min",
            "entropy", "max_score_gap", "n_methods_confident",
        ]

    def extract_features(self, method_scores: dict[str, pd.Series]) -> pd.DataFrame:
        """Extract disagreement features for each gene across methods."""
        score_df = pd.DataFrame(method_scores)
        rank_df = score_df.rank(ascending=False)

        features = pd.DataFrame(index=score_df.index)
        features["rank_variance"] = rank_df.var(axis=1)
        features["rank_range"] = rank_df.max(axis=1) - rank_df.min(axis=1)

        # Top-1 agreement: fraction of methods placing this gene in top 3
        features["top1_agreement"] = (rank_df <= 1).sum(axis=1) / rank_df.shape[1]
        features["top3_overlap"] = (rank_df <= 3).sum(axis=1) / rank_df.shape[1]

        # Coding vs regulatory method discord
        coding_methods = [c for c in ["distance", "vep", "finemap"] if c in score_df.columns]
        reg_methods = [c for c in ["twas", "coloc", "abc", "eqtl"] if c in score_df.columns]
        if coding_methods and reg_methods:
            coding_ranks = rank_df[coding_methods].mean(axis=1)
            reg_ranks = rank_df[reg_methods].mean(axis=1)
            features["coding_vs_regulatory_discord"] = (coding_ranks - reg_ranks).abs()
        else:
            features["coding_vs_regulatory_discord"] = 0.0

        # Pairwise rank correlation (Kendall's tau approximation via rank diff)
        from itertools import combinations
        methods = list(score_df.columns)
        if len(methods) >= 2:
            taus = []
            for m1, m2 in combinations(methods, 2):
                tau = rank_df[m1].corr(rank_df[m2], method="kendall")
                taus.append(tau if not np.isnan(tau) else 0.0)
            features["pairwise_tau_mean"] = np.mean(taus)
            features["pairwise_tau_min"] = np.min(taus)
        else:
            features["pairwise_tau_mean"] = 1.0
            features["pairwise_tau_min"] = 1.0

        # Entropy of normalized scores
        norm_scores = score_df.div(score_df.sum(axis=0), axis=1).clip(lower=1e-10)
        features["entropy"] = -(norm_scores * np.log(norm_scores)).sum(axis=1)

        features["max_score_gap"] = score_df.max(axis=1) - score_df.median(axis=1)
        features["n_methods_confident"] = (score_df > 0.5).sum(axis=1)

        return features[self.feature_names]

    def train(self, method_scores: dict[str, pd.Series], labels: pd.Series):
        """Train disagreement model on labeled loci."""
        X = self.extract_features(method_scores)
        self.model.fit(X, labels)
        self.trained = True

    def predict(self, method_scores: dict[str, pd.Series]) -> pd.Series:
        """Predict causal gene probability using disagreement features."""
        X = self.extract_features(method_scores)
        probs = self.model.predict_proba(X)[:, 1]
        return pd.Series(probs, index=X.index, name="disagreement_score")

    def feature_importance(self) -> pd.Series:
        if not self.trained:
            raise ValueError("Model not trained")
        return pd.Series(
            self.model.feature_importances_, index=self.feature_names
        ).sort_values(ascending=False)
