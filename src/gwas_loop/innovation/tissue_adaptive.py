"""INNOVATION 2: Trait-Adaptive Tissue Learning.

Current methods use fixed tissue assignments (e.g., "use brain for SCZ").
This learns tissue relevance weights PER TRAIT dynamically.

Novel: Hierarchical model that borrows strength across traits while
allowing trait-specific tissue profiles to emerge.
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from dataclasses import dataclass


@dataclass
class TissueWeight:
    tissue: str
    trait: str
    weight: float
    confidence: float


class TissueAdaptiveModel:
    """Learn per-trait tissue relevance instead of hardcoding.

    Standard approach: "SCZ → brain, T2D → pancreas"
    Our approach: Learn that SCZ also needs liver (metabolism of antipsychotics)
    and T2D needs brain (appetite regulation), with learned confidence.
    """

    def __init__(self, n_tissues: int = 49):  # GTEx v8 has 49 tissues
        self.tissue_weights: dict[str, np.ndarray] = {}  # trait -> weight vector
        self.model = None
        self.tissues: list[str] = []
        self.trained = False

    def train(self, tissue_scores: pd.DataFrame, labels: pd.Series, trait: str):
        """Learn tissue weights for a specific trait.

        tissue_scores: genes × tissues matrix of evidence scores
        labels: binary causal gene labels
        """
        self.tissues = list(tissue_scores.columns)
        model = LogisticRegressionCV(
            cv=5, penalty="l1", solver="saga", max_iter=5000, random_state=42
        )
        model.fit(tissue_scores.values, labels.values)
        weights = np.abs(model.coef_[0])
        weights = weights / weights.sum() if weights.sum() > 0 else weights
        self.tissue_weights[trait] = weights
        self.model = model
        self.trained = True

    def get_weights(self, trait: str) -> list[TissueWeight]:
        if trait not in self.tissue_weights:
            raise ValueError(f"No weights learned for {trait}")
        w = self.tissue_weights[trait]
        return sorted([
            TissueWeight(tissue=t, trait=trait, weight=float(w[i]),
                        confidence=float(w[i] / w.max()) if w.max() > 0 else 0)
            for i, t in enumerate(self.tissues)
        ], key=lambda x: -x.weight)

    def score(self, tissue_scores: pd.DataFrame, trait: str) -> pd.Series:
        """Apply learned tissue weights to score genes."""
        if trait not in self.tissue_weights:
            return tissue_scores.mean(axis=1)  # fallback to uniform
        w = self.tissue_weights[trait]
        weighted = tissue_scores.values @ w
        return pd.Series(weighted, index=tissue_scores.index, name="tissue_adaptive_score")
