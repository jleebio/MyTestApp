"""INNOVATION 5: Cross-Trait Transfer Learning.

Train on multiple traits, test on unseen diseases.
Novel: Learns trait-invariant features for causal gene discovery,
enabling prediction for diseases without known causal genes.
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import LeaveOneGroupOut
from dataclasses import dataclass


@dataclass
class TransferResult:
    source_traits: list[str]
    target_trait: str
    auroc: float
    precision_at_5: float
    n_genes_tested: int
    top_predictions: list[str]


class CrossTraitTransfer:
    """Learn generalizable causal gene features across traits.

    Key idea: Some features predict causal genes regardless of trait
    (e.g., gene constraint, regulatory density). Learn these shared
    features, then apply to new traits with zero labeled data.
    """

    def __init__(self):
        self.model = GradientBoostingClassifier(
            n_estimators=300, max_depth=6, learning_rate=0.05,
            subsample=0.8, random_state=42
        )
        self.trained = False
        self.results: list[TransferResult] = []

    def train_leave_one_trait_out(
        self, features: pd.DataFrame, labels: pd.Series, traits: pd.Series
    ) -> list[TransferResult]:
        """Leave-one-trait-out cross-validation.

        features: gene × feature matrix (all traits pooled)
        labels: binary causal labels
        traits: trait assignment per gene
        """
        logo = LeaveOneGroupOut()
        results = []

        for train_idx, test_idx in logo.split(features, labels, traits):
            X_train, X_test = features.iloc[train_idx], features.iloc[test_idx]
            y_train, y_test = labels.iloc[train_idx], labels.iloc[test_idx]
            target_trait = traits.iloc[test_idx].iloc[0]
            source_traits = list(traits.iloc[train_idx].unique())

            self.model.fit(X_train, y_train)
            y_pred = self.model.predict_proba(X_test)[:, 1]

            from sklearn.metrics import roc_auc_score
            auroc = roc_auc_score(y_test, y_pred) if len(np.unique(y_test)) > 1 else 0.0
            top5_genes = features.iloc[test_idx].index[np.argsort(y_pred)[-5:]].tolist()
            p5 = (y_test.iloc[np.argsort(y_pred)[-5:]] == 1).mean()

            result = TransferResult(
                source_traits=source_traits, target_trait=target_trait,
                auroc=auroc, precision_at_5=float(p5),
                n_genes_tested=len(y_test), top_predictions=top5_genes,
            )
            results.append(result)

        self.results = results
        self.trained = True
        return results

    def predict_new_trait(self, features: pd.DataFrame) -> pd.Series:
        """Predict causal genes for a completely new trait."""
        if not self.trained:
            raise ValueError("Train first with train_leave_one_trait_out")
        probs = self.model.predict_proba(features)[:, 1]
        return pd.Series(probs, index=features.index, name="transfer_score")
