"""Detect failure patterns across method predictions."""
import numpy as np
import pandas as pd
from dataclasses import dataclass


@dataclass
class LocusDifficulty:
    locus_id: str
    disagreement_index: float
    uncertainty_score: float
    classification: str  # easy, moderate, hard


class FailureDetector:
    """Identify loci where methods disagree or predictions are unreliable."""

    @staticmethod
    def disagreement_index(predictions: dict[str, pd.Series]) -> pd.Series:
        """Compute per-gene disagreement across methods."""
        df = pd.DataFrame(predictions)
        ranks = df.rank(ascending=False)
        return ranks.std(axis=1) / ranks.mean(axis=1)

    @staticmethod
    def classify_loci(disagreement: pd.Series, thresholds: tuple[float, float] = (0.3, 0.7)) -> pd.Series:
        low, high = thresholds
        return disagreement.map(
            lambda d: "easy" if d < low else ("hard" if d > high else "moderate")
        )

    def detect(self, predictions: dict[str, pd.Series]) -> list[LocusDifficulty]:
        dis = self.disagreement_index(predictions)
        cls = self.classify_loci(dis)
        return [
            LocusDifficulty(locus_id=gene, disagreement_index=dis[gene],
                          uncertainty_score=dis[gene], classification=cls[gene])
            for gene in dis.index
        ]
