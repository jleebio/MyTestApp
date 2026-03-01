"""Distance-based gene prioritization."""
import numpy as np
import pandas as pd
from .base import BaseMethod


class DistanceMethod(BaseMethod):
    """Prioritize genes by proximity to lead variant.

    Simple but effective baseline — nearest gene is causal ~30-40% of the time.
    Uses inverse-distance weighting with configurable decay.
    """
    name = "distance"

    def __init__(self, decay: float = 1e-5, max_distance: int = 1_000_000):
        self.decay = decay
        self.max_distance = max_distance

    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        """Score genes by distance to lead SNP.

        Expects columns: gene_id, distance_to_lead (bp).
        """
        dist = locus_data["distance_to_lead"].astype(float).clip(lower=1)
        scores = np.exp(-self.decay * dist)
        scores[dist > self.max_distance] = 0.0
        return pd.Series(scores.values, index=locus_data["gene_id"], name="distance_score")


class TSSSDistanceMethod(BaseMethod):
    """Prioritize by distance from lead variant to gene TSS."""
    name = "tss_distance"

    def __init__(self, decay: float = 5e-6):
        self.decay = decay

    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        dist = locus_data["tss_distance"].astype(float).clip(lower=1)
        scores = np.exp(-self.decay * dist)
        return pd.Series(scores.values, index=locus_data["gene_id"], name="tss_score")
