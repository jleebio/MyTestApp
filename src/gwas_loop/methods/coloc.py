"""Colocalization-based gene prioritization.

Tests whether GWAS and eQTL signals share the same causal variant.
This implements a scoring wrapper for external coloc results
(e.g., coloc, eCAVIAR, HyPrColoc).
"""
import numpy as np
import pandas as pd
from .base import BaseMethod


class ColocMethod(BaseMethod):
    """Score genes using colocalization posterior probabilities.

    Input locus_data columns:
        gene_id, pp_h4 (posterior prob of shared causal variant),
        tissue (optional), pp_h3 (optional, independent signals)
    """
    name = "coloc"

    def __init__(self, h4_threshold: float = 0.5, aggregate: str = "max"):
        self.h4_threshold = h4_threshold
        self.aggregate = aggregate  # max, mean across tissues

    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        df = locus_data.copy()

        if "tissue" in df.columns:
            if self.aggregate == "max":
                df = df.loc[df.groupby("gene_id")["pp_h4"].idxmax()]
            elif self.aggregate == "mean":
                df = df.groupby("gene_id").agg({"pp_h4": "mean"}).reset_index()

        scores = df["pp_h4"].clip(lower=0, upper=1)
        return pd.Series(scores.values, index=df["gene_id"].values, name="coloc_score")
