"""TWAS-inspired gene prioritization.

Transcriptome-wide association: prioritize genes whose predicted expression
correlates with trait. We implement a simplified scoring framework that
can consume external TWAS results (e.g., FUSION, PrediXcan, UTMOST).
"""
import numpy as np
import pandas as pd
from .base import BaseMethod


class TWASMethod(BaseMethod):
    """Score genes using TWAS z-scores across tissues.

    Input locus_data must have columns:
        gene_id, twas_z, twas_p, tissue (optional)
    If multiple tissues, takes the best per gene (most significant).
    """
    name = "twas"

    def __init__(self, aggregate: str = "best"):
        self.aggregate = aggregate  # best, mean, stouffer

    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        df = locus_data.copy()
        df["abs_z"] = df["twas_z"].abs()

        if "tissue" in df.columns and self.aggregate == "best":
            # Best tissue per gene
            idx = df.groupby("gene_id")["abs_z"].idxmax()
            df = df.loc[idx]
        elif "tissue" in df.columns and self.aggregate == "mean":
            df = df.groupby("gene_id").agg({"abs_z": "mean", "twas_p": "min"}).reset_index()
        elif "tissue" in df.columns and self.aggregate == "stouffer":
            grouped = df.groupby("gene_id")
            z_combined = grouped["twas_z"].apply(lambda x: x.sum() / np.sqrt(len(x)))
            df = z_combined.reset_index()
            df.columns = ["gene_id", "abs_z"]
            df["abs_z"] = df["abs_z"].abs()

        # Convert to 0-1 score using z-score scaling
        max_z = df["abs_z"].max()
        if max_z > 0:
            scores = df["abs_z"] / max_z
        else:
            scores = pd.Series(0.0, index=df.index)

        return pd.Series(scores.values, index=df["gene_id"].values, name="twas_score")
