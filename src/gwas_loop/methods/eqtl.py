"""eQTL-based gene prioritization."""
import numpy as np
import pandas as pd
from .base import BaseMethod


class eQTLMethod(BaseMethod):
    """Direct eQTL mapping: genes whose expression is regulated by GWAS variants."""
    name = "eqtl"

    def __init__(self, aggregate: str = "best"):
        self.aggregate = aggregate

    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        """Score by eQTL effect size and significance.

        Expects: gene_id, eqtl_beta, eqtl_p, tissue (optional)
        """
        df = locus_data.copy()
        df["abs_beta"] = df["eqtl_beta"].abs()

        if "tissue" in df.columns and self.aggregate == "best":
            idx = df.groupby("gene_id")["abs_beta"].idxmax()
            df = df.loc[idx]

        gene_scores = df.groupby("gene_id")["abs_beta"].max()
        scores = gene_scores / gene_scores.max() if gene_scores.max() > 0 else gene_scores
        return scores.rename("eqtl_score")
