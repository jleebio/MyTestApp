"""Chromatin interaction methods (Hi-C, ABC model)."""
import numpy as np
import pandas as pd
from .base import BaseMethod


class ABCMethod(BaseMethod):
    """Activity-by-Contact model: enhancer-gene regulatory links.

    Scores genes by ABC score connecting GWAS variants to gene promoters
    via active enhancers.
    """
    name = "abc"

    def __init__(self, score_threshold: float = 0.015):
        self.score_threshold = score_threshold

    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        """Score genes by ABC enhancer-gene connections.

        Expects: gene_id, abc_score, is_enhancer_linked (bool)
        """
        df = locus_data[locus_data["abc_score"] >= self.score_threshold].copy()
        gene_scores = df.groupby("gene_id")["abc_score"].max()
        all_genes = locus_data["gene_id"].unique()
        return gene_scores.reindex(all_genes, fill_value=0.0).rename("abc_score")


class HiCMethod(BaseMethod):
    """Hi-C chromatin contact based gene prioritization."""
    name = "hic"

    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        """Expects: gene_id, hic_contact_frequency"""
        freq = locus_data.groupby("gene_id")["hic_contact_frequency"].max()
        scores = freq / freq.max() if freq.max() > 0 else freq
        return scores.rename("hic_score")
