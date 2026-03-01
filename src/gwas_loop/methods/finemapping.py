"""Fine-mapping based gene prioritization (SuSiE/FINEMAP wrapper)."""
import numpy as np
import pandas as pd
from .base import BaseMethod


class FinemapMethod(BaseMethod):
    """Prioritize genes using fine-mapping posterior inclusion probabilities.

    Aggregates variant-level PIPs to gene-level scores by summing PIPs
    of variants mapped to each gene (coding, intronic, or regulatory).
    """
    name = "finemap"

    def __init__(self, pip_threshold: float = 0.01, method: str = "susie"):
        self.pip_threshold = pip_threshold
        self.method = method  # susie, finemap, polyfun

    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        """Gene-level PIP aggregation.

        Expects: gene_id, variant_pip, variant_gene_mapping (e.g. coding/regulatory).
        Multiple variants can map to one gene — sum their PIPs (capped at 1.0).
        """
        df = locus_data[locus_data["variant_pip"] >= self.pip_threshold].copy()
        gene_pip = df.groupby("gene_id")["variant_pip"].sum().clip(upper=1.0)
        # Fill missing genes with 0
        all_genes = locus_data["gene_id"].unique()
        return gene_pip.reindex(all_genes, fill_value=0.0).rename("finemap_score")


class CredibleSetMethod(BaseMethod):
    """Score genes by credible set membership and size."""
    name = "credible_set"

    def __init__(self, coverage: float = 0.95):
        self.coverage = coverage

    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        """Genes containing credible set variants score higher.
        Smaller credible sets = higher confidence.

        Expects: gene_id, in_credible_set (bool), credible_set_size
        """
        df = locus_data.copy()
        in_cs = df.groupby("gene_id")["in_credible_set"].any().astype(float)
        avg_cs_size = df.groupby("gene_id")["credible_set_size"].mean().clip(lower=1)
        # Score: in credible set, penalized by set size
        scores = in_cs / np.log2(avg_cs_size + 1)
        scores = scores / scores.max() if scores.max() > 0 else scores
        return scores.rename("credible_set_score")
