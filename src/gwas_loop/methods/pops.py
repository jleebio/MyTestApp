"""PoPS — Polygenic Priority Score."""
import numpy as np
import pandas as pd
from .base import BaseMethod


class PoPSMethod(BaseMethod):
    """Polygenic Priority Score: gene prioritization using gene-level features.

    Combines MAGMA gene associations with biological features
    (gene expression, chromatin, protein networks) via penalized regression.
    """
    name = "pops"

    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        """Score using pre-computed PoPS results.

        Expects: gene_id, pops_score
        """
        return locus_data.set_index("gene_id")["pops_score"].rename("pops_score")
