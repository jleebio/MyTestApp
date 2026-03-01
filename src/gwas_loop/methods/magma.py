"""MAGMA gene-level association scoring."""
import numpy as np
import pandas as pd
from .base import BaseMethod


class MAGMAMethod(BaseMethod):
    """Gene-level association from GWAS summary stats via MAGMA.

    Converts MAGMA gene-level p-values to prioritization scores.
    Accounts for gene size and LD structure.
    """
    name = "magma"

    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        """Score from MAGMA gene-level p-values.

        Expects: gene_id, magma_p, magma_z (optional)
        """
        df = locus_data.copy()
        if "magma_z" in df.columns:
            z = df.set_index("gene_id")["magma_z"].abs()
        else:
            from scipy.stats import norm
            p = df.set_index("gene_id")["magma_p"].clip(lower=1e-300)
            z = norm.isf(p).abs()

        scores = z / z.max() if z.max() > 0 else z
        return scores.rename("magma_score")
