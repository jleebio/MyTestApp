"""Variant effect prediction scoring (CADD, PolyPhen, SIFT)."""
import numpy as np
import pandas as pd
from .base import BaseMethod


class VEPMethod(BaseMethod):
    """Prioritize genes by variant functional impact.

    Uses CADD scores, PolyPhen predictions, and consequence severity
    to identify genes harboring damaging variants.
    """
    name = "vep"

    CONSEQUENCE_WEIGHTS = {
        "frameshift_variant": 1.0,
        "stop_gained": 1.0,
        "splice_donor_variant": 0.95,
        "splice_acceptor_variant": 0.95,
        "missense_variant": 0.7,
        "inframe_deletion": 0.6,
        "inframe_insertion": 0.6,
        "regulatory_region_variant": 0.3,
        "3_prime_UTR_variant": 0.2,
        "5_prime_UTR_variant": 0.2,
        "synonymous_variant": 0.05,
        "intron_variant": 0.02,
        "intergenic_variant": 0.01,
    }

    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        """Combine CADD + consequence into gene scores.

        Expects: gene_id, cadd_phred (optional), consequence (optional)
        """
        df = locus_data.copy()
        scores = pd.Series(0.0, index=df["gene_id"].unique())

        if "cadd_phred" in df.columns:
            cadd = df.groupby("gene_id")["cadd_phred"].max() / 40.0  # Normalize to 0-1
            scores = scores.add(cadd, fill_value=0)

        if "consequence" in df.columns:
            df["csq_weight"] = df["consequence"].map(self.CONSEQUENCE_WEIGHTS).fillna(0.01)
            csq = df.groupby("gene_id")["csq_weight"].max()
            scores = scores.add(csq, fill_value=0)

        scores = scores / scores.max() if scores.max() > 0 else scores
        return scores.rename("vep_score")
