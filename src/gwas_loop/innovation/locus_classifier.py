"""INNOVATION 3: Locus Mechanism Classification.

Automatically classify GWAS loci into mechanism types and apply
DIFFERENT prioritization strategies per class.

No existing tool does this — they all apply one strategy everywhere.
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from dataclasses import dataclass
from enum import Enum


class LocusMechanism(str, Enum):
    CODING = "coding"           # Driven by protein-coding variants
    REGULATORY = "regulatory"   # Driven by gene regulation (enhancers, eQTLs)
    NETWORK = "network"         # Signal propagates through gene networks
    POLYGENIC = "polygenic"     # Multiple small-effect genes, no clear driver


@dataclass
class ClassifiedLocus:
    locus_id: str
    mechanism: LocusMechanism
    confidence: float
    recommended_methods: list[str]

    # Method recommendations per mechanism
    METHOD_MAP = {
        LocusMechanism.CODING: ["vep", "finemap", "distance"],
        LocusMechanism.REGULATORY: ["coloc", "twas", "abc", "eqtl"],
        LocusMechanism.NETWORK: ["network", "pops", "magma"],
        LocusMechanism.POLYGENIC: ["bayesian_ensemble", "pops"],
    }


class LocusClassifier:
    """Classify loci by mechanism and route to optimal methods.

    Features used:
    - Fraction of PIP in coding variants
    - eQTL signal strength
    - Number of genes in credible set
    - Network connectivity of top genes
    - LD structure complexity
    """

    FEATURE_NAMES = [
        "coding_pip_fraction", "max_eqtl_effect", "credible_set_size",
        "n_genes_in_locus", "top_gene_network_degree", "ld_complexity",
        "has_missense", "regulatory_element_count", "gene_density",
    ]

    def __init__(self):
        self.model = RandomForestClassifier(
            n_estimators=200, max_depth=8, random_state=42, class_weight="balanced"
        )
        self.trained = False

    def extract_features(self, locus_data: pd.DataFrame) -> np.ndarray:
        """Extract mechanism-informative features from a locus."""
        feats = {}
        feats["coding_pip_fraction"] = locus_data.get("coding_pip", pd.Series(0)).sum()
        feats["max_eqtl_effect"] = locus_data.get("eqtl_beta", pd.Series(0)).abs().max()
        feats["credible_set_size"] = locus_data.get("credible_set_size", pd.Series(1)).iloc[0]
        feats["n_genes_in_locus"] = locus_data["gene_id"].nunique()
        feats["top_gene_network_degree"] = locus_data.get("network_degree", pd.Series(0)).max()
        feats["ld_complexity"] = locus_data.get("ld_r2_mean", pd.Series(0.5)).mean()
        feats["has_missense"] = float(locus_data.get("has_missense", pd.Series(False)).any())
        feats["regulatory_element_count"] = locus_data.get("n_enhancers", pd.Series(0)).sum()
        feats["gene_density"] = feats["n_genes_in_locus"] / max(
            locus_data.get("locus_size_kb", pd.Series(1000)).iloc[0], 1
        )
        return np.array([[feats[f] for f in self.FEATURE_NAMES]])

    def train(self, loci_features: np.ndarray, labels: np.ndarray):
        self.model.fit(loci_features, labels)
        self.trained = True

    def classify(self, locus_data: pd.DataFrame, locus_id: str = "") -> ClassifiedLocus:
        """Classify a locus and recommend methods."""
        X = self.extract_features(locus_data)
        if self.trained:
            mech = LocusMechanism(self.model.predict(X)[0])
            conf = float(self.model.predict_proba(X).max())
        else:
            # Heuristic fallback before training
            mech = self._heuristic_classify(locus_data)
            conf = 0.5

        return ClassifiedLocus(
            locus_id=locus_id, mechanism=mech, confidence=conf,
            recommended_methods=ClassifiedLocus.METHOD_MAP[mech],
        )

    def _heuristic_classify(self, locus_data: pd.DataFrame) -> LocusMechanism:
        has_coding = locus_data.get("has_missense", pd.Series(False)).any()
        strong_eqtl = (locus_data.get("eqtl_beta", pd.Series(0)).abs() > 0.5).any()
        many_genes = locus_data["gene_id"].nunique() > 10

        if has_coding:
            return LocusMechanism.CODING
        elif strong_eqtl:
            return LocusMechanism.REGULATORY
        elif many_genes:
            return LocusMechanism.POLYGENIC
        else:
            return LocusMechanism.NETWORK
