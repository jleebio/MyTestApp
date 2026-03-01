"""MODULE 8: Genetic Therapeutic Confidence Score (GTCS).

Composite score integrating all evidence layers for a therapeutic hypothesis.
"""
import numpy as np
from dataclasses import dataclass


@dataclass
class ConfidenceScore:
    gene_id: str
    disease: str
    total_score: float
    components: dict[str, float]
    tier: str  # "high", "medium", "low"


class GeneticTherapeuticScorer:
    """Compute Genetic Therapeutic Confidence Score.

    Components:
    1. Fine-mapping certainty (PIP)
    2. Colocalization support (PP.H4)
    3. Direction consistency (genetic vs drug)
    4. Tissue relevance (adaptive tissue weight)
    5. Pathway coherence (biological plausibility)
    6. Human LoF evidence (constraint, LoF intolerance)
    """

    WEIGHTS = {
        "finemapping": 0.20,
        "colocalization": 0.15,
        "direction": 0.20,
        "tissue": 0.15,
        "pathway": 0.15,
        "lof_evidence": 0.15,
    }

    def score(
        self,
        gene_id: str,
        disease: str,
        finemapping_pip: float = 0.0,
        coloc_h4: float = 0.0,
        direction_concordant: bool = False,
        tissue_weight: float = 0.0,
        pathway_score: float = 0.0,
        lof_intolerance: float = 0.0,
    ) -> ConfidenceScore:

        components = {
            "finemapping": min(finemapping_pip, 1.0),
            "colocalization": min(coloc_h4, 1.0),
            "direction": 1.0 if direction_concordant else 0.2,
            "tissue": min(tissue_weight, 1.0),
            "pathway": min(pathway_score, 1.0),
            "lof_evidence": min(lof_intolerance, 1.0),
        }

        total = sum(components[k] * self.WEIGHTS[k] for k in self.WEIGHTS)

        if total >= 0.7:
            tier = "high"
        elif total >= 0.4:
            tier = "medium"
        else:
            tier = "low"

        return ConfidenceScore(
            gene_id=gene_id, disease=disease,
            total_score=round(total, 3), components=components, tier=tier,
        )
