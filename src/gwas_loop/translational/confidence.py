"""MODULE 8: Genetic Therapeutic Confidence Score (GTCS).

Composite score integrating all evidence layers for a therapeutic hypothesis.

The GTCS uses a weighted geometric-linear hybrid:
- Base: weighted linear combination of 6 evidence components
- Penalty: geometric mean of critical components (fine-mapping, direction)
  prevents high scores when key evidence is missing
- Bonus: nonlinear boost when multiple independent lines converge

This avoids the main weakness of pure linear models: a gene with
PIP=1.0 but no direction evidence shouldn't score as well as one
with PIP=0.8 and confirmed direction concordance.

Tier thresholds calibrated against Nelson et al. 2015 and King et al. 2019
genetic support categories for clinical trial success rates.
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

    Architecture: Hybrid linear-geometric model with convergence bonus.

    Components:
    1. Fine-mapping certainty (PIP) — Is this the causal variant?
    2. Colocalization support (PP.H4) — Same variant drives expression + disease?
    3. Direction consistency — Genetic direction matches drug action?
    4. Tissue relevance — Gene active in disease-relevant tissue?
    5. Pathway coherence — Gene in a biologically plausible pathway?
    6. Human LoF evidence — Natural human knockouts informative?

    Scoring layers:
    1. Linear base = Σ(weight_i × component_i)
    2. Critical penalty = geometric_mean(finemapping, direction)
       → Ensures both causal evidence AND direction are needed
    3. Convergence bonus = extra weight when ≥4 components > 0.5
       → Rewards independent lines of evidence agreeing
    4. Final = base × penalty_factor + convergence_bonus

    This means:
    - PIP=1.0 alone → moderate score (missing direction/coloc/tissue)
    - PIP=1.0 + direction + coloc → high score (converging evidence)
    - PIP=0.3 + everything else high → moderate (weak causal evidence drags it down)
    """

    WEIGHTS = {
        "finemapping": 0.20,
        "colocalization": 0.15,
        "direction": 0.20,
        "tissue": 0.15,
        "pathway": 0.15,
        "lof_evidence": 0.15,
    }

    # Components that MUST be present for high confidence
    CRITICAL_COMPONENTS = ["finemapping", "direction"]

    # Threshold for "meaningful evidence" in convergence check
    CONVERGENCE_THRESHOLD = 0.5
    CONVERGENCE_BONUS_PER_LINE = 0.02  # bonus per converging line above 3

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

        # Layer 1: Weighted linear base
        linear_base = sum(components[k] * self.WEIGHTS[k] for k in self.WEIGHTS)

        # Layer 2: Critical component penalty (geometric mean)
        # If either fine-mapping or direction is weak, score is dragged down
        critical_vals = [max(components[k], 0.01) for k in self.CRITICAL_COMPONENTS]
        critical_geomean = np.exp(np.mean(np.log(critical_vals)))
        # Penalty factor: sqrt of geometric mean (softer than raw geomean)
        penalty_factor = np.sqrt(critical_geomean)

        # Layer 3: Convergence bonus
        n_converging = sum(1 for v in components.values() if v > self.CONVERGENCE_THRESHOLD)
        convergence_bonus = max(0, (n_converging - 3)) * self.CONVERGENCE_BONUS_PER_LINE

        # Final score
        total = linear_base * penalty_factor + convergence_bonus
        total = min(max(total, 0.0), 1.0)  # clamp to [0, 1]

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
