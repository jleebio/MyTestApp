"""MODULE 7: Drug Repurposing Engine.

Core logic:
  Gene activity ↑ → disease risk → Drug that ↓ gene → Therapeutic benefit
  Gene activity ↓ → disease risk → Drug that ↑ gene → Therapeutic benefit
"""
import numpy as np
import pandas as pd
from dataclasses import dataclass


@dataclass
class RepurposingHypothesis:
    gene_id: str
    drug_name: str
    disease: str
    genetic_direction: str      # "gain" or "loss" of function → risk
    drug_action: str            # "inhibitor", "agonist", etc.
    therapeutic_rationale: str  # Why this should work
    concordance: bool           # Does drug action oppose genetic direction?
    score: float                # Composite repurposing score
    original_indication: str
    evidence: list[str]


class RepurposingEngine:
    """Predict drug repurposing opportunities from genetic evidence.

    Innovation: Combines genetic direction of effect with drug mechanism
    to identify concordant therapeutic opportunities.
    """

    # Drug actions that REDUCE gene activity
    INHIBITORY_ACTIONS = {"inhibitor", "antagonist", "blocker", "negative_modulator", "antisense"}
    # Drug actions that INCREASE gene activity
    ACTIVATING_ACTIONS = {"agonist", "activator", "positive_modulator", "potentiator"}

    def generate_hypotheses(
        self,
        gene_id: str,
        disease: str,
        genetic_direction: str,
        drugs: list[dict],
        causal_confidence: float = 0.5,
    ) -> list[RepurposingHypothesis]:
        """Generate repurposing hypotheses for a gene-disease pair."""
        hypotheses = []

        for drug in drugs:
            moa = drug.get("mechanism_of_action", "unknown")
            concordance = self._check_concordance(genetic_direction, moa)
            rationale = self._build_rationale(gene_id, disease, genetic_direction, drug["drug_name"], moa, concordance)
            score = self._compute_score(concordance, causal_confidence, drug.get("max_phase", 0))

            hypotheses.append(RepurposingHypothesis(
                gene_id=gene_id,
                drug_name=drug["drug_name"],
                disease=disease,
                genetic_direction=genetic_direction,
                drug_action=moa,
                therapeutic_rationale=rationale,
                concordance=concordance,
                score=score,
                original_indication=drug.get("indication", "unknown"),
                evidence=drug.get("evidence", []),
            ))

        return sorted(hypotheses, key=lambda h: -h.score)

    def _check_concordance(self, genetic_dir: str, drug_action: str) -> bool:
        """Check if drug action opposes genetic risk direction."""
        if genetic_dir == "gain" and drug_action in self.INHIBITORY_ACTIONS:
            return True
        if genetic_dir == "loss" and drug_action in self.ACTIVATING_ACTIONS:
            return True
        return False

    def _build_rationale(self, gene, disease, direction, drug, moa, concordant):
        if concordant:
            return (f"Genetic evidence: {direction}-of-function in {gene} increases {disease} risk. "
                    f"{drug} ({moa}) opposes this direction → predicted therapeutic benefit.")
        else:
            return (f"Genetic evidence: {direction}-of-function in {gene} increases {disease} risk. "
                    f"{drug} ({moa}) does NOT clearly oppose genetic direction → lower confidence.")

    def _compute_score(self, concordant, causal_conf, max_phase):
        score = 0.0
        score += 0.4 * causal_conf                          # Genetic confidence
        score += 0.3 * (1.0 if concordant else 0.2)         # Direction concordance
        score += 0.3 * (max_phase / 4.0)                    # Clinical advancement
        return round(min(score, 1.0), 3)
