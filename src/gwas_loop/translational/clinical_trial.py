"""HIGH-IMPACT EXTENSION: Genetics → Clinical Trial Success Prediction.

Predict probability of clinical trial success based on strength
of genetic evidence supporting the target.

Key insight: Genetically supported targets are 2-3x more likely
to succeed in clinical trials (Nelson et al. 2015, King et al. 2019).
This module quantifies that probability.
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import cross_val_score
from dataclasses import dataclass


@dataclass
class TrialPrediction:
    gene_id: str
    disease: str
    drug: str
    predicted_success_probability: float
    phase_specific: dict[str, float]  # P(success) per phase
    genetic_support_tier: str
    key_drivers: list[str]
    comparable_successes: list[str]  # Similar targets that succeeded


class ClinicalTrialPredictor:
    """Predict clinical trial success from genetic evidence strength.

    Features derived from:
    1. GWAS effect size and significance
    2. Fine-mapping confidence
    3. Colocalization support
    4. Direction concordance (genetic → drug)
    5. LoF human evidence (natural experiments)
    6. Target tractability
    7. Mechanism clarity
    8. Cross-trait genetic support

    Baseline success rates (historical):
    - Phase I → II:  ~52%
    - Phase II → III: ~29%
    - Phase III → Approval: ~58%
    - Overall (I → Approval): ~9%

    With genetic support:
    - Overall: ~2-3x higher (~18-25%)
    """

    FEATURE_NAMES = [
        "gwas_pvalue_log10", "finemapping_pip", "coloc_h4",
        "direction_concordant", "lof_intolerance", "tractability_score",
        "mechanism_confidence", "n_supporting_traits", "tissue_relevance",
        "gene_constraint_oe", "has_human_lof_phenotype",
        "druggable_protein_class", "existing_drug_phase",
    ]

    HISTORICAL_BASE_RATES = {
        "phase1_to_phase2": 0.52,
        "phase2_to_phase3": 0.29,
        "phase3_to_approval": 0.58,
        "overall": 0.09,
    }

    GENETIC_SUPPORT_MULTIPLIERS = {
        "high": 2.5,    # Strong genetic evidence
        "medium": 1.8,  # Moderate evidence
        "low": 1.2,     # Weak evidence
        "none": 1.0,    # No genetic support
    }

    def __init__(self):
        self.model = GradientBoostingClassifier(
            n_estimators=200, max_depth=5, learning_rate=0.05,
            subsample=0.8, random_state=42
        )
        self.trained = False

    def extract_features(self, target_data: pd.DataFrame) -> pd.DataFrame:
        """Extract clinical success prediction features."""
        features = pd.DataFrame(index=target_data.index)
        for feat in self.FEATURE_NAMES:
            if feat in target_data.columns:
                features[feat] = target_data[feat]
            elif feat == "direction_concordant":
                features[feat] = target_data.get(feat, 0).astype(float)
            else:
                features[feat] = 0.0
        return features

    def train(self, features: pd.DataFrame, outcomes: pd.Series):
        """Train on historical trial outcomes."""
        X = self.extract_features(features) if set(features.columns) != set(self.FEATURE_NAMES) else features
        self.model.fit(X, outcomes)
        self.trained = True

    def predict(self, target_data: pd.DataFrame) -> list[TrialPrediction]:
        """Predict trial success for each target."""
        results = []
        for idx, row in target_data.iterrows():
            genetic_tier = self._classify_genetic_support(row)
            multiplier = self.GENETIC_SUPPORT_MULTIPLIERS[genetic_tier]

            # Phase-specific predictions
            phase_probs = {
                phase: min(base * multiplier, 0.95)
                for phase, base in self.HISTORICAL_BASE_RATES.items()
            }

            # Overall prediction
            if self.trained:
                X = self.extract_features(row.to_frame().T)
                p_success = float(self.model.predict_proba(X)[0, 1])
            else:
                p_success = phase_probs["overall"]

            drivers = self._identify_key_drivers(row)

            results.append(TrialPrediction(
                gene_id=str(row.get("gene_id", idx)),
                disease=str(row.get("disease", "unknown")),
                drug=str(row.get("drug", "unknown")),
                predicted_success_probability=round(p_success, 3),
                phase_specific=phase_probs,
                genetic_support_tier=genetic_tier,
                key_drivers=drivers,
                comparable_successes=[],
            ))

        return results

    def _classify_genetic_support(self, row) -> str:
        """Classify genetic support tier.

        Uses pre-computed GTCS tier if available, otherwise computes from components.
        """
        # Use pre-computed tier from GTCS if available
        tier = str(row.get("genetic_evidence_tier", ""))
        if tier in ("high", "medium", "low"):
            return tier

        # Use GTCS score directly if available
        gtcs = float(row.get("gtcs_score", 0))
        if gtcs > 0:
            if gtcs >= 0.7:
                return "high"
            elif gtcs >= 0.4:
                return "medium"
            elif gtcs > 0:
                return "low"

        # Fall back to component-level computation
        pip = float(row.get("finemapping_pip", row.get("pip", 0)))
        h4 = float(row.get("coloc_h4", 0))
        concordant = bool(row.get("direction_concordant", False))
        lof = float(row.get("lof_intolerance", 0))

        score = pip * 0.3 + h4 * 0.25 + (0.25 if concordant else 0) + lof * 0.2
        if score >= 0.6:
            return "high"
        elif score >= 0.35:
            return "medium"
        elif score >= 0.15:
            return "low"
        return "none"

    def _identify_key_drivers(self, row) -> list[str]:
        drivers = []
        if float(row.get("finemapping_pip", 0)) > 0.5:
            drivers.append("strong_finemapping")
        if float(row.get("coloc_h4", 0)) > 0.5:
            drivers.append("colocalization_support")
        if bool(row.get("direction_concordant", False)):
            drivers.append("direction_concordant")
        if float(row.get("has_human_lof_phenotype", 0)) > 0:
            drivers.append("human_LoF_evidence")
        if float(row.get("lof_intolerance", 0)) > 0.5:
            drivers.append("constrained_gene")
        return drivers

    def summary(self, predictions: list[TrialPrediction]) -> pd.DataFrame:
        return pd.DataFrame([
            {"gene": p.gene_id, "disease": p.disease, "drug": p.drug,
             "P(success)": p.predicted_success_probability,
             "tier": p.genetic_support_tier, "drivers": ", ".join(p.key_drivers)}
            for p in predictions
        ]).sort_values("P(success)", ascending=False)
