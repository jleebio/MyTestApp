"""MODULE 4: Biological Mechanism Inference.

Classify HOW a causal gene drives disease: coding disruption,
expression regulation, splicing, chromatin, or network perturbation.
"""
import numpy as np
import pandas as pd
from dataclasses import dataclass
from enum import Enum


class MechanismType(str, Enum):
    CODING_DISRUPTION = "coding_disruption"
    EXPRESSION_REGULATION = "expression_regulation"
    SPLICING_ALTERATION = "splicing_alteration"
    CHROMATIN_REGULATION = "chromatin_regulation"
    NETWORK_PERTURBATION = "network_perturbation"
    UNKNOWN = "unknown"


@dataclass
class MechanismAnnotation:
    gene_id: str
    mechanism: MechanismType
    confidence: float
    tissue: str
    cell_type: str
    direction_of_effect: str  # "gain", "loss", "unknown"
    evidence: list[str]


class MechanismInference:
    """Infer biological mechanism for each causal gene.

    Uses variant annotations, eQTL direction, splicing QTLs,
    chromatin marks, and network topology to classify mechanism.
    """

    CONSEQUENCE_TO_MECHANISM = {
        "frameshift_variant": MechanismType.CODING_DISRUPTION,
        "stop_gained": MechanismType.CODING_DISRUPTION,
        "missense_variant": MechanismType.CODING_DISRUPTION,
        "splice_donor_variant": MechanismType.SPLICING_ALTERATION,
        "splice_acceptor_variant": MechanismType.SPLICING_ALTERATION,
        "splice_region_variant": MechanismType.SPLICING_ALTERATION,
        "regulatory_region_variant": MechanismType.CHROMATIN_REGULATION,
        "TF_binding_site_variant": MechanismType.CHROMATIN_REGULATION,
    }

    def infer(self, gene_data: pd.DataFrame) -> list[MechanismAnnotation]:
        """Infer mechanism for each gene.

        Expects columns: gene_id, consequence (optional), eqtl_beta (optional),
        sqtl_p (optional), abc_score (optional), network_degree (optional),
        tissue (optional), cell_type (optional)
        """
        results = []
        for gene_id, group in gene_data.groupby("gene_id"):
            mechanism, confidence, evidence = self._classify_gene(group)
            direction = self._infer_direction(group)
            tissue = group["tissue"].mode().iloc[0] if "tissue" in group.columns and len(group["tissue"].dropna()) > 0 else "unknown"
            cell_type = group["cell_type"].iloc[0] if "cell_type" in group.columns and len(group["cell_type"].dropna()) > 0 else "unknown"

            results.append(MechanismAnnotation(
                gene_id=str(gene_id), mechanism=mechanism, confidence=confidence,
                tissue=tissue, cell_type=cell_type,
                direction_of_effect=direction, evidence=evidence,
            ))
        return results

    def _classify_gene(self, group: pd.DataFrame) -> tuple[MechanismType, float, list[str]]:
        evidence = []
        scores = {m: 0.0 for m in MechanismType}

        # Coding evidence
        if "consequence" in group.columns:
            for csq in group["consequence"].dropna():
                if csq in self.CONSEQUENCE_TO_MECHANISM:
                    mech = self.CONSEQUENCE_TO_MECHANISM[csq]
                    scores[mech] += 0.4
                    evidence.append(f"variant_consequence:{csq}")

        # eQTL evidence → expression regulation
        if "eqtl_beta" in group.columns:
            max_eqtl = group["eqtl_beta"].abs().max()
            if max_eqtl > 0.1:
                scores[MechanismType.EXPRESSION_REGULATION] += min(max_eqtl, 1.0) * 0.3
                evidence.append(f"eqtl_effect:{max_eqtl:.3f}")

        # Splicing QTL
        if "sqtl_p" in group.columns:
            min_sqtl = group["sqtl_p"].min()
            if min_sqtl < 1e-5:
                scores[MechanismType.SPLICING_ALTERATION] += 0.3
                evidence.append(f"sqtl_p:{min_sqtl:.2e}")

        # Chromatin / ABC
        if "abc_score" in group.columns:
            max_abc = group["abc_score"].max()
            if max_abc > 0.015:
                scores[MechanismType.CHROMATIN_REGULATION] += min(max_abc * 10, 0.4)
                evidence.append(f"abc_score:{max_abc:.3f}")

        # Network
        if "network_degree" in group.columns:
            max_deg = group["network_degree"].max()
            if max_deg > 10:
                scores[MechanismType.NETWORK_PERTURBATION] += 0.2
                evidence.append(f"network_degree:{max_deg}")

        if all(v == 0 for v in scores.values()):
            return MechanismType.UNKNOWN, 0.0, evidence

        best = max(scores, key=scores.get)
        total = sum(scores.values())
        confidence = scores[best] / total if total > 0 else 0.0
        return best, confidence, evidence

    def _infer_direction(self, group: pd.DataFrame) -> str:
        if "eqtl_beta" in group.columns:
            mean_beta = group["eqtl_beta"].mean()
            if mean_beta > 0.05:
                return "gain"
            elif mean_beta < -0.05:
                return "loss"
        return "unknown"
