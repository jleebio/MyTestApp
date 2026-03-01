"""MODULE 5: Target Tractability Analysis.

Assess whether a causal gene encodes a druggable protein.
"""
import pandas as pd
from dataclasses import dataclass
from enum import Enum


class TractabilityBucket(str, Enum):
    CLINICAL = "clinical_candidate"       # Known drug target in clinical trials
    DISCOVERY = "discovery_tractable"     # Druggable protein class, no clinical drug yet
    PREDICTED = "predicted_tractable"     # ML-predicted druggability
    INTRACTABLE = "currently_intractable" # No known druggability


class ProteinClass(str, Enum):
    KINASE = "kinase"
    GPCR = "gpcr"
    ION_CHANNEL = "ion_channel"
    NUCLEAR_RECEPTOR = "nuclear_receptor"
    ENZYME = "enzyme"
    TRANSPORTER = "transporter"
    SECRETED = "secreted"
    RECEPTOR = "receptor"
    TRANSCRIPTION_FACTOR = "transcription_factor"
    OTHER = "other"


@dataclass
class TractabilityResult:
    gene_id: str
    bucket: TractabilityBucket
    score: float  # 0-1, higher = more druggable
    protein_class: ProteinClass
    has_crystal_structure: bool
    has_active_site: bool
    is_secreted: bool
    antibody_accessible: bool
    evidence: list[str]


class TractabilityAnalyzer:
    """Evaluate druggability of gene targets.

    Integrates:
    - Protein class (kinase, GPCR, etc.)
    - ChEMBL target annotations
    - UniProt features (signal peptide, TM domains)
    - PDB structure availability
    - Open Targets tractability assessments
    """

    # Protein classes with established druggability
    DRUGGABLE_CLASSES = {
        ProteinClass.KINASE: 0.9,
        ProteinClass.GPCR: 0.85,
        ProteinClass.ION_CHANNEL: 0.8,
        ProteinClass.NUCLEAR_RECEPTOR: 0.85,
        ProteinClass.ENZYME: 0.7,
        ProteinClass.TRANSPORTER: 0.6,
        ProteinClass.RECEPTOR: 0.7,
        ProteinClass.SECRETED: 0.5,  # Antibody accessible
    }

    def analyze(self, gene_data: pd.DataFrame) -> list[TractabilityResult]:
        """Assess tractability for each gene.

        Expects: gene_id, protein_class (optional), has_structure (optional),
        chembl_targets (optional), is_secreted (optional), n_drugs (optional)
        """
        results = []
        for _, row in gene_data.drop_duplicates("gene_id").iterrows():
            gene_id = row["gene_id"]
            pclass = ProteinClass(row.get("protein_class", "other"))
            has_struct = bool(row.get("has_structure", False))
            is_secreted = bool(row.get("is_secreted", False))
            n_drugs = int(row.get("n_drugs", 0))

            score, bucket, evidence = self._score_gene(
                pclass, has_struct, is_secreted, n_drugs
            )

            results.append(TractabilityResult(
                gene_id=gene_id, bucket=bucket, score=score,
                protein_class=pclass, has_crystal_structure=has_struct,
                has_active_site=has_struct,  # Simplified
                is_secreted=is_secreted,
                antibody_accessible=is_secreted or pclass == ProteinClass.RECEPTOR,
                evidence=evidence,
            ))
        return results

    def _score_gene(self, pclass, has_struct, is_secreted, n_drugs):
        score = 0.0
        evidence = []

        if n_drugs > 0:
            score += 0.5
            evidence.append(f"existing_drugs:{n_drugs}")

        base = self.DRUGGABLE_CLASSES.get(pclass, 0.1)
        score += base * 0.3
        evidence.append(f"protein_class:{pclass.value}")

        if has_struct:
            score += 0.15
            evidence.append("has_crystal_structure")

        if is_secreted:
            score += 0.05
            evidence.append("secreted_protein")

        score = min(score, 1.0)

        if n_drugs > 0:
            bucket = TractabilityBucket.CLINICAL
        elif base >= 0.6:
            bucket = TractabilityBucket.DISCOVERY
        elif base >= 0.3:
            bucket = TractabilityBucket.PREDICTED
        else:
            bucket = TractabilityBucket.INTRACTABLE

        return score, bucket, evidence
