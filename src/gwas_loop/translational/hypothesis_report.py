"""MODULE 9: Therapeutic Hypothesis Report Generator.

Produces structured, publication-ready hypothesis outputs.
"""
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from datetime import datetime, timezone


@dataclass
class TherapeuticHypothesis:
    trait: str
    locus_id: str
    causal_gene: str
    mechanism: str
    direction_of_effect: str
    drug_candidate: str
    drug_action: str
    repurposing_rationale: str
    confidence_score: float
    confidence_tier: str
    tissue: str
    supporting_evidence: list[str]
    timestamp: str = ""

    def __post_init__(self):
        if not self.timestamp:
            self.timestamp = datetime.now(timezone.utc).isoformat()


class HypothesisReporter:
    """Generate and store therapeutic hypothesis reports."""

    def __init__(self, output_dir: str | Path = "reports"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.hypotheses: list[TherapeuticHypothesis] = []

    def add(self, hypothesis: TherapeuticHypothesis):
        self.hypotheses.append(hypothesis)

    def to_dataframe(self):
        import pandas as pd
        return pd.DataFrame([asdict(h) for h in self.hypotheses])

    def save_json(self, filename: str = "hypotheses.json"):
        path = self.output_dir / filename
        path.write_text(json.dumps([asdict(h) for h in self.hypotheses], indent=2))
        return path

    def save_tsv(self, filename: str = "hypotheses.tsv"):
        df = self.to_dataframe()
        path = self.output_dir / filename
        df.to_csv(path, sep="\t", index=False)
        return path

    def summary(self) -> dict:
        if not self.hypotheses:
            return {"total": 0}
        tiers = {}
        for h in self.hypotheses:
            tiers[h.confidence_tier] = tiers.get(h.confidence_tier, 0) + 1
        traits = list(set(h.trait for h in self.hypotheses))
        drugs = list(set(h.drug_candidate for h in self.hypotheses))
        return {
            "total": len(self.hypotheses),
            "by_tier": tiers,
            "traits_covered": traits,
            "unique_drugs": len(drugs),
            "top_scored": sorted(self.hypotheses, key=lambda h: -h.confidence_score)[:5],
        }

    def format_hypothesis(self, h: TherapeuticHypothesis) -> str:
        return (
            f"━━━ Therapeutic Hypothesis ━━━\n"
            f"Trait:       {h.trait}\n"
            f"Locus:       {h.locus_id}\n"
            f"Gene:        {h.causal_gene}\n"
            f"Mechanism:   {h.mechanism} ({h.direction_of_effect})\n"
            f"Drug:        {h.drug_candidate} ({h.drug_action})\n"
            f"Rationale:   {h.repurposing_rationale}\n"
            f"Tissue:      {h.tissue}\n"
            f"Confidence:  {h.confidence_score:.2f} [{h.confidence_tier}]\n"
            f"Evidence:    {', '.join(h.supporting_evidence)}\n"
        )
