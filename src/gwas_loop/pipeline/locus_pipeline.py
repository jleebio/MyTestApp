"""Real end-to-end pipeline: GWAS locus → SuSiE fine-mapping → therapeutic hypothesis."""
import pandas as pd
import numpy as np

from ..methods.susie_real import SuSiERunner
from ..translational.confidence import GeneticTherapeuticScorer
from ..translational.clinical_trial import ClinicalTrialPredictor
from ..data.validated_genes import VALIDATED_GENES
from ..data.drug_database import get_drugs_for_gene


class LocusPipeline:
    """Process a single GWAS locus through fine-mapping → drugs."""

    def __init__(self, ld_ref_dir: str = "data/ld_ref", **kwargs):
        self.susie = SuSiERunner(ld_ref_dir=ld_ref_dir)
        self.scorer = GeneticTherapeuticScorer()
        self.trial_predictor = ClinicalTrialPredictor()

    def process_locus(self, sumstats, chrom, lead_bp, lead_snp, trait_id,
                      gene_map=None, candidate_genes=None):
        result = {"locus": f"chr{chrom}:{lead_bp}", "lead_snp": lead_snp,
                  "trait": trait_id, "stages": {}}

        # Stage 1: Fine-mapping with real SuSiE
        fm = self.susie.finemapping_locus(sumstats, chrom, lead_bp, gene_map=gene_map)
        result["stages"]["finemapping"] = {
            "n_variants": fm.get("n_variants", 0),
            "n_credible_sets": fm.get("n_credible_sets", 0),
            "converged": fm.get("converged", False),
            "top_variants": fm.get("top_variants", [])[:5],
        }

        # Use fine-mapping gene scores or provided candidates
        if candidate_genes is None and fm.get("gene_scores"):
            candidate_genes = [{"gene_id": g, "pip": p}
                               for g, p in sorted(fm["gene_scores"].items(), key=lambda x: -x[1])[:20]]
        if not candidate_genes:
            result["error"] = "No candidate genes identified"
            return result

        top_gene = candidate_genes[0]
        gid = top_gene["gene_id"]

        # Stage 2: Drug lookup
        drugs = get_drugs_for_gene(gid)
        drug_list = drugs.to_dict(orient="records") if len(drugs) > 0 else []
        result["stages"]["drugs"] = drug_list

        # Stage 3: Confidence scoring
        conf = self.scorer.score(
            gene_id=gid, disease=trait_id,
            finemapping_pip=top_gene.get("pip", 0),
            coloc_h4=top_gene.get("coloc_score", 0),
            direction_concordant=top_gene.get("direction_concordant", False),
            tissue_weight=top_gene.get("tissue_weight", 0.5),
            pathway_score=top_gene.get("pathway_score", 0.5),
            lof_intolerance=top_gene.get("lof_score", 0),
        )

        # Stage 4: Trial prediction
        trial = self.trial_predictor.predict(conf)

        validated = VALIDATED_GENES.get(trait_id)
        is_val = gid in (validated.genes if validated else set())

        result["hypothesis"] = {
            "causal_gene": gid,
            "drugs": [d["drug_name"] for d in drug_list[:3]],
            "confidence_score": conf.total_score,
            "confidence_tier": conf.tier,
            "p_trial_success": trial.get("p_success", 0.09),
            "is_validated": is_val,
        }
        return result


class FullPipeline:
    """Run pipeline across all significant loci for a trait."""

    def __init__(self, **kwargs):
        self.locus_pipeline = LocusPipeline(**kwargs)

    def identify_loci(self, sumstats, p_threshold=5e-8, window_kb=500):
        sig = sumstats[sumstats["P"] < p_threshold].copy().sort_values("P")
        loci, used = [], set()
        for _, row in sig.iterrows():
            key = f"{row['CHR']}:{int(row['BP']) // (window_kb * 1000)}"
            if key in used:
                continue
            used.add(key)
            loci.append({"chrom": row["CHR"], "lead_bp": int(row["BP"]),
                         "lead_snp": row.get("SNP", f"chr{row['CHR']}:{int(row['BP'])}"),
                         "lead_p": float(row["P"])})
        return loci

    def run_trait(self, sumstats, trait_id, gene_map=None, max_loci=50):
        loci = self.identify_loci(sumstats)[:max_loci]
        print(f"\n  {trait_id}: {len(loci)} genome-wide significant loci")
        results = []
        for i, loc in enumerate(loci):
            print(f"  [{i+1}/{len(loci)}] {loc['lead_snp']} (P={loc['lead_p']:.2e})")
            try:
                r = self.locus_pipeline.process_locus(
                    sumstats, loc["chrom"], loc["lead_bp"], loc["lead_snp"],
                    trait_id, gene_map)
                results.append(r)
            except Exception as e:
                results.append({"locus": f"chr{loc['chrom']}:{loc['lead_bp']}", "error": str(e)})

        n_ok = sum(1 for r in results if "hypothesis" in r)
        n_val = sum(1 for r in results if r.get("hypothesis", {}).get("is_validated"))
        print(f"  Done: {n_ok}/{len(results)} hypotheses, {n_val} validated")
        return results
