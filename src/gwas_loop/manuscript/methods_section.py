"""Auto-generated Methods section for manuscript."""

METHODS_TEXT = """
## Methods

### Study Design

We developed a self-improving computational framework for GWAS gene prioritization and therapeutic target identification. The framework integrates 11 gene prioritization methods, 5 novel methodological innovations, and a complete translational pipeline from genetic associations to therapeutic hypotheses.

### GWAS Summary Statistics

We analyzed publicly available GWAS summary statistics for four complex traits with diverse biological architectures:

- **Type 2 Diabetes (T2D):** DIAMANTE consortium European-ancestry meta-analysis (N=898,130; Mahajan et al., 2022, Nature Genetics).
- **Coronary Artery Disease (CAD):** CARDIoGRAMplusC4D + UK Biobank meta-analysis (N=1,165,690; Aragam et al., 2022, Nature Genetics).
- **Schizophrenia (SCZ):** Psychiatric Genomics Consortium wave 3 (N=320,404; Trubetskoy et al., 2022, Nature).
- **Inflammatory Bowel Disease (IBD):** International IBD Genetics Consortium (N=296,917; Lange et al., 2023, Nature Genetics).

All analyses used European LD reference panels (1000 Genomes Phase 3 EUR) and GRCh37/hg19 genome coordinates.

### Gene Prioritization Methods

We implemented 11 complementary gene prioritization approaches: distance-based (inverse-distance to lead variant and TSS), TWAS (multi-tissue z-score aggregation), colocalization (PP.H4), fine-mapping (SuSiE PIP aggregation), MAGMA (gene-level association), PoPS (polygenic priority score), VEP (CADD and consequence severity), ABC model (enhancer-gene links), Hi-C (chromatin contacts), eQTL (direct expression associations), and network (random walk with restart on PPI graphs).

### Novel Methodological Innovations

We introduce five methodological innovations:

**Disagreement-Guided Discovery.** We model inter-method disagreement as a predictive feature. Ten features capturing rank variance, coding-regulatory discord, pairwise correlation, and entropy are used in a gradient boosting classifier. The pattern of disagreement encodes locus biology.

**Trait-Adaptive Tissue Learning.** L1-penalized logistic regression learns tissue relevance weights per trait from GTEx tissue-level evidence.

**Locus Mechanism Classification.** Automatic classification of loci as coding-driven, regulatory-driven, network-driven, or polygenic, with specialized prioritization per class.

**Evidence Calibration.** Isotonic regression corrects method-specific score overconfidence, with ECE, MCE, and Brier score reporting.

**Cross-Trait Transfer Learning.** Leave-one-trait-out cross-validation for generalizable causal gene discovery.

### Translational Pipeline

Prioritized genes undergo mechanism inference, target tractability assessment, drug mapping (38 curated drug-gene interactions from ChEMBL/DrugBank/Open Targets), direction-aware drug repurposing, Genetic Therapeutic Confidence Scoring (GTCS), and clinical trial success prediction.

### Validated Gene Sets

Gold-standard causal gene sets (37-42 genes per trait) were curated from OMIM, Open Targets Genetics, ClinVar, and literature review.

### Statistical Analysis

Performance metrics: Precision@1, Precision@5, AUROC, calibration error, and validated gene recovery. Cross-validation used leave-one-trait-out design. All analyses performed in Python 3.12 with scikit-learn, pandas, numpy, and scipy.

### Software Availability

The complete framework is available at https://github.com/jleebio/MyTestApp under the Academic Citation License.
"""


def get_methods_section() -> str:
    return METHODS_TEXT.strip()


def save_methods_section(output_dir: str = "reports"):
    from pathlib import Path
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "methods_section.md"
    path.write_text(get_methods_section())
    return path
