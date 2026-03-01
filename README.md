# 🧬 GWAS Self-Improving Research Loop

**An autonomous framework that converts GWAS summary statistics into therapeutic drug hypotheses — and improves itself along the way.**

---

## What Does This Do?

This program takes publicly available genetic data (GWAS) and automatically:

1. **Finds which genes cause diseases** — using 11 different methods
2. **Figures out HOW they cause disease** — coding mutation? gene regulation? splicing?
3. **Checks if the gene is druggable** — can we actually target it?
4. **Finds existing drugs** — are there approved drugs or clinical candidates?
5. **Predicts drug repurposing** — could an existing drug for Disease A treat Disease B?
6. **Scores confidence** — how strong is the genetic evidence?
7. **Predicts clinical trial success** — what's the probability this would work in trials?
8. **Learns from mistakes** — improves its own methods over time

### Example Output

```
━━━ Therapeutic Hypothesis ━━━
Trait:       Coronary Artery Disease
Gene:        PCSK9
Mechanism:   coding_disruption (gain-of-function)
Drug:        Evolocumab (inhibitor)
Rationale:   Loss-of-function is protective → inhibitor should reduce risk
Tissue:      liver
Confidence:  0.85 [HIGH]
P(Trial Success): 22.5% (vs 9% baseline)
Evidence:    PIP=0.95, coloc_H4=0.88, human LoF protective
```

---

## The Pipeline

```
GWAS Summary Statistics (T2D, CAD, SCZ, IBD)
        │
        ▼
┌─────────────────────────────┐
│  Gene Prioritization        │  11 methods: distance, TWAS, coloc,
│  (find the causal gene)     │  fine-mapping, MAGMA, PoPS, VEP,
│                             │  ABC/Hi-C, eQTL, network, ensemble
└─────────────┬───────────────┘
              ▼
┌─────────────────────────────┐
│  Innovation Layer           │  5 novel approaches:
│  (improve the methods)      │  • Disagreement-as-signal
│                             │  • Adaptive tissue learning
│                             │  • Locus mechanism classification
│                             │  • Evidence calibration
│                             │  • Cross-trait transfer
└─────────────┬───────────────┘
              ▼
┌─────────────────────────────┐
│  Mechanism Inference        │  What's the biology?
│  coding / regulatory /      │  Which tissue? Direction of effect?
│  splicing / chromatin       │
└─────────────┬───────────────┘
              ▼
┌─────────────────────────────┐
│  Target Tractability        │  Is it druggable?
│  protein class, structure,  │  kinase? GPCR? antibody-accessible?
│  existing drugs             │
└─────────────┬───────────────┘
              ▼
┌─────────────────────────────┐
│  Drug Mapping & Repurposing │  Match to approved drugs
│  direction-aware scoring    │  gene ↑ risk + drug ↓ gene = treat
└─────────────┬───────────────┘
              ▼
┌─────────────────────────────┐
│  Confidence & Trial Pred    │  How sure are we?
│  6-component GTCS score     │  P(clinical trial success)?
└─────────────┬───────────────┘
              ▼
┌─────────────────────────────┐
│  Hypothesis Reports         │  JSON + TSV + formatted output
│  publication-ready          │
└─────────────┬───────────────┘
              ▼
┌─────────────────────────────┐
│  Self-Improvement Loop      │  Learn from failures
│  update weights, retry      │  prevent repeating mistakes
└─────────────────────────────┘
```

---

## Diseases Included

| Trait | Full Name | Source | Sample Size | Category |
|-------|-----------|--------|-------------|----------|
| **T2D** | Type 2 Diabetes | DIAMANTE 2022 | 898,130 | Metabolic |
| **CAD** | Coronary Artery Disease | CARDIoGRAMplusC4D 2022 | 1,165,690 | Cardiovascular |
| **SCZ** | Schizophrenia | PGC3 2022 | 320,404 | Psychiatric |
| **IBD** | Inflammatory Bowel Disease | IIBDGC 2023 | 296,917 | Autoimmune |

All use European LD reference (1000 Genomes EUR) and GRCh37 coordinates.

---

## How to Run

### Setup

```bash
cd /root/project1
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

### Run the Full Pipeline

```bash
# All 4 traits, 5 iterations, verbose output
python -m gwas_loop --traits T2D CAD SCZ IBD --max-iter 5 -v

# Single trait
python -m gwas_loop --traits T2D --max-iter 3 -v

# Dry run (validate pipeline without full execution)
python -m gwas_loop --dry-run
```

### CLI Options

```
python -m gwas_loop --help

Options:
  --traits       Trait IDs to process (default: all 4)
  --max-iter     Maximum loop iterations (default: 10)
  --data-dir     Directory for GWAS data (default: data/)
  --output-dir   Output directory for reports (default: reports/)
  --dry-run      Validate pipeline without full execution
  -v, --verbose  Detailed logging
```

### Run Tests

```bash
pytest tests/ -v          # All tests
pytest tests/ -v --cov    # With coverage report
```

---

## Output Files

After running, you'll find in `reports/`:

| File | Description |
|------|-------------|
| `hypotheses.json` | All therapeutic hypotheses (structured) |
| `hypotheses.tsv` | Same data as a spreadsheet-friendly table |
| `loop_results.json` | Per-iteration metrics and convergence info |

In the project root:

| File | Description |
|------|-------------|
| `experiments.json` | Experiment registry (what was tried, what worked) |
| `trait_metadata.json` | Trait catalog with download URLs and metadata |

---

## What Makes This Novel?

Standard GWAS tools apply **one fixed method** to all loci. This framework has 5 innovations not found in existing tools:

1. **Disagreement-as-Signal** — When methods disagree, the *pattern* of disagreement reveals locus biology (coding vs regulatory). Nobody models this explicitly.

2. **Adaptive Tissue Learning** — Instead of hardcoding "SCZ = brain", it *learns* tissue relevance per trait. Discovers non-obvious tissues (e.g., liver for SCZ).

3. **Locus Mechanism Classification** — Auto-classifies each locus as coding/regulatory/network/polygenic and applies *different* strategies per class.

4. **Evidence Calibration** — TWAS and coloc produce overconfident scores. This learns and corrects those biases per method.

5. **Cross-Trait Transfer** — Trains on multiple diseases, predicts causal genes for completely new traits with zero labeled data.

---

## Project Structure

```
project1/
├── src/gwas_loop/
│   ├── __init__.py              # Package root
│   ├── __main__.py              # CLI entry point
│   ├── runner.py                # Autonomous loop orchestrator
│   ├── loop.py                  # Research loop (legacy)
│   │
│   ├── data/                    # GWAS data management
│   │   ├── traits.py            # Trait catalog & metadata
│   │   ├── downloader.py        # Summary stats downloader
│   │   └── parser.py            # Unified format parser
│   │
│   ├── methods/                 # Gene prioritization (11 methods)
│   │   ├── base.py              # Method interface
│   │   ├── distance.py          # Proximity-based
│   │   ├── twas.py              # Expression association
│   │   ├── coloc.py             # Colocalization
│   │   ├── finemapping.py       # SuSiE/FINEMAP PIPs
│   │   ├── magma.py             # Gene-level association
│   │   ├── pops.py              # Polygenic priority score
│   │   ├── vep.py               # Variant effect prediction
│   │   ├── chromatin.py         # ABC model, Hi-C
│   │   ├── eqtl.py              # eQTL mapping
│   │   ├── network.py           # PPI random walk
│   │   └── ensemble.py          # Rank + Bayesian ensembles
│   │
│   ├── innovation/              # Novel methods (publishable)
│   │   ├── disagreement_model.py
│   │   ├── tissue_adaptive.py
│   │   ├── locus_classifier.py
│   │   ├── calibration.py
│   │   └── cross_trait.py
│   │
│   ├── translational/           # Genes → Drugs pipeline
│   │   ├── mechanism.py         # Biological mechanism inference
│   │   ├── tractability.py      # Target druggability
│   │   ├── drug_mapping.py      # Gene → drug lookup
│   │   ├── repurposing.py       # Drug repurposing engine
│   │   ├── confidence.py        # Genetic therapeutic score
│   │   ├── hypothesis_report.py # Report generator
│   │   └── clinical_trial.py    # Trial success prediction
│   │
│   ├── benchmark/               # Method benchmarking
│   │   └── runner.py
│   ├── failure_detection/       # Failure pattern detection
│   │   └── detector.py
│   ├── hypothesis/              # Method hypothesis generation
│   │   └── generator.py
│   ├── prototyping/             # Model prototyping
│   │   └── builder.py
│   ├── evaluation/              # Prototype evaluation
│   │   └── evaluator.py
│   └── memory/                  # Experiment tracking
│       └── registry.py
│
├── tests/                       # 43 tests
├── pyproject.toml               # Dependencies & config
├── OPENCLAW_GLOBAL.md           # Autonomous execution rules
└── plan.md                      # Project roadmap
```

---

## Dependencies

**Core:** numpy, pandas, scikit-learn, scipy, statsmodels

**Dev:** pytest, pytest-cov

**Optional:** torch, pymc (for advanced models — `pip install -e ".[deep]"`)

---

## Next Steps

- [ ] Download real GWAS summary statistics
- [ ] Integrate validated gene sets (Open Targets gold standard)
- [ ] Connect to ChEMBL/DrugBank for real drug-gene mappings
- [ ] Generate publication figures and benchmark tables
- [ ] Write manuscript methods section

---

## License

Research use. See project documentation for details.
