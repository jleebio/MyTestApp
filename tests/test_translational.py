"""Tests for translational pipeline modules."""
import numpy as np
import pandas as pd
from gwas_loop.translational.mechanism import MechanismInference, MechanismType
from gwas_loop.translational.tractability import TractabilityAnalyzer, TractabilityBucket, ProteinClass
from gwas_loop.translational.drug_mapping import DrugMapper, DrugGeneLink
from gwas_loop.translational.repurposing import RepurposingEngine
from gwas_loop.translational.confidence import GeneticTherapeuticScorer
from gwas_loop.translational.hypothesis_report import HypothesisReporter, TherapeuticHypothesis


def test_mechanism_coding():
    data = pd.DataFrame({
        "gene_id": ["TCF7L2", "TCF7L2"],
        "consequence": ["missense_variant", "intron_variant"],
        "eqtl_beta": [0.0, 0.0],
    })
    mi = MechanismInference()
    results = mi.infer(data)
    assert len(results) == 1
    assert results[0].mechanism == MechanismType.CODING_DISRUPTION


def test_mechanism_regulatory():
    data = pd.DataFrame({
        "gene_id": ["GENE1"] * 3,
        "consequence": ["intron_variant"] * 3,
        "eqtl_beta": [0.5, -0.3, 0.2],
        "tissue": ["brain", "liver", "blood"],
    })
    mi = MechanismInference()
    results = mi.infer(data)
    assert results[0].mechanism == MechanismType.EXPRESSION_REGULATION
    assert results[0].direction_of_effect == "gain"


def test_tractability_kinase():
    data = pd.DataFrame({
        "gene_id": ["JAK2"],
        "protein_class": ["kinase"],
        "has_structure": [True],
        "is_secreted": [False],
        "n_drugs": [3],
    })
    ta = TractabilityAnalyzer()
    results = ta.analyze(data)
    assert results[0].bucket == TractabilityBucket.CLINICAL
    assert results[0].score > 0.7


def test_tractability_intractable():
    data = pd.DataFrame({
        "gene_id": ["MYSTERY"],
        "protein_class": ["other"],
        "has_structure": [False],
        "is_secreted": [False],
        "n_drugs": [0],
    })
    ta = TractabilityAnalyzer()
    results = ta.analyze(data)
    assert results[0].bucket == TractabilityBucket.INTRACTABLE


def test_drug_mapping():
    drug_db = pd.DataFrame({
        "gene_id": ["JAK2", "JAK2", "PCSK9"],
        "drug_name": ["Ruxolitinib", "Fedratinib", "Evolocumab"],
        "drug_type": ["small_molecule", "small_molecule", "antibody"],
        "mechanism_of_action": ["inhibitor", "inhibitor", "inhibitor"],
        "max_phase": [4, 4, 4],
        "source": ["chembl"] * 3,
        "indication": ["myelofibrosis", "myelofibrosis", "hypercholesterolemia"],
    })
    mapper = DrugMapper()
    mapper.load_drug_database(drug_db)
    links = mapper.map_direct(["JAK2"])
    assert len(links) == 2
    assert all(l.gene_id == "JAK2" for l in links)


def test_repurposing_concordant():
    engine = RepurposingEngine()
    drugs = [{"drug_name": "Ruxolitinib", "mechanism_of_action": "inhibitor",
              "max_phase": 4, "indication": "myelofibrosis"}]
    hyps = engine.generate_hypotheses("JAK2", "CAD", "gain", drugs, causal_confidence=0.8)
    assert len(hyps) == 1
    assert hyps[0].concordance is True  # Gain-of-function + inhibitor = concordant
    assert hyps[0].score > 0.5


def test_repurposing_discordant():
    engine = RepurposingEngine()
    drugs = [{"drug_name": "SomeDrug", "mechanism_of_action": "agonist",
              "max_phase": 2, "indication": "other"}]
    hyps = engine.generate_hypotheses("GENE1", "T2D", "gain", drugs, causal_confidence=0.6)
    assert hyps[0].concordance is False  # Gain + agonist = wrong direction


def test_confidence_scorer():
    scorer = GeneticTherapeuticScorer()
    result = scorer.score(
        "PCSK9", "CAD",
        finemapping_pip=0.95, coloc_h4=0.88,
        direction_concordant=True, tissue_weight=0.9,
        pathway_score=0.7, lof_intolerance=0.95,
    )
    assert result.tier == "high"
    assert result.total_score > 0.7


def test_confidence_low():
    scorer = GeneticTherapeuticScorer()
    result = scorer.score("UNKNOWN", "T2D", finemapping_pip=0.1, coloc_h4=0.05)
    assert result.tier == "low"


def test_hypothesis_report(tmp_path):
    reporter = HypothesisReporter(output_dir=tmp_path)
    reporter.add(TherapeuticHypothesis(
        trait="CAD", locus_id="chr1:100000", causal_gene="PCSK9",
        mechanism="coding_disruption", direction_of_effect="gain",
        drug_candidate="Evolocumab", drug_action="inhibitor",
        repurposing_rationale="LoF protective, inhibitor should help",
        confidence_score=0.85, confidence_tier="high",
        tissue="liver", supporting_evidence=["PIP=0.95", "coloc_H4=0.88"],
    ))
    path = reporter.save_json()
    assert path.exists()
    path_tsv = reporter.save_tsv()
    assert path_tsv.exists()
    summary = reporter.summary()
    assert summary["total"] == 1
    assert summary["by_tier"]["high"] == 1
