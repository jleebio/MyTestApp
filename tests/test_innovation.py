"""Tests for novel innovation modules."""
import numpy as np
import pandas as pd
from gwas_loop.innovation.disagreement_model import DisagreementModel
from gwas_loop.innovation.tissue_adaptive import TissueAdaptiveModel
from gwas_loop.innovation.locus_classifier import LocusClassifier, LocusMechanism
from gwas_loop.innovation.calibration import EvidenceCalibrator
from gwas_loop.innovation.cross_trait import CrossTraitTransfer


def test_disagreement_features():
    genes = [f"G{i}" for i in range(20)]
    scores = {
        "distance": pd.Series(np.random.rand(20), index=genes),
        "twas": pd.Series(np.random.rand(20), index=genes),
        "coloc": pd.Series(np.random.rand(20), index=genes),
        "vep": pd.Series(np.random.rand(20), index=genes),
    }
    dm = DisagreementModel()
    feats = dm.extract_features(scores)
    assert feats.shape == (20, 10)
    assert "coding_vs_regulatory_discord" in feats.columns


def test_disagreement_train_predict():
    n = 50
    genes = [f"G{i}" for i in range(n)]
    scores = {
        "distance": pd.Series(np.random.rand(n), index=genes),
        "twas": pd.Series(np.random.rand(n), index=genes),
        "coloc": pd.Series(np.random.rand(n), index=genes),
    }
    labels = pd.Series(np.random.randint(0, 2, n), index=genes)
    dm = DisagreementModel()
    dm.train(scores, labels)
    preds = dm.predict(scores)
    assert len(preds) == n
    assert dm.feature_importance().shape[0] == 10


def test_tissue_adaptive():
    genes = [f"G{i}" for i in range(30)]
    tissues = ["brain", "liver", "pancreas", "blood"]
    tissue_scores = pd.DataFrame(
        np.random.rand(30, 4), index=genes, columns=tissues
    )
    labels = pd.Series(np.random.randint(0, 2, 30), index=genes)
    tam = TissueAdaptiveModel(n_tissues=4)
    tam.train(tissue_scores, labels, trait="T2D")
    weights = tam.get_weights("T2D")
    assert len(weights) == 4
    scores = tam.score(tissue_scores, "T2D")
    assert len(scores) == 30


def test_locus_classifier_heuristic():
    locus = pd.DataFrame({
        "gene_id": ["A", "B", "C"],
        "has_missense": [True, False, False],
    })
    lc = LocusClassifier()
    result = lc.classify(locus, "locus1")
    assert result.mechanism == LocusMechanism.CODING
    assert "vep" in result.recommended_methods


def test_calibrator():
    cal = EvidenceCalibrator()
    scores = np.random.rand(100)
    labels = (scores > 0.5).astype(int)  # Perfectly calibrated
    cal.fit("twas", scores, labels)
    calibrated = cal.calibrate("twas", scores)
    assert len(calibrated) == 100
    report = cal.summary()
    assert "twas" in report["method"].values


def test_cross_trait():
    n_per_trait = 30
    traits = []
    features_list = []
    labels_list = []
    for trait in ["T2D", "CAD", "SCZ"]:
        feats = pd.DataFrame(
            np.random.rand(n_per_trait, 5),
            index=[f"{trait}_{i}" for i in range(n_per_trait)],
            columns=[f"f{j}" for j in range(5)],
        )
        labs = pd.Series(np.random.randint(0, 2, n_per_trait), index=feats.index)
        features_list.append(feats)
        labels_list.append(labs)
        traits.extend([trait] * n_per_trait)

    all_feats = pd.concat(features_list)
    all_labels = pd.concat(labels_list)
    all_traits = pd.Series(traits, index=all_feats.index)

    ct = CrossTraitTransfer()
    results = ct.train_leave_one_trait_out(all_feats, all_labels, all_traits)
    assert len(results) == 3
    assert all(r.target_trait in ["T2D", "CAD", "SCZ"] for r in results)
