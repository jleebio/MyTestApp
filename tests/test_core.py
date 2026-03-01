"""Core functionality tests."""
import numpy as np
import pandas as pd
from gwas_loop.failure_detection.detector import FailureDetector
from gwas_loop.hypothesis.generator import HypothesisGenerator
from gwas_loop.prototyping.builder import PrototypeBuilder
from gwas_loop.evaluation.evaluator import Evaluator
from gwas_loop.memory.registry import ExperimentRegistry, Experiment
from gwas_loop.loop import ResearchLoop

def test_failure_detector():
    preds = {
        "method_a": pd.Series([0.9, 0.1, 0.5], index=["G1", "G2", "G3"]),
        "method_b": pd.Series([0.1, 0.9, 0.5], index=["G1", "G2", "G3"]),
    }
    det = FailureDetector()
    results = det.detect(preds)
    assert len(results) == 3
    assert all(r.classification in ("easy", "moderate", "hard") for r in results)

def test_hypothesis_generator():
    gen = HypothesisGenerator()
    hyps = gen.generate()
    assert len(hyps) >= 3
    assert all(h.status == "proposed" for h in hyps)

def test_prototype_builder():
    builder = PrototypeBuilder()
    proto = builder.build("test_hypothesis")
    assert not proto.trained
    X = np.random.rand(50, 5)
    y = (X[:, 0] > 0.5).astype(int)
    proto = builder.train(proto, X, y)
    assert proto.trained

def test_evaluator():
    ev = Evaluator(improvement_threshold=0.05)
    y_true = np.array([0, 0, 1, 1, 0, 1, 0, 0, 1, 0])
    y_scores = np.array([0.1, 0.2, 0.8, 0.7, 0.3, 0.9, 0.2, 0.15, 0.6, 0.1])
    result = ev.evaluate(y_true, y_scores, baseline_auroc=0.5, prototype_name="test", trait="T2D")
    assert result.auroc > 0.5
    assert result.promoted

def test_registry(tmp_path):
    reg = ExperimentRegistry(tmp_path / "exp.json")
    reg.log(Experiment(name="test1", version="0.1", auroc=0.8, promoted=True))
    reg.log(Experiment(name="bad_idea", version="0.1", auroc=0.4, promoted=False))
    assert reg.was_tried("test1")
    assert not reg.was_tried("never_tried")
    assert "bad_idea" in reg.failed_strategies()

def test_research_loop(tmp_path):
    loop = ResearchLoop(registry_path=str(tmp_path / "exp.json"))
    status = loop.status()
    assert status["methods_registered"] == 0
    assert status["experiments_total"] == 0
