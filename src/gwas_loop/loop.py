"""Main self-improving research loop orchestrator."""
from .benchmark.runner import BenchmarkRunner
from .failure_detection.detector import FailureDetector
from .hypothesis.generator import HypothesisGenerator
from .prototyping.builder import PrototypeBuilder
from .evaluation.evaluator import Evaluator
from .memory.registry import ExperimentRegistry, Experiment


class ResearchLoop:
    """Orchestrate the autonomous GWAS research improvement cycle."""

    def __init__(self, registry_path: str = "experiments.json"):
        self.benchmark = BenchmarkRunner()
        self.detector = FailureDetector()
        self.generator = HypothesisGenerator()
        self.builder = PrototypeBuilder()
        self.evaluator = Evaluator()
        self.registry = ExperimentRegistry(registry_path)

    def status(self) -> dict:
        return {
            "methods_registered": len(self.benchmark.methods),
            "experiments_total": len(self.registry.experiments),
            "failed_strategies": self.registry.failed_strategies(),
        }
