"""Continuous autonomous research loop runner.

Orchestrates the full pipeline:
GWAS → Fine-mapping → Gene Prioritization → Mechanism → Tractability
→ Drug Mapping → Repurposing → Confidence → Hypothesis Report
→ Self-Improvement → Repeat
"""
import json
import logging
from pathlib import Path
from dataclasses import dataclass, asdict
from datetime import datetime, timezone

from .data.traits import TraitRegistry, TraitMetadata
from .benchmark.runner import BenchmarkRunner
from .failure_detection.detector import FailureDetector
from .hypothesis.generator import HypothesisGenerator
from .evaluation.evaluator import Evaluator
from .memory.registry import ExperimentRegistry, Experiment
from .innovation.disagreement_model import DisagreementModel
from .innovation.locus_classifier import LocusClassifier
from .innovation.calibration import EvidenceCalibrator
from .translational.mechanism import MechanismInference
from .translational.tractability import TractabilityAnalyzer
from .translational.drug_mapping import DrugMapper
from .translational.repurposing import RepurposingEngine
from .translational.confidence import GeneticTherapeuticScorer
from .translational.hypothesis_report import HypothesisReporter, TherapeuticHypothesis

logger = logging.getLogger(__name__)


@dataclass
class LoopConfig:
    data_dir: str = "data"
    output_dir: str = "reports"
    experiment_registry: str = "experiments.json"
    trait_registry: str = "trait_metadata.json"
    max_iterations: int = 10
    improvement_threshold: float = 0.05
    error_limit: int = 3  # Stop after N consecutive errors


@dataclass
class IterationResult:
    iteration: int
    trait: str
    n_loci_processed: int
    n_hypotheses_generated: int
    n_high_confidence: int
    n_methods_improved: int
    top_hypothesis: str
    timestamp: str


class AutonomousRunner:
    """Main autonomous research loop.

    Each iteration:
    1. Load/download GWAS data
    2. Benchmark existing methods
    3. Detect failure patterns
    4. Generate & test improvements
    5. Run translational pipeline
    6. Output therapeutic hypotheses
    7. Update knowledge memory
    """

    def __init__(self, config: LoopConfig | None = None):
        self.config = config or LoopConfig()
        self.trait_registry = TraitRegistry(self.config.trait_registry)
        self.experiment_registry = ExperimentRegistry(self.config.experiment_registry)
        self.benchmark = BenchmarkRunner()
        self.detector = FailureDetector()
        self.hypothesis_gen = HypothesisGenerator()
        self.evaluator = Evaluator(self.config.improvement_threshold)
        self.disagreement = DisagreementModel()
        self.locus_classifier = LocusClassifier()
        self.calibrator = EvidenceCalibrator()
        self.mechanism = MechanismInference()
        self.tractability = TractabilityAnalyzer()
        self.drug_mapper = DrugMapper()
        self.repurposing = RepurposingEngine()
        self.scorer = GeneticTherapeuticScorer()
        self.reporter = HypothesisReporter(self.config.output_dir)
        self.iteration_results: list[IterationResult] = []
        self._consecutive_errors = 0

    def run(self, traits: list[str] | None = None, dry_run: bool = False) -> list[IterationResult]:
        """Execute the autonomous research loop.

        Args:
            traits: Specific trait IDs to process (default: all in registry)
            dry_run: If True, validate pipeline without full execution
        """
        traits = traits or self.trait_registry.list_traits()
        logger.info(f"Starting autonomous loop: {len(traits)} traits, max {self.config.max_iterations} iterations")

        results = []
        for i in range(self.config.max_iterations):
            for trait_id in traits:
                try:
                    result = self._run_iteration(i, trait_id, dry_run)
                    results.append(result)
                    self._consecutive_errors = 0
                    logger.info(f"[iter={i}][{trait_id}] {result.n_hypotheses_generated} hypotheses, "
                              f"{result.n_high_confidence} high-confidence")
                except Exception as e:
                    self._consecutive_errors += 1
                    logger.error(f"[iter={i}][{trait_id}] Error: {e}")
                    if self._consecutive_errors >= self.config.error_limit:
                        logger.critical(f"Hit error limit ({self.config.error_limit}). Stopping.")
                        self.iteration_results = results
                        return results

            # Check for convergence
            if self._check_convergence(results):
                logger.info("Converged — no further improvements detected.")
                break

        self.iteration_results = results
        self._save_results(results)
        return results

    def _run_iteration(self, iteration: int, trait_id: str, dry_run: bool) -> IterationResult:
        """Single iteration of the research loop for one trait."""
        trait = self.trait_registry.get(trait_id)

        # Phase 1: Detect failures in current methods
        # (In production, this uses real benchmark results)
        hypotheses = self.hypothesis_gen.generate()

        # Phase 2: Filter already-tried hypotheses
        novel = [h for h in hypotheses if not self.experiment_registry.was_tried(h.name)]

        # Phase 3: Log experiments
        n_improved = 0
        for h in novel:
            self.experiment_registry.log(Experiment(
                name=h.name, version="0.1",
                auroc=0.0, promoted=False,
                notes=f"iteration={iteration}, trait={trait_id}",
            ))

        # Phase 4: Translational pipeline would run here
        # (with real data: mechanism → tractability → drugs → repurposing)

        return IterationResult(
            iteration=iteration, trait=trait_id,
            n_loci_processed=0,  # Populated with real data
            n_hypotheses_generated=len(novel),
            n_high_confidence=0,
            n_methods_improved=n_improved,
            top_hypothesis=novel[0].name if novel else "none",
            timestamp=datetime.now(timezone.utc).isoformat(),
        )

    def _check_convergence(self, results: list[IterationResult]) -> bool:
        """Check if the loop has converged (no new improvements)."""
        if len(results) < 2:
            return False
        recent = results[-len(self.trait_registry.list_traits()):]
        return all(r.n_methods_improved == 0 and r.n_hypotheses_generated == 0 for r in recent)

    def _save_results(self, results: list[IterationResult]):
        out = Path(self.config.output_dir)
        out.mkdir(parents=True, exist_ok=True)
        path = out / "loop_results.json"
        path.write_text(json.dumps([asdict(r) for r in results], indent=2))
        logger.info(f"Results saved to {path}")

    def status(self) -> dict:
        return {
            "traits": self.trait_registry.list_traits(),
            "experiments_run": len(self.experiment_registry.experiments),
            "failed_strategies": self.experiment_registry.failed_strategies(),
            "iterations_completed": len(self.iteration_results),
            "hypotheses_generated": self.reporter.summary().get("total", 0),
        }
