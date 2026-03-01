"""Evaluate prototypes against baselines."""
import numpy as np
from sklearn.metrics import roc_auc_score, precision_score
from dataclasses import dataclass


@dataclass
class EvalResult:
    prototype_name: str
    trait: str
    precision_at_1: float
    precision_at_5: float
    auroc: float
    calibration_error: float
    promoted: bool


class Evaluator:
    """Evaluate and decide whether to promote a prototype."""

    def __init__(self, improvement_threshold: float = 0.05):
        self.threshold = improvement_threshold

    def evaluate(self, y_true: np.ndarray, y_scores: np.ndarray,
                 baseline_auroc: float, prototype_name: str, trait: str) -> EvalResult:
        auroc = roc_auc_score(y_true, y_scores) if len(np.unique(y_true)) > 1 else 0.0
        top1 = (y_true[np.argsort(y_scores)[-1:]] == 1).mean()
        top5 = (y_true[np.argsort(y_scores)[-5:]] == 1).mean()
        promoted = (auroc - baseline_auroc) >= self.threshold
        return EvalResult(
            prototype_name=prototype_name, trait=trait,
            precision_at_1=float(top1), precision_at_5=float(top5),
            auroc=auroc, calibration_error=0.0, promoted=promoted,
        )
