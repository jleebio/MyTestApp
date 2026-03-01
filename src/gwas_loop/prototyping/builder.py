"""Convert hypotheses into executable model prototypes."""
import numpy as np
from sklearn.ensemble import GradientBoostingClassifier
from dataclasses import dataclass, field
from typing import Any


@dataclass
class Prototype:
    hypothesis_name: str
    version: str
    model: Any = None
    trained: bool = False


class PrototypeBuilder:
    """Build sklearn-based prototypes from hypotheses."""

    def build(self, hypothesis_name: str, version: str = "0.1") -> Prototype:
        model = GradientBoostingClassifier(
            n_estimators=100, max_depth=4, learning_rate=0.1, random_state=42
        )
        return Prototype(hypothesis_name=hypothesis_name, version=version, model=model)

    def train(self, prototype: Prototype, X: np.ndarray, y: np.ndarray) -> Prototype:
        prototype.model.fit(X, y)
        prototype.trained = True
        return prototype
