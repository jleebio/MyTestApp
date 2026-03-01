"""Benchmark runner for GWAS gene prioritization methods."""
import pandas as pd
import numpy as np
from dataclasses import dataclass, field
from typing import Protocol, Any


class PrioritizationMethod(Protocol):
    """Interface for any gene prioritization method."""
    name: str
    def score(self, locus_data: pd.DataFrame) -> pd.Series: ...


@dataclass
class BenchmarkResult:
    method_name: str
    trait: str
    precision_at_1: float
    precision_at_5: float
    auroc: float
    calibration_error: float
    validated_recovery: float
    predictions: pd.DataFrame = field(default_factory=pd.DataFrame)


class BenchmarkRunner:
    """Run and compare multiple prioritization methods on GWAS loci."""

    def __init__(self, methods: list[PrioritizationMethod] | None = None):
        self.methods = methods or []
        self.results: list[BenchmarkResult] = []

    def add_method(self, method: PrioritizationMethod):
        self.methods.append(method)

    def run(self, loci: pd.DataFrame, validated_genes: set[str], trait: str) -> list[BenchmarkResult]:
        results = []
        for method in self.methods:
            scores = method.score(loci)
            ranked = scores.sort_values(ascending=False)
            top1 = set(ranked.head(1).index)
            top5 = set(ranked.head(5).index)
            p1 = len(top1 & validated_genes) / max(len(top1), 1)
            p5 = len(top5 & validated_genes) / max(len(top5), 1)
            result = BenchmarkResult(
                method_name=method.name, trait=trait,
                precision_at_1=p1, precision_at_5=p5,
                auroc=0.0, calibration_error=0.0,
                validated_recovery=len(set(ranked.index) & validated_genes) / max(len(validated_genes), 1),
                predictions=ranked.to_frame("score"),
            )
            results.append(result)
        self.results.extend(results)
        return results
