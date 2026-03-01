"""Smoke tests for all modules."""
import importlib

def test_all_imports():
    modules = [
        "gwas_loop", "gwas_loop.benchmark.runner",
        "gwas_loop.failure_detection.detector",
        "gwas_loop.hypothesis.generator",
        "gwas_loop.prototyping.builder",
        "gwas_loop.evaluation.evaluator",
        "gwas_loop.memory.registry",
        "gwas_loop.loop",
    ]
    for mod in modules:
        importlib.import_module(mod)

def test_version():
    from gwas_loop import __version__
    assert __version__ == "0.1.0"
