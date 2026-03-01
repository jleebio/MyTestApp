"""Smoke tests to verify project structure."""
import importlib

def test_package_imports():
    modules = [
        "gwas_loop",
        "gwas_loop.benchmark",
        "gwas_loop.failure_detection",
        "gwas_loop.hypothesis",
        "gwas_loop.prototyping",
        "gwas_loop.evaluation",
        "gwas_loop.memory",
    ]
    for mod in modules:
        importlib.import_module(mod)

def test_version():
    from gwas_loop import __version__
    assert __version__ == "0.1.0"
