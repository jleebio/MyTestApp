"""Tests for the autonomous loop runner."""
from gwas_loop.runner import AutonomousRunner, LoopConfig


def test_runner_init(tmp_path):
    config = LoopConfig(
        data_dir=str(tmp_path / "data"),
        output_dir=str(tmp_path / "reports"),
        experiment_registry=str(tmp_path / "exp.json"),
        trait_registry=str(tmp_path / "traits.json"),
    )
    runner = AutonomousRunner(config)
    status = runner.status()
    assert len(status["traits"]) == 4
    assert status["experiments_run"] == 0


def test_runner_dry_run(tmp_path):
    config = LoopConfig(
        data_dir=str(tmp_path / "data"),
        output_dir=str(tmp_path / "reports"),
        experiment_registry=str(tmp_path / "exp.json"),
        trait_registry=str(tmp_path / "traits.json"),
        max_iterations=1,
    )
    runner = AutonomousRunner(config)
    results = runner.run(traits=["T2D"], dry_run=True)
    assert len(results) == 1
    assert results[0].trait == "T2D"


def test_runner_convergence(tmp_path):
    config = LoopConfig(
        data_dir=str(tmp_path / "data"),
        output_dir=str(tmp_path / "reports"),
        experiment_registry=str(tmp_path / "exp.json"),
        trait_registry=str(tmp_path / "traits.json"),
        max_iterations=5,
    )
    runner = AutonomousRunner(config)
    results = runner.run(traits=["T2D"])
    # Should converge after hypotheses are exhausted
    assert len(results) >= 1
