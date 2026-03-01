"""Tests for figure and table generation."""
from gwas_loop.figures.benchmark_tables import method_comparison_table, therapeutic_hypothesis_table
from gwas_loop.figures.plots import figure1_pipeline_overview, save_figure_data


def test_method_comparison_empty():
    df = method_comparison_table([])
    assert df.empty

def test_method_comparison(tmp_path):
    results = [
        {"method_name": "distance", "trait": "T2D", "precision_at_1": 0.35, "precision_at_5": 0.20, "auroc": 0.65},
        {"method_name": "twas", "trait": "T2D", "precision_at_1": 0.45, "precision_at_5": 0.30, "auroc": 0.72},
        {"method_name": "distance", "trait": "CAD", "precision_at_1": 0.30, "precision_at_5": 0.18, "auroc": 0.60},
        {"method_name": "twas", "trait": "CAD", "precision_at_1": 0.50, "precision_at_5": 0.35, "auroc": 0.75},
    ]
    df = method_comparison_table(results, str(tmp_path))
    assert "distance" in df.index
    assert "twas" in df.index

def test_figure1_data():
    data = figure1_pipeline_overview()
    assert len(data["nodes"]) >= 8
    assert len(data["edges"]) >= 8

def test_save_figure(tmp_path):
    data = figure1_pipeline_overview()
    path = save_figure_data(data, "fig1", str(tmp_path))
    assert path.exists()
