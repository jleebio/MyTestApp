"""Generate publication figures.

Note: matplotlib is optional. These functions generate data for plots;
actual rendering can use matplotlib, plotly, or R/ggplot2.
"""
import json
from pathlib import Path


def figure1_pipeline_overview() -> dict:
    """Figure 1: Pipeline architecture diagram data."""
    return {
        "title": "Self-improving GWAS-to-therapeutics pipeline",
        "nodes": [
            {"id": "gwas", "label": "GWAS Summary\nStatistics", "layer": 0},
            {"id": "finemap", "label": "Fine-mapping\n& Variant QC", "layer": 1},
            {"id": "prioritize", "label": "Gene\nPrioritization\n(11 methods)", "layer": 2},
            {"id": "innovate", "label": "Innovation\nLayer\n(5 modules)", "layer": 2, "style": "highlight"},
            {"id": "mechanism", "label": "Mechanism\nInference", "layer": 3},
            {"id": "tractability", "label": "Target\nTractability", "layer": 4},
            {"id": "drugs", "label": "Drug Mapping\n& Repurposing", "layer": 5},
            {"id": "confidence", "label": "Confidence\nScoring", "layer": 6},
            {"id": "trial", "label": "Clinical Trial\nPrediction", "layer": 7},
            {"id": "feedback", "label": "Self-improvement\nFeedback", "layer": 8},
        ],
        "edges": [
            ("gwas", "finemap"), ("finemap", "prioritize"),
            ("prioritize", "innovate"), ("innovate", "prioritize"),
            ("prioritize", "mechanism"), ("mechanism", "tractability"),
            ("tractability", "drugs"), ("drugs", "confidence"),
            ("confidence", "trial"), ("trial", "feedback"),
            ("feedback", "prioritize"),
        ],
    }


def figure2_method_comparison_data(results: list[dict]) -> dict:
    """Figure 2: Heatmap data — methods × traits × metrics."""
    return {
        "title": "Gene prioritization performance across traits",
        "type": "heatmap",
        "data": results,
        "x_axis": "trait",
        "y_axis": "method_name",
        "value": "auroc",
    }


def figure3_disagreement_analysis_data(disagreement_scores: dict) -> dict:
    """Figure 3: Method disagreement reveals locus biology."""
    return {
        "title": "Disagreement index correlates with locus mechanism",
        "type": "scatter",
        "data": disagreement_scores,
        "x_axis": "disagreement_index",
        "y_axis": "mechanism_type",
    }


def figure4_repurposing_network_data(hypotheses: list[dict]) -> dict:
    """Figure 4: Drug repurposing network — traits connected by shared drug targets."""
    return {
        "title": "Cross-trait drug repurposing opportunities",
        "type": "network",
        "data": hypotheses,
    }


def save_figure_data(figure_data: dict, name: str, output_dir: str = "reports/figures"):
    """Save figure data as JSON for rendering."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{name}.json"
    path.write_text(json.dumps(figure_data, indent=2, default=str))
    return path
