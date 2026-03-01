"""Generate benchmark tables for publication."""
import pandas as pd
from pathlib import Path


def method_comparison_table(results: list[dict], output_dir: str = "reports") -> pd.DataFrame:
    """Table 1: Method performance across traits.

    Columns: Method | T2D P@1 | T2D P@5 | T2D AUROC | CAD P@1 | ... | Mean
    """
    df = pd.DataFrame(results)
    if df.empty:
        return df

    pivot = df.pivot_table(
        index="method_name",
        columns="trait",
        values=["precision_at_1", "precision_at_5", "auroc"],
        aggfunc="mean",
    )
    pivot.columns = [f"{trait}_{metric}" for metric, trait in pivot.columns]

    # Add mean across traits
    for metric in ["precision_at_1", "precision_at_5", "auroc"]:
        cols = [c for c in pivot.columns if metric in c]
        pivot[f"mean_{metric}"] = pivot[cols].mean(axis=1)

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    pivot.round(3).to_csv(out / "table1_method_comparison.tsv", sep="\t")
    return pivot.round(3)


def innovation_impact_table(baseline: dict, innovation: dict, output_dir: str = "reports") -> pd.DataFrame:
    """Table 2: Impact of each innovation module.

    Shows: Baseline → +Disagreement → +Tissue → +Calibration → +All
    """
    rows = [{"configuration": "baseline", **baseline}]
    for name, metrics in innovation.items():
        rows.append({"configuration": f"+{name}", **metrics})

    df = pd.DataFrame(rows).set_index("configuration")
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    df.round(3).to_csv(out / "table2_innovation_impact.tsv", sep="\t")
    return df.round(3)


def therapeutic_hypothesis_table(hypotheses: list[dict], output_dir: str = "reports") -> pd.DataFrame:
    """Table 3: Top therapeutic hypotheses.

    Columns: Trait | Gene | Mechanism | Drug | Confidence | P(Trial Success)
    """
    df = pd.DataFrame(hypotheses)
    if df.empty:
        return df

    cols = ["trait", "causal_gene", "mechanism", "drug_candidate",
            "confidence_score", "confidence_tier"]
    available = [c for c in cols if c in df.columns]
    df = df[available].sort_values("confidence_score", ascending=False)

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "table3_therapeutic_hypotheses.tsv", sep="\t", index=False)
    return df


def cross_trait_transfer_table(results: list[dict], output_dir: str = "reports") -> pd.DataFrame:
    """Table 4: Cross-trait transfer learning results.

    Shows leave-one-trait-out AUROC and top predictions.
    """
    df = pd.DataFrame(results)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "table4_cross_trait_transfer.tsv", sep="\t", index=False)
    return df
