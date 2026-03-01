"""Network-based gene prioritization.

Uses protein-protein interaction and pathway proximity to
known disease genes for guilt-by-association scoring.
"""
import numpy as np
import pandas as pd
from collections import defaultdict
from .base import BaseMethod


class NetworkMethod(BaseMethod):
    """Score genes by network proximity to seed (known causal) genes.

    Uses a simple random walk with restart on a gene interaction graph.
    """
    name = "network"

    def __init__(self, restart_prob: float = 0.15, max_iter: int = 100, tol: float = 1e-6):
        self.restart_prob = restart_prob
        self.max_iter = max_iter
        self.tol = tol

    def _build_adjacency(self, edges: pd.DataFrame, genes: list[str]) -> np.ndarray:
        """Build normalized adjacency matrix from edge list."""
        gene_idx = {g: i for i, g in enumerate(genes)}
        n = len(genes)
        A = np.zeros((n, n))
        for _, row in edges.iterrows():
            g1, g2 = row["gene1"], row["gene2"]
            if g1 in gene_idx and g2 in gene_idx:
                w = row.get("weight", 1.0)
                A[gene_idx[g1], gene_idx[g2]] = w
                A[gene_idx[g2], gene_idx[g1]] = w
        # Row-normalize
        row_sums = A.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        return A / row_sums

    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        """Random walk with restart from seed genes.

        locus_data columns: gene_id, is_seed (bool)
        Also requires 'edges' in locus_data.attrs or a default empty graph.
        """
        genes = locus_data["gene_id"].tolist()
        seeds = locus_data[locus_data.get("is_seed", pd.Series(False, index=locus_data.index)) == True]["gene_id"].tolist()

        if not seeds or "edges" not in locus_data.attrs:
            return pd.Series(0.0, index=genes, name="network_score")

        edges = locus_data.attrs["edges"]
        A = self._build_adjacency(edges, genes)
        n = len(genes)

        # Restart vector: uniform over seeds
        r = np.zeros(n)
        seed_idx = [genes.index(s) for s in seeds if s in genes]
        if seed_idx:
            r[seed_idx] = 1.0 / len(seed_idx)

        # Random walk
        p = r.copy()
        for _ in range(self.max_iter):
            p_new = (1 - self.restart_prob) * A.T @ p + self.restart_prob * r
            if np.abs(p_new - p).sum() < self.tol:
                break
            p = p_new

        return pd.Series(p, index=genes, name="network_score")
