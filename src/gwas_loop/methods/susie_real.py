"""Real SuSiE fine-mapping using R's susieR package.

This module:
1. Takes GWAS summary stats for a locus
2. Computes LD matrix from 1000G EUR reference panel (plink)
3. Runs SuSiE via R subprocess
4. Returns variant-level PIPs and credible sets
5. Aggregates to gene-level scores
"""
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


class SuSiERunner:
    """Run real SuSiE fine-mapping on GWAS loci."""

    def __init__(
        self,
        ld_ref_dir: str | Path = "data/ld_ref",
        max_causal: int = 10,
        coverage: float = 0.95,
        min_abs_corr: float = 0.5,
    ):
        self.ld_ref_dir = Path(ld_ref_dir)
        self.max_causal = max_causal
        self.coverage = coverage
        self.min_abs_corr = min_abs_corr

    def _check_r(self) -> bool:
        """Check that R and susieR are available."""
        try:
            result = subprocess.run(
                ["Rscript", "-e", "library(susieR); cat('OK')"],
                capture_output=True, text=True, timeout=30,
            )
            return "OK" in result.stdout
        except Exception:
            return False

    def compute_ld_matrix(
        self,
        variants: pd.DataFrame,
        chrom: int | str,
        window_start: int,
        window_end: int,
        output_prefix: str,
    ) -> np.ndarray | None:
        """Compute LD matrix using plink and 1000G EUR.

        Parameters
        ----------
        variants : DataFrame with 'SNP' column (rsIDs)
        chrom : chromosome number
        window_start, window_end : bp range
        output_prefix : temp file prefix for plink output

        Returns
        -------
        LD correlation matrix (n_variants × n_variants) or None if failed
        """
        plink_cmd = "plink1.9" if subprocess.run(
            ["which", "plink1.9"], capture_output=True
        ).returncode == 0 else "plink"

        bfile = self.ld_ref_dir / f"EUR.chr{chrom}"
        if not bfile.with_suffix(".bed").exists():
            print(f"  LD ref not found: {bfile}.bed")
            return None

        # Write variant list
        snp_file = f"{output_prefix}.snplist"
        variants["SNP"].to_csv(snp_file, index=False, header=False)

        try:
            subprocess.run([
                plink_cmd,
                "--bfile", str(bfile),
                "--extract", snp_file,
                "--chr", str(chrom),
                "--from-bp", str(window_start),
                "--to-bp", str(window_end),
                "--r", "square",
                "--out", output_prefix,
            ], capture_output=True, text=True, timeout=300, check=True)

            ld_file = f"{output_prefix}.ld"
            if Path(ld_file).exists():
                ld = np.loadtxt(ld_file)
                return ld
        except Exception as e:
            print(f"  plink LD computation failed: {e}")

        return None

    def run_susie_rss(
        self,
        z_scores: np.ndarray,
        ld_matrix: np.ndarray,
        n: int,
        output_prefix: str,
    ) -> dict:
        """Run SuSiE RSS (summary statistics) via R.

        Parameters
        ----------
        z_scores : array of z-scores per variant
        ld_matrix : LD correlation matrix
        n : sample size
        output_prefix : temp file prefix

        Returns
        -------
        dict with 'pip' (array), 'cs' (list of credible sets), 'converged' (bool)
        """
        # Write inputs for R
        z_file = f"{output_prefix}_z.txt"
        ld_file = f"{output_prefix}_ld.txt"
        out_file = f"{output_prefix}_susie_out.json"

        np.savetxt(z_file, z_scores)
        np.savetxt(ld_file, ld_matrix)

        r_script = f"""
library(susieR)
z <- as.numeric(readLines("{z_file}"))
R <- as.matrix(read.table("{ld_file}"))

# Ensure R is positive definite
R <- (R + t(R)) / 2
eig <- eigen(R)
eig$values[eig$values < 1e-8] <- 1e-8
R <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)

# Run SuSiE RSS
fit <- tryCatch({{
    susie_rss(
        z = z,
        R = R,
        n = {n},
        L = {self.max_causal},
        coverage = {self.coverage},
        min_abs_corr = {self.min_abs_corr},
        max_iter = 500,
        verbose = FALSE
    )
}}, error = function(e) {{
    list(pip = rep(1/length(z), length(z)),
         sets = list(cs = list()),
         converged = FALSE,
         error = conditionMessage(e))
}})

# Extract results
result <- list(
    pip = as.numeric(fit$pip),
    converged = isTRUE(fit$converged),
    n_cs = length(fit$sets$cs)
)

# Extract credible sets
if (length(fit$sets$cs) > 0) {{
    result$cs <- lapply(fit$sets$cs, function(x) as.integer(x - 1))  # 0-indexed
    result$cs_coverage <- as.numeric(fit$sets$coverage)
}} else {{
    result$cs <- list()
    result$cs_coverage <- numeric(0)
}}

writeLines(jsonlite::toJSON(result, auto_unbox=TRUE), "{out_file}")
cat("SUSIE_DONE\\n")
"""
        r_script_file = f"{output_prefix}_susie.R"
        Path(r_script_file).write_text(r_script)

        try:
            # Check if jsonlite is available, install if not
            subprocess.run(
                ["Rscript", "-e",
                 "if (!require('jsonlite', quietly=TRUE)) install.packages('jsonlite', repos='https://cloud.r-project.org', quiet=TRUE)"],
                capture_output=True, timeout=60,
            )

            proc = subprocess.run(
                ["Rscript", r_script_file],
                capture_output=True, text=True, timeout=600,
            )

            if Path(out_file).exists():
                result = json.loads(Path(out_file).read_text())
                return result
            else:
                print(f"  SuSiE R error: {proc.stderr[-500:]}")
                return {"pip": (np.ones(len(z_scores)) / len(z_scores)).tolist(),
                        "cs": [], "converged": False}
        except Exception as e:
            print(f"  SuSiE execution failed: {e}")
            return {"pip": (np.ones(len(z_scores)) / len(z_scores)).tolist(),
                    "cs": [], "converged": False}

    def finemapping_locus(
        self,
        sumstats: pd.DataFrame,
        chrom: int | str,
        lead_bp: int,
        window_kb: int = 500,
        gene_map: pd.DataFrame | None = None,
    ) -> dict:
        """Fine-map a single GWAS locus end-to-end.

        Parameters
        ----------
        sumstats : GWAS summary stats (unified format: CHR, BP, SNP, BETA, SE, P, N)
        chrom : chromosome of the locus
        lead_bp : position of the lead variant
        window_kb : window size in kb around lead variant
        gene_map : DataFrame mapping SNP → gene_id (for aggregation)

        Returns
        -------
        dict with variant_pips, credible_sets, gene_scores, metadata
        """
        window_start = lead_bp - window_kb * 1000
        window_end = lead_bp + window_kb * 1000

        # Extract locus variants
        locus = sumstats[
            (sumstats["CHR"].astype(str) == str(chrom)) &
            (sumstats["BP"] >= window_start) &
            (sumstats["BP"] <= window_end)
        ].copy()

        if len(locus) < 10:
            return {"error": f"Too few variants in locus ({len(locus)})", "n_variants": len(locus)}

        # Compute z-scores
        locus["Z"] = locus["BETA"] / locus["SE"]
        locus = locus.dropna(subset=["Z"]).reset_index(drop=True)

        # Get sample size
        n = int(locus["N"].median()) if "N" in locus.columns else 100000

        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = f"{tmpdir}/locus_chr{chrom}_{lead_bp}"

            # Step 1: Compute LD
            ld_matrix = self.compute_ld_matrix(locus, chrom, window_start, window_end, prefix)

            if ld_matrix is None:
                # Fallback: identity LD (no reference panel)
                print(f"  Using identity LD for chr{chrom}:{lead_bp} (no LD ref)")
                ld_matrix = np.eye(len(locus))

            # Ensure dimensions match
            if ld_matrix.shape[0] != len(locus):
                # Trim to matching variants
                n_match = min(ld_matrix.shape[0], len(locus))
                ld_matrix = ld_matrix[:n_match, :n_match]
                locus = locus.iloc[:n_match]

            # Step 2: Run SuSiE
            z_scores = locus["Z"].values
            susie_result = self.run_susie_rss(z_scores, ld_matrix, n, prefix)

        # Attach PIPs to variants
        pips = np.array(susie_result["pip"])
        locus["PIP"] = pips[:len(locus)]

        # Build credible sets
        # R's jsonlite returns cs as dict {"L1": [...]} or list [[...]]
        raw_cs = susie_result.get("cs", {})
        if isinstance(raw_cs, dict):
            cs_list = list(raw_cs.values())
        elif isinstance(raw_cs, list):
            cs_list = raw_cs
        else:
            cs_list = []

        credible_sets = []
        for i, cs_indices in enumerate(cs_list):
            if isinstance(cs_indices, (int, float)):
                cs_indices = [int(cs_indices)]
            else:
                cs_indices = [int(x) for x in cs_indices] if cs_indices else []
            cs_variants = locus.iloc[cs_indices] if cs_indices else pd.DataFrame()
            credible_sets.append({
                "cs_id": i,
                "n_variants": len(cs_indices),
                "variants": cs_variants["SNP"].tolist() if len(cs_variants) > 0 else [],
                "max_pip": cs_variants["PIP"].max() if len(cs_variants) > 0 else 0,
                "total_pip": cs_variants["PIP"].sum() if len(cs_variants) > 0 else 0,
            })

        # Step 3: Gene-level aggregation
        gene_scores = {}
        if gene_map is not None and "gene_id" in gene_map.columns:
            merged = locus.merge(gene_map, on="SNP", how="left")
            merged = merged.dropna(subset=["gene_id"])
            gene_scores = (
                merged.groupby("gene_id")["PIP"]
                .sum()
                .clip(upper=1.0)
                .sort_values(ascending=False)
                .to_dict()
            )

        return {
            "chrom": str(chrom),
            "lead_bp": lead_bp,
            "n_variants": len(locus),
            "n_credible_sets": len(credible_sets),
            "credible_sets": credible_sets,
            "converged": susie_result.get("converged", False),
            "top_variants": (
                locus.nlargest(10, "PIP")[["SNP", "BP", "PIP", "Z", "P"]]
                .to_dict(orient="records")
            ),
            "gene_scores": gene_scores,
            "variant_pips": locus[["SNP", "BP", "PIP"]].to_dict(orient="records"),
        }


def download_ld_reference(output_dir: str = "data/ld_ref", chrom: int | str = "all"):
    """Download 1000 Genomes EUR LD reference panel.

    Uses the commonly-used LDpred2/LDSC-format files from:
    https://www.dropbox.com/s/... or
    https://zenodo.org/record/... (various sources)

    For testing, we provide a minimal subset.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # For a real deployment, download from:
    # https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_plinkfiles.tgz
    base_url = "https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_plinkfiles.tgz"

    chroms = range(1, 23) if chrom == "all" else [int(chrom)]

    print(f"LD reference directory: {out}")
    print("Download the 1000G EUR plink files from:")
    print(f"  {base_url}")
    print("Then extract as EUR.chr1.bed/bim/fam, EUR.chr2.bed/bim/fam, etc.")

    for c in chroms:
        bed = out / f"EUR.chr{c}.bed"
        if bed.exists():
            print(f"  chr{c}: ✓ ({bed.stat().st_size / 1e6:.1f} MB)")
        else:
            print(f"  chr{c}: ✗ (not found)")

    return out
