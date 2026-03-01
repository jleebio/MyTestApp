"""Tests for real SuSiE fine-mapping and pipeline integration."""
import numpy as np
import pandas as pd
import pytest
from gwas_loop.methods.susie_real import SuSiERunner


def test_susie_runner_init():
    runner = SuSiERunner()
    assert runner.max_causal == 10
    assert runner.coverage == 0.95


def test_susie_check_r():
    runner = SuSiERunner()
    assert runner._check_r() is True


def test_susie_finemapping_too_few_variants():
    runner = SuSiERunner()
    sumstats = pd.DataFrame({
        "CHR": [1, 1, 1],
        "BP": [100000, 100100, 100200],
        "SNP": ["rs1", "rs2", "rs3"],
        "BETA": [0.1, 0.05, -0.02],
        "SE": [0.01, 0.01, 0.01],
        "P": [1e-10, 1e-5, 0.5],
        "N": [100000, 100000, 100000],
    })
    result = runner.finemapping_locus(sumstats, chrom=1, lead_bp=100100)
    assert "error" in result


def test_susie_finemapping_identity_ld():
    """Test SuSiE with identity LD (no LD ref panel — uses fallback)."""
    runner = SuSiERunner(ld_ref_dir="/nonexistent")
    np.random.seed(42)
    n_vars = 50
    betas = np.random.normal(0, 0.01, n_vars)
    betas[25] = 0.15  # strong causal signal

    sumstats = pd.DataFrame({
        "CHR": [1] * n_vars,
        "BP": np.arange(1000000, 1000000 + n_vars * 100, 100),
        "SNP": [f"rs{i}" for i in range(n_vars)],
        "BETA": betas,
        "SE": [0.01] * n_vars,
        "P": [0.5] * n_vars,
        "N": [100000] * n_vars,
    })
    from scipy import stats
    z = sumstats["BETA"] / sumstats["SE"]
    sumstats["P"] = 2 * stats.norm.sf(np.abs(z))

    result = runner.finemapping_locus(sumstats, chrom=1, lead_bp=1002500, window_kb=10)

    assert result["n_variants"] == n_vars
    assert len(result["top_variants"]) > 0
    # The strong signal (rs25) should have high PIP
    top_snps = [v["SNP"] for v in result["top_variants"][:5]]
    assert "rs25" in top_snps, f"Causal variant rs25 not in top 5: {top_snps}"
    # Check credible sets exist
    assert result["n_credible_sets"] >= 1 or result["converged"]


def test_susie_pip_sum_reasonable():
    """PIPs should be between 0 and 1."""
    runner = SuSiERunner(ld_ref_dir="/nonexistent")
    np.random.seed(123)
    n_vars = 30
    betas = np.random.normal(0, 0.01, n_vars)
    betas[10] = 0.2

    sumstats = pd.DataFrame({
        "CHR": [2] * n_vars,
        "BP": np.arange(5000000, 5000000 + n_vars * 200, 200),
        "SNP": [f"rs{100+i}" for i in range(n_vars)],
        "BETA": betas,
        "SE": [0.01] * n_vars,
        "P": [0.5] * n_vars,
        "N": [200000] * n_vars,
    })
    from scipy import stats
    z = sumstats["BETA"] / sumstats["SE"]
    sumstats["P"] = 2 * stats.norm.sf(np.abs(z))

    result = runner.finemapping_locus(sumstats, chrom=2, lead_bp=5003000, window_kb=5)
    pips = [v["PIP"] for v in result["variant_pips"]]
    assert all(0 <= p <= 1 for p in pips), f"Invalid PIPs: {pips}"


def test_pipeline_init():
    from gwas_loop.pipeline.locus_pipeline import LocusPipeline
    pipeline = LocusPipeline(ld_ref_dir="/nonexistent")
    assert pipeline.susie is not None
    assert pipeline.scorer is not None
    assert pipeline.trial_predictor is not None


def test_full_pipeline_identify_loci():
    from gwas_loop.pipeline.locus_pipeline import FullPipeline
    pipeline = FullPipeline(ld_ref_dir="/nonexistent")
    sumstats = pd.DataFrame({
        "CHR": [1, 1, 2, 2, 3],
        "BP": [1000000, 1000500, 5000000, 50000000, 8000000],
        "SNP": ["rs1", "rs2", "rs3", "rs4", "rs5"],
        "P": [1e-10, 1e-8, 1e-12, 1e-9, 0.5],
    })
    loci = pipeline.identify_loci(sumstats)
    # Should find 3 loci (chr1 clumps together, chr2 has 2, chr3 ns)
    assert len(loci) >= 2
    assert loci[0]["lead_p"] < loci[-1]["lead_p"]  # sorted by p
