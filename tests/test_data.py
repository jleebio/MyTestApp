"""Tests for GWAS data modules."""
from gwas_loop.data.traits import TraitRegistry, TRAIT_CATALOG, TraitMetadata

def test_catalog_has_four_traits():
    assert len(TRAIT_CATALOG) >= 4
    assert set(TRAIT_CATALOG.keys()) == {"T2D", "CAD", "SCZ", "IBD"}

def test_trait_metadata_fields():
    for tid, t in TRAIT_CATALOG.items():
        assert t.sample_size > 100_000, f"{tid} N < 100k"
        assert t.population == "European"
        assert t.genome_build == "GRCh37"
        assert t.download_url.startswith("http")
        assert t.study_pmid

def test_trait_diversity():
    cats = {t.category for t in TRAIT_CATALOG.values()}
    assert len(cats) >= 4, f"Need diverse architectures, got: {cats}"

def test_registry_init(tmp_path):
    reg = TraitRegistry(tmp_path / "traits.json")
    assert set(reg.list_traits()) == {"T2D", "CAD", "SCZ", "IBD"}
    summary = reg.summary()
    assert len(summary) == 4
    assert all(s["N"] > 100_000 for s in summary)

def test_registry_add(tmp_path):
    reg = TraitRegistry(tmp_path / "traits.json")
    reg.add(TraitMetadata(
        trait_id="BMI", name="Body Mass Index", category="anthropometric",
        sample_size=700000, n_cases=0, n_controls=700000,
        population="European", study_pmid="999", first_author="Test",
        year=2023, source="GIANT", download_url="https://example.com/bmi.gz",
    ))
    assert "BMI" in reg.list_traits()
    # Reload from disk
    reg2 = TraitRegistry(tmp_path / "traits.json")
    assert "BMI" in reg2.list_traits()
