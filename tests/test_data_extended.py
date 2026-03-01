"""Tests for validated genes and drug database."""
from gwas_loop.data.validated_genes import VALIDATED_GENES, ValidationDB
from gwas_loop.data.drug_database import get_drug_dataframe, get_drugs_for_gene, summary


def test_validated_genes_four_traits():
    assert set(VALIDATED_GENES.keys()) == {"T2D", "CAD", "SCZ", "IBD"}

def test_validated_genes_size():
    for tid, vs in VALIDATED_GENES.items():
        assert len(vs.genes) >= 30, f"{tid} has too few validated genes"
        assert vs.confidence == "gold"

def test_validated_genes_known():
    assert "TCF7L2" in VALIDATED_GENES["T2D"].genes
    assert "PCSK9" in VALIDATED_GENES["CAD"].genes
    assert "DRD2" in VALIDATED_GENES["SCZ"].genes
    assert "NOD2" in VALIDATED_GENES["IBD"].genes

def test_validation_db(tmp_path):
    db = ValidationDB(tmp_path / "val.json")
    genes = db.get("T2D")
    assert "TCF7L2" in genes
    summary = db.summary()
    assert len(summary) == 4

def test_drug_database():
    df = get_drug_dataframe()
    assert len(df) >= 30
    assert set(df.columns) == {"gene_id", "drug_name", "drug_type", "mechanism_of_action", "max_phase", "source", "indication"}

def test_drug_lookup():
    pcsk9 = get_drugs_for_gene("PCSK9")
    assert len(pcsk9) >= 3  # Evolocumab, Alirocumab, Inclisiran
    assert "Evolocumab" in pcsk9["drug_name"].values

def test_drug_summary():
    s = summary()
    assert "PCSK9" in s.index
    assert s.loc["PCSK9", "n_drugs"] >= 3
