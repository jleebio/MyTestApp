"""Tests for manuscript generation."""
from gwas_loop.manuscript.methods_section import get_methods_section, save_methods_section


def test_methods_section_content():
    text = get_methods_section()
    assert "DIAMANTE" in text
    assert "Psychiatric Genomics Consortium" in text
    assert "Disagreement" in text
    assert "Clinical Trial" in text.lower() or "clinical trial" in text.lower()
    assert len(text) > 2000


def test_save_methods(tmp_path):
    path = save_methods_section(str(tmp_path))
    assert path.exists()
    assert path.stat().st_size > 2000
