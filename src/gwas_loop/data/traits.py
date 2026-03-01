"""Trait metadata and GWAS summary statistics catalog."""
from dataclasses import dataclass, field, asdict
import json
from pathlib import Path


@dataclass
class TraitMetadata:
    trait_id: str
    name: str
    category: str
    sample_size: int
    n_cases: int
    n_controls: int
    population: str
    study_pmid: str
    first_author: str
    year: int
    source: str
    download_url: str
    genome_build: str = "GRCh37"
    ld_reference: str = "1000G_EUR"
    n_variants: int = 0
    notes: str = ""


# Curated catalog of high-quality public GWAS datasets
TRAIT_CATALOG = {
    "T2D": TraitMetadata(
        trait_id="T2D",
        name="Type 2 Diabetes",
        category="metabolic",
        sample_size=898130,
        n_cases=180834,
        n_controls=1159055,
        population="European",
        study_pmid="36333501",
        first_author="Mahajan",
        year=2022,
        source="DIAMANTE",
        download_url="https://diagram-consortium.org/downloads/DIAMANTE-EUR.sumstat.txt.gz",
        genome_build="GRCh37",
        notes="DIAMANTE trans-ancestry, EUR subset. Nature Genetics 2022.",
    ),
    "CAD": TraitMetadata(
        trait_id="CAD",
        name="Coronary Artery Disease",
        category="cardiovascular",
        sample_size=1165690,
        n_cases=181522,
        n_controls=984168,
        population="European",
        study_pmid="36474045",
        first_author="Aragam",
        year=2022,
        source="CARDIoGRAMplusC4D",
        download_url="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132314/GCST90132314_buildGRCh37.tsv.gz",
        genome_build="GRCh37",
        notes="CARDIoGRAMplusC4D + UK Biobank meta-analysis. Nature Genetics 2022.",
    ),
    "SCZ": TraitMetadata(
        trait_id="SCZ",
        name="Schizophrenia",
        category="psychiatric",
        sample_size=320404,
        n_cases=76755,
        n_controls=243649,
        population="European",
        study_pmid="35396580",
        first_author="Trubetskoy",
        year=2022,
        source="PGC3",
        download_url="https://figshare.com/ndownloader/files/34517828",
        genome_build="GRCh37",
        notes="PGC3 wave 3. Nature 2022. 287 loci.",
    ),
    "IBD": TraitMetadata(
        trait_id="IBD",
        name="Inflammatory Bowel Disease",
        category="autoimmune",
        sample_size=296917,
        n_cases=25042,
        n_controls=271875,
        population="European",
        study_pmid="36038634",
        first_author="Lange",
        year=2023,
        source="IIBDGC",
        download_url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90244001-GCST90245000/GCST90244763/GCST90244763_buildGRCh37.tsv.gz",
        genome_build="GRCh37",
        notes="International IBD Genetics Consortium. Nature Genetics 2023.",
    ),
}


class TraitRegistry:
    """Manage trait metadata with persistence."""

    def __init__(self, path: str | Path = "trait_metadata.json"):
        self.path = Path(path)
        self.traits: dict[str, TraitMetadata] = {}
        if self.path.exists():
            data = json.loads(self.path.read_text())
            self.traits = {k: TraitMetadata(**v) for k, v in data.items()}
        else:
            self.traits = {k: TraitMetadata(**asdict(v)) for k, v in TRAIT_CATALOG.items()}
            self._save()

    def get(self, trait_id: str) -> TraitMetadata:
        return self.traits[trait_id]

    def list_traits(self) -> list[str]:
        return list(self.traits.keys())

    def summary(self) -> list[dict]:
        return [
            {"id": t.trait_id, "name": t.name, "category": t.category,
             "N": t.sample_size, "source": t.source, "year": t.year}
            for t in self.traits.values()
        ]

    def add(self, trait: TraitMetadata):
        self.traits[trait.trait_id] = trait
        self._save()

    def _save(self):
        self.path.write_text(json.dumps(
            {k: asdict(v) for k, v in self.traits.items()}, indent=2
        ))
