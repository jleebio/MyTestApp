"""Validated causal gene sets for benchmarking.

Sources:
- Open Targets Genetics gold standard (L2G training set)
- OMIM confirmed disease genes
- ClinVar pathogenic genes
- Drug target genes (approved drugs = validated biology)

These serve as ground truth for evaluating gene prioritization methods.
"""
import json
import pandas as pd
from pathlib import Path
from dataclasses import dataclass


@dataclass
class ValidationSet:
    trait_id: str
    source: str
    genes: set[str]
    confidence: str  # "gold", "silver", "bronze"
    notes: str


# Curated gold-standard causal genes (from literature + Open Targets + OMIM)
VALIDATED_GENES = {
    "T2D": ValidationSet(
        trait_id="T2D",
        source="OMIM + Open Targets L2G + literature",
        confidence="gold",
        notes="Well-established T2D causal genes from GWAS follow-up studies",
        genes={
            "TCF7L2", "SLC30A8", "KCNJ11", "ABCC8", "PPARG", "HNF1A", "HNF4A",
            "GCK", "GLIS3", "CDKAL1", "CDKN2A", "CDKN2B", "IGF2BP2", "FTO",
            "HHEX", "IDE", "TCF2", "WFS1", "JAZF1", "KCNQ1", "MTNR1B",
            "PROX1", "DGKB", "ADCY5", "GIPR", "ARAP1", "HNF1B", "PAX4",
            "INS", "IRS1", "IRS2", "PCSK1", "PDX1", "NEUROD1", "KLF11",
            "CEL", "PAM", "ANK1", "ZMIZ1", "TLE4", "GRB14", "THADA",
        },
    ),
    "CAD": ValidationSet(
        trait_id="CAD",
        source="OMIM + CARDIoGRAMplusC4D + literature",
        confidence="gold",
        notes="Established CAD genes from GWAS, Mendelian genetics, and drug targets",
        genes={
            "PCSK9", "LDLR", "APOB", "APOE", "APOA1", "CETP", "LPA", "HMGCR",
            "NPC1L1", "SORT1", "PHACTR1", "ADAMTS7", "COL4A1", "COL4A2",
            "CXCL12", "FLT1", "IL6R", "GUCY1A1", "EDNRA", "SH2B3",
            "CDKN2A", "CDKN2B", "TCF21", "ZC3HC1", "ABO", "FN1",
            "MIA3", "SMAD3", "LIPA", "PLG", "LPL", "ANGPTL4", "ANGPTL3",
            "ASGR1", "TRIB1", "PPAP2B", "REST", "KCNE2", "SLC22A3",
        },
    ),
    "SCZ": ValidationSet(
        trait_id="SCZ",
        source="PGC3 + OMIM + literature",
        confidence="gold",
        notes="Established SCZ genes from GWAS fine-mapping and rare variant studies",
        genes={
            "GRIN2A", "SP4", "TRIO", "CACNA1C", "DRD2", "GRM3", "SETD1A",
            "CUL1", "XKR6", "FOXP1", "TCF4", "NRGN", "ZNF804A", "MIR137",
            "CACNA1I", "RIMS1", "SNAP91", "CLCN3", "GRIA1", "GRIN2B",
            "SRR", "FES", "FURIN", "TSNARE1", "PLCL1", "MAPK3",
            "CSMD1", "MHC", "NCAN", "DPYD", "GNL3", "GIGYF2",
            "SLC39A8", "VRK2", "AS3MT", "IGSF9B", "CNKSR2",
        },
    ),
    "IBD": ValidationSet(
        trait_id="IBD",
        source="IIBDGC + OMIM + literature",
        confidence="gold",
        notes="Established IBD genes from GWAS, Crohn's/UC follow-up, and rare variants",
        genes={
            "NOD2", "IL23R", "ATG16L1", "IRGM", "IL10", "IL10RA", "IL10RB",
            "CARD9", "PTPN2", "TNFSF15", "JAK2", "STAT3", "TYK2",
            "MST1", "NKX2-3", "SMAD3", "ERAP2", "FUT2", "LRRK2",
            "PTPN22", "TNFAIP3", "REL", "ICOSLG", "ORMDL3", "GSDMB",
            "CCR6", "IL2RA", "IL12B", "IL18RAP", "IFIH1", "SLC22A5",
            "PTGER4", "MUC19", "HNF4A", "CDH1", "ITLN1", "RIPK2",
            "XIAP", "GP2", "ADAM17", "RTEL1", "SBNO2",
        },
    ),
}


class ValidationDB:
    """Manage validated gene sets for benchmarking."""

    def __init__(self, path: str | Path = "validated_genes.json"):
        self.path = Path(path)
        self.sets: dict[str, ValidationSet] = {}
        if self.path.exists():
            data = json.loads(self.path.read_text())
            self.sets = {
                k: ValidationSet(genes=set(v["genes"]), **{kk: vv for kk, vv in v.items() if kk != "genes"})
                for k, v in data.items()
            }
        else:
            self.sets = {k: ValidationSet(**{**v.__dict__}) for k, v in VALIDATED_GENES.items()}
            self._save()

    def get(self, trait_id: str) -> set[str]:
        return self.sets[trait_id].genes

    def summary(self) -> pd.DataFrame:
        return pd.DataFrame([
            {"trait": s.trait_id, "n_genes": len(s.genes), "confidence": s.confidence, "source": s.source}
            for s in self.sets.values()
        ])

    def _save(self):
        data = {}
        for k, v in self.sets.items():
            d = v.__dict__.copy()
            d["genes"] = sorted(list(d["genes"]))
            data[k] = d
        self.path.write_text(json.dumps(data, indent=2))
