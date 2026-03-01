"""Parse GWAS summary statistics into a unified format."""
import pandas as pd
from pathlib import Path
from .traits import TraitMetadata

# Standard column schema
UNIFIED_COLUMNS = [
    "CHR", "BP", "SNP", "A1", "A2", "BETA", "SE", "P", "N", "FREQ"
]

# Known column mappings for different sources
COLUMN_MAPS = {
    "DIAMANTE": {"Chr": "CHR", "Pos": "BP", "rsID": "SNP", "EA": "A1", "NEA": "A2",
                 "Beta": "BETA", "SE": "SE", "Pvalue": "P", "Neff": "N", "EAF": "FREQ"},
    "CARDIoGRAMplusC4D": {"chromosome": "CHR", "base_pair_location": "BP",
                          "variant_id": "SNP", "effect_allele": "A1",
                          "other_allele": "A2", "beta": "BETA", "standard_error": "SE",
                          "p_value": "P", "effect_allele_frequency": "FREQ"},
    "PGC3": {"CHR": "CHR", "BP": "BP", "SNP": "SNP", "A1": "A1", "A2": "A2",
             "OR": "_OR", "SE": "SE", "P": "P", "Nca": "_NCA", "Nco": "_NCO", "FRQ_A_71460": "FREQ"},
    "IIBDGC": {"chromosome": "CHR", "base_pair_location": "BP",
               "variant_id": "SNP", "effect_allele": "A1",
               "other_allele": "A2", "beta": "BETA", "standard_error": "SE",
               "p_value": "P", "effect_allele_frequency": "FREQ"},
}


class SumstatParser:
    """Parse diverse GWAS formats into a unified schema."""

    def parse(self, path: Path, trait: TraitMetadata, nrows: int | None = None) -> pd.DataFrame:
        sep = "\t" if ".tsv" in path.name else None
        df = pd.read_csv(path, sep=sep, nrows=nrows, comment="#",
                         low_memory=False, na_values=["NA", ".", ""])
        col_map = COLUMN_MAPS.get(trait.source, {})
        df = df.rename(columns=col_map)

        # Handle OR -> BETA conversion (SCZ uses OR)
        if "_OR" in df.columns and "BETA" not in df.columns:
            import numpy as np
            df["BETA"] = np.log(df["_OR"].astype(float))

        # Handle split N
        if "_NCA" in df.columns and "_NCO" in df.columns and "N" not in df.columns:
            df["N"] = df["_NCA"].astype(float) + df["_NCO"].astype(float)

        # Keep unified columns that exist
        available = [c for c in UNIFIED_COLUMNS if c in df.columns]
        df = df[available].copy()

        # Basic QC
        if "P" in df.columns:
            df["P"] = pd.to_numeric(df["P"], errors="coerce")
            df = df.dropna(subset=["P"])

        df.attrs["trait_id"] = trait.trait_id
        df.attrs["source"] = trait.source
        return df
