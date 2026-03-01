"""MODULE 6: Drug Mapping.

Map causal genes to approved drugs, clinical candidates, and tool compounds.
"""
import pandas as pd
from dataclasses import dataclass


@dataclass
class DrugGeneLink:
    gene_id: str
    drug_name: str
    drug_type: str          # small_molecule, antibody, antisense, etc.
    mechanism_of_action: str  # inhibitor, agonist, antagonist, modulator
    max_phase: int          # 0=preclinical, 1-3=trials, 4=approved
    source: str             # chembl, drugbank, opentargets
    indication: str


class DrugMapper:
    """Map genes to known drugs and compounds.

    Sources: ChEMBL, DrugBank, Open Targets, DGIdb.
    When no direct drug exists, performs pathway expansion:
    finds druggable genes in same pathway.
    """

    def __init__(self):
        self.drug_db: pd.DataFrame = pd.DataFrame()

    def load_drug_database(self, path_or_df):
        """Load drug-gene interaction database.

        Expected columns: gene_id, drug_name, drug_type, mechanism_of_action,
        max_phase, source, indication
        """
        if isinstance(path_or_df, pd.DataFrame):
            self.drug_db = path_or_df
        else:
            self.drug_db = pd.read_csv(path_or_df, sep="\t")

    def map_direct(self, gene_ids: list[str]) -> list[DrugGeneLink]:
        """Find drugs directly targeting these genes."""
        if self.drug_db.empty:
            return []
        hits = self.drug_db[self.drug_db["gene_id"].isin(gene_ids)]
        return [DrugGeneLink(**row.to_dict()) for _, row in hits.iterrows()]

    def map_pathway_expansion(self, gene_ids: list[str],
                               pathway_db: pd.DataFrame) -> list[DrugGeneLink]:
        """Find druggable genes in the same pathway as target genes.

        pathway_db: gene_id, pathway_id, pathway_name
        """
        if self.drug_db.empty or pathway_db.empty:
            return []

        # Find pathways containing our genes
        target_pathways = pathway_db[pathway_db["gene_id"].isin(gene_ids)]["pathway_id"].unique()
        # Find all genes in those pathways
        expanded_genes = pathway_db[pathway_db["pathway_id"].isin(target_pathways)]["gene_id"].unique()
        # Find drugs for expanded gene set (excluding original targets)
        new_genes = set(expanded_genes) - set(gene_ids)
        return self.map_direct(list(new_genes))

    def summary(self, links: list[DrugGeneLink]) -> pd.DataFrame:
        if not links:
            return pd.DataFrame()
        return pd.DataFrame([
            {"gene": l.gene_id, "drug": l.drug_name, "type": l.drug_type,
             "moa": l.mechanism_of_action, "phase": l.max_phase, "indication": l.indication}
            for l in links
        ]).sort_values("phase", ascending=False)
