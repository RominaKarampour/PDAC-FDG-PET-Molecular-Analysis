
from pathlib import Path
import pandas as pd

# STEP 1: Find our project folders
project = Path.home() / "Desktop" / "PDAC-FDG-PET-Molecular-Analysis"
data_folder = project / "data"
results_folder = project / "results"

# STEP 2: Read the selected PDAC cell-line IDs
pdac_ids = pd.read_csv(
    results_folder / "pdac_model_ids.csv"
)["ModelID"]

# STEP 3: Read the original RNA expression dataset
rna = pd.read_csv(
    data_folder / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
)

# STEP 4: Keep only default RNA profiles
rna = rna[
    rna["IsDefaultEntryForModel"]
    .astype(str)
    .str.strip()
    .str.lower()
    .isin(["yes", "true", "1"])
].copy()

# STEP 5: Keep the 52 selected PDAC cell lines
rna = rna[rna["ModelID"].isin(pdac_ids)].copy()

# Check that each cell line appears exactly once
if len(rna) != 52 or rna["ModelID"].duplicated().any():
    raise ValueError("Expected exactly 52 unique PDAC RNA profiles")

# STEP 6: Read the Hallmark Glycolysis gene set
gmt_file = data_folder / "h.all.v2026.2.Hs.symbols.gmt"

glycolysis_genes = None

with open(gmt_file) as file:
    for line in file:
        parts = line.strip().split("\t")

        if parts[0] == "HALLMARK_GLYCOLYSIS":
            glycolysis_genes = parts[2:]
            break

if glycolysis_genes is None:
    raise ValueError("HALLMARK_GLYCOLYSIS not found")

# STEP 7: Match gene symbols to RNA columns
gene_columns = {}

for column in rna.columns:
    if " (" in column:
        gene_symbol = column.split(" (")[0]
        gene_columns[gene_symbol] = column

matched_genes = [
    gene for gene in glycolysis_genes
    if gene in gene_columns
]

columns = [gene_columns[gene] for gene in matched_genes]

print("Genes in Hallmark Glycolysis:", len(glycolysis_genes))
print("Genes matched to DepMap RNA:", len(matched_genes))

# STEP 8: Identify the 50 reference cell lines
target_ids = ["ACH-000320", "ACH-000270"]

reference = rna[~rna["ModelID"].isin(target_ids)]

if len(reference) != 50:
    raise ValueError("Expected 50 reference cell lines")

# STEP 9: Calculate the reference mean and standard deviation
reference_values = reference[columns].astype(float)

gene_means = reference_values.mean(axis=0)
gene_std = reference_values.std(axis=0, ddof=0)

# Exclude genes that cannot be standardized
valid_columns = gene_std[
    gene_std.notna() & (gene_std > 0)
].index.tolist()

print("Genes with usable variation:", len(valid_columns))

# STEP 10: Calculate z-scores for all 52 cell lines
expression = rna[valid_columns].astype(float)

z_scores = (
    expression - gene_means[valid_columns]
) / gene_std[valid_columns]

# STEP 11: Average the gene z-scores
rna["GlycolysisScore"] = z_scores.mean(axis=1)

# STEP 12: Display PSN1 and HPAC scores
targets = rna[rna["ModelID"].isin(target_ids)]

print("\nHALLMARK GLYCOLYSIS SCORES")

for name, model_id in [
    ("PSN1", "ACH-000320"),
    ("HPAC", "ACH-000270")
]:
    score = targets.loc[
        targets["ModelID"] == model_id,
        "GlycolysisScore"
    ].iloc[0]

    print(f"{name}: {score:.6f}")

# STEP 13: Save scores for all 52 PDAC cell lines
rna[["ModelID", "GlycolysisScore"]].to_csv(
    results_folder / "hallmark_glycolysis_scores.csv",
    index=False
)

print("\nSaved: results/hallmark_glycolysis_scores.csv")

