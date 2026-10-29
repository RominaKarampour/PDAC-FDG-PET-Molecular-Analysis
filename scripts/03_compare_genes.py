
from pathlib import Path
import pandas as pd

# STEP 1: Find our project folder
project = Path(__file__).resolve().parent.parent
data_folder = project / "data"
results_folder = project / "results"

# STEP 2: Open the RNA expression dataset
rna = pd.read_csv(
    data_folder / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
)

# STEP 3: Keep only default RNA profiles
rna = rna[
    rna["IsDefaultEntryForModel"]
    .astype(str)
    .str.strip()
    .str.lower()
    .isin(["yes", "true", "1"])
].copy()

# STEP 4: Define the cell lines we want to study
cell_lines = {
    "PSN1": "ACH-000320",
    "HPAC": "ACH-000270"
}

# STEP 5: Define the genes we want to study
genes = ["SLC2A1", "SLC2A3", "HK1", "HK2"]

# STEP 6: Match gene symbols to column names
gene_columns = {}

for gene in genes:
    matches = [
        column for column in rna.columns
        if column.startswith(gene + " (")
    ]

    if len(matches) != 1:
        raise ValueError(
            f"Expected one column for {gene}, found {len(matches)}"
        )

    gene_columns[gene] = matches[0]

# STEP 7: Extract RNA expression values
results = []

for cell_name, model_id in cell_lines.items():
    cell_data = rna[rna["ModelID"] == model_id]

    if len(cell_data) != 1:
        raise ValueError(
            f"Expected one default RNA profile for {cell_name}, "
            f"found {len(cell_data)}"
        )

    for gene in genes:
        value = cell_data.iloc[0][gene_columns[gene]]

        results.append({
            "CellLine": cell_name,
            "ModelID": model_id,
            "Gene": gene,
            "RNA_log2_TPM_plus_1": value
        })

# STEP 8: Create a results table
results_df = pd.DataFrame(results)

# STEP 9: Print the table
print("\nRNA expression results:")
print(results_df.to_string(index=False))

# STEP 10: Save the results as a CSV file
results_df.to_csv(
    results_folder / "psn1_hpac_rna_expression.csv",
    index=False
)

print("\nSaved: results/psn1_hpac_rna_expression.csv")

