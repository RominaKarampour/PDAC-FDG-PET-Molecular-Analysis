
from pathlib import Path
import pandas as pd
from scipy.stats import linregress

# STEP 1: Locate the project folders
project = Path(__file__).resolve().parent.parent
data = project / "data"
results = project / "results"

# STEP 2: Load E2F and G2/M scores
scores = pd.read_csv(
    results / "e2f_g2m_scores.csv"
)

# STEP 3: Load RNA expression
rna = pd.read_csv(
    data / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
)

# STEP 4: Keep default RNA profiles
rna = rna[
    rna["IsDefaultEntryForModel"]
    .astype(str)
    .str.strip()
    .str.lower()
    .isin(["yes", "true", "1"])
].copy()

# STEP 5: Find GLUT1 and HK1 gene columns
genes = {
    "GLUT1": "SLC2A1",
    "HK2": "HK2"
}

gene_columns = {}

for label, gene in genes.items():
    matches = [
        col for col in rna.columns
        if col.startswith(gene + " (")
    ]

    if len(matches) != 1:
        raise ValueError(
            f"Expected one RNA column for {gene}"
        )

    gene_columns[label] = matches[0]

# STEP 6: Combine pathway scores and RNA data
selected = rna[
    ["ModelID"] + list(gene_columns.values())
].copy()

selected = selected.rename(
    columns={value: key for key, value in gene_columns.items()}
)

df = scores.merge(
    selected,
    on="ModelID",
    validate="one_to_one"
)

# STEP 7: Define the two cell lines
targets = {
    "PSN1": "ACH-000320",
    "HPAC": "ACH-000270"
}

# STEP 8: Select 50 reference cell lines
reference = df[
    ~df["ModelID"].isin(targets.values())
]

if len(reference) != 50:
    raise ValueError("Expected 50 reference cell lines")

# STEP 9: Define the pathways
pathways = {
    "E2F": "HALLMARK_E2F_TARGETS",
    "G2M": "HALLMARK_G2M_CHECKPOINT"
}

# STEP 10: Calculate residuals
output = []

for pathway_name, pathway_col in pathways.items():

    for gene in ["GLUT1", "HK2"]:

        # Fit a regression using the 50 reference lines
        fit = linregress(
            reference[pathway_col],
            reference[gene]
        )

        print(
            f"{pathway_name} vs {gene}: "
            f"reference r = {fit.rvalue:.3f}, "
            f"P = {fit.pvalue:.4g}"
        )

        # Calculate predicted expression and residuals
        for cell_name, model_id in targets.items():

            row = df[df["ModelID"] == model_id]

            if len(row) != 1:
                raise ValueError(
                    f"Expected one row for {cell_name}"
                )

            pathway_score = row[pathway_col].iloc[0]
            observed = row[gene].iloc[0]

            predicted = (
                fit.intercept
                + fit.slope * pathway_score
            )

            residual = observed - predicted

            output.append({
                "Pathway": pathway_name,
                "Gene": gene,
                "CellLine": cell_name,
                "ObservedRNA": observed,
                "PredictedRNA": predicted,
                "Residual": residual
            })

# STEP 11: Save the results
output_df = pd.DataFrame(output)

print("\nADJUSTED RNA RESIDUALS")
print(output_df.to_string(index=False))

output_df.to_csv(
    results / "e2f_g2m_adjusted_residuals.csv",
    index=False
)

print("\nSaved: results/e2f_g2m_adjusted_residuals.csv")

