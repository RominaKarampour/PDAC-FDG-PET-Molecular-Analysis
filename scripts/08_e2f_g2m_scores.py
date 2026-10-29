
from pathlib import Path
import pandas as pd

# 1. Find our project folders
project = Path.home() / "Desktop" / "PDAC-FDG-PET-Molecular-Analysis"
data = project / "data"
results = project / "results"

# 2. Read the 52 selected PDAC model IDs
pdac_ids = pd.read_csv(
    results / "pdac_model_ids.csv"
)["ModelID"]

# 3. Read the RNA dataset
rna = pd.read_csv(
    data / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
)

# 4. Keep default RNA profiles
rna = rna[
    rna["IsDefaultEntryForModel"]
    .astype(str)
    .str.strip()
    .str.lower()
    .isin(["yes", "true", "1"])
].copy()

# 5. Keep our 52 PDAC cell lines
rna = rna[rna["ModelID"].isin(pdac_ids)].copy()

if len(rna) != 52 or rna["ModelID"].duplicated().any():
    raise ValueError("Expected 52 unique PDAC RNA profiles")

# 6. Define the two Hallmark pathways
pathways = [
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT"
]

# 7. Read gene sets from the MSigDB GMT file
gene_sets = {}

with open(data / "h.all.v2026.2.Hs.symbols.gmt") as file:
    for line in file:
        parts = line.strip().split("\t")
        if parts[0] in pathways:
            gene_sets[parts[0]] = parts[2:]

# 8. Match gene symbols to RNA column names
gene_columns = {}

for column in rna.columns:
    if " (" in column:
        gene = column.split(" (")[0]
        gene_columns[gene] = column

# 9. Identify PSN1, HPAC, and reference lines
target_ids = ["ACH-000320", "ACH-000270"]

reference = rna[
    ~rna["ModelID"].isin(target_ids)
]

if len(reference) != 50:
    raise ValueError("Expected 50 reference cell lines")

# 10. Calculate each pathway score
for pathway in pathways:

    if pathway not in gene_sets:
        raise ValueError(f"Missing pathway: {pathway}")

    genes = gene_sets[pathway]

    matched = [
        gene for gene in genes
        if gene in gene_columns
    ]

    columns = [
        gene_columns[gene]
        for gene in matched
    ]

    # Calculate gene means and standard deviations
    # using only the 50 reference cell lines
    ref_values = reference[columns].astype(float)

    means = ref_values.mean(axis=0)
    stds = ref_values.std(axis=0, ddof=0)

    valid = stds[
        stds.notna() & (stds > 0)
    ].index.tolist()

    # Standardize all 52 cell lines
    z = (
        rna[valid].astype(float) - means[valid]
    ) / stds[valid]

    # Average standardized expression per cell line
    rna[pathway] = z.mean(axis=1)

    print("\n", pathway)
    print("Genes in set:", len(genes))
    print("Genes matched:", len(matched))
    print("Genes usable:", len(valid))

    for name, model_id in [
        ("PSN1", "ACH-000320"),
        ("HPAC", "ACH-000270")
    ]:
        score = rna.loc[
            rna["ModelID"] == model_id,
            pathway
        ].iloc[0]

        print(f"{name}: {score:.6f}")

# 11. Save both pathway scores for all 52 cell lines
rna[
    ["ModelID"] + pathways
].to_csv(
    results / "e2f_g2m_scores.csv",
    index=False
)

print("\nSaved: results/e2f_g2m_scores.csv")

