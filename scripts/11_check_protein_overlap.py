
from pathlib import Path
import pandas as pd

project = Path(__file__).resolve().parent.parent
results = project / "results"
data = project / "data"

# Load the 52 PDAC cell-line identifiers
models = pd.read_csv(results / "pdac_model_ids.csv")

# Load the proteomics measurements
protein = pd.read_csv(
    data / "harmonized_Sanger_MS_2022.csv"
)

# Identify the column containing DepMap model IDs
protein = protein.rename(
    columns={"Unnamed: 0": "ModelID"}
)

# GLUT1 and HK2 UniProt identifiers
protein_genes = {
    "GLUT1": "P11166",
    "HK2": "P52789"
}

# Find the ModelID column in the PDAC model list
if "ModelID" not in models.columns:
    raise ValueError("ModelID column not found in PDAC model list")

# Keep PDAC models with proteomics measurements
merged = models[["ModelID"]].merge(
    protein,
    on="ModelID",
    how="inner"
)

print("PDAC cell lines with proteomics:", len(merged))

# Exclude PSN1 and HPAC from the reference group
reference = merged[
    ~merged["ModelID"].isin(["ACH-000320", "ACH-000270"])
]

print("Reference PDAC lines with proteomics:", len(reference))

for name, uniprot in protein_genes.items():
    count = reference[uniprot].notna().sum()
    print(f"{name}: {count} reference lines with measurements")

for name, model_id in [
    ("PSN1", "ACH-000320"),
    ("HPAC", "ACH-000270")
]:
    rows = merged[merged["ModelID"] == model_id]
    print(f"{name} present:", len(rows) == 1)

