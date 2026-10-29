
from pathlib import Path
import pandas as pd

# STEP 1: Find the proteomics file
project = Path.home() / "Desktop" / "PDAC-FDG-PET-Molecular-Analysis"
protein_file = project / "data" / "harmonized_Sanger_MS_2022.csv"

# STEP 2: Open the proteomics table
proteins = pd.read_csv(protein_file)

# STEP 3: Give the first column a meaningful name
proteins = proteins.rename(columns={"Unnamed: 0": "ModelID"})

# STEP 4: Define our two cell lines
cell_lines = {
    "PSN1": "ACH-000320",
    "HPAC": "ACH-000270"
}

# STEP 5: Define proteins and their UniProt identifiers
protein_ids = {
    "GLUT1": "P11166",
    "GLUT3": "P11169",
    "HK1": "P19367",
    "HK2": "P52789",
    "LDHA": "P00338",
    "MCT4": "O15427"
}

# STEP 6: Check whether each cell line exists
print("CELL-LINE CHECK")

for name, model_id in cell_lines.items():
    found = model_id in proteins["ModelID"].values
    print(name, model_id, "found:", found)

# STEP 7: Check whether each protein column exists
print("\nPROTEIN CHECK")

for name, uniprot_id in protein_ids.items():
    found = uniprot_id in proteins.columns
    print(name, uniprot_id, "found:", found)

