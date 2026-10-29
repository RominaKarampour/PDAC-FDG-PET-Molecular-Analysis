
from pathlib import Path
import pandas as pd

# STEP 1: Locate the project
project = Path(__file__).resolve().parent.parent

# STEP 2: Open the proteomics dataset
proteins = pd.read_csv(
    project / "data" / "harmonized_Sanger_MS_2022.csv"
)

# STEP 3: Rename the cell-line identifier column
proteins = proteins.rename(
    columns={"Unnamed: 0": "ModelID"}
)

# STEP 4: Define the two cell lines
cell_lines = {
    "PSN1": "ACH-000320",
    "HPAC": "ACH-000270"
}

# STEP 5: Define the five available proteins
protein_ids = {
    "GLUT1": "P11166",
    "HK1": "P19367",
    "HK2": "P52789",
    "LDHA": "P00338",
    "MCT4": "O15427"
}

# STEP 6: Collect the measurements
results = []

for cell_name, model_id in cell_lines.items():

    cell_data = proteins[
        proteins["ModelID"] == model_id
    ]

    if len(cell_data) != 1:
        raise ValueError(
            f"Expected one row for {cell_name}, found {len(cell_data)}"
        )

    for protein_name, uniprot_id in protein_ids.items():

        value = cell_data.iloc[0][uniprot_id]

        results.append({
            "CellLine": cell_name,
            "ModelID": model_id,
            "Protein": protein_name,
            "UniProtID": uniprot_id,
            "ProteinAbundance": value
        })

# STEP 7: Create a table
results_df = pd.DataFrame(results)

# STEP 8: Display the measurements
print("\nPROTEIN ABUNDANCE RESULTS")
print(results_df.to_string(index=False))

# STEP 9: Save the table
results_df.to_csv(
    project / "results" / "psn1_hpac_protein_abundance.csv",
    index=False
)

print("\nSaved: results/psn1_hpac_protein_abundance.csv")


