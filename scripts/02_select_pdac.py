
from pathlib import Path
import pandas as pd

# 1. Locate our project and data folder
project = Path(__file__).resolve().parent.parent
data_folder = project / "data"

# 2. Open the cell-line information file
models = pd.read_csv(data_folder / "Model.csv")

# 3. Show the number of models in the database
print("Total model records:", len(models))

# 4. Select pancreatic adenocarcinoma cell lines
pdac = models[
    (models["OncotreeCode"] == "PAAD") &
    (models["ModelType"] == "Cell Line")
].copy()

# 5. Count the selected models
print("PDAC models in metadata:", len(pdac))

# 6. Open the RNA dataset
rna = pd.read_csv(
    data_folder / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
)

# 7. Keep the default RNA profile for each model
rna = rna[
    rna["IsDefaultEntryForModel"]
    .astype(str)
    .str.strip()
    .str.lower()
    .isin(["yes", "true", "1"])
].copy()

# 8. Match PDAC model information to RNA measurements
pdac_rna = pdac.merge(
    rna,
    on="ModelID",
    how="inner",
    validate="one_to_one"
)

# 9. Display the final number of PDAC cell lines
print("PDAC cell lines with default RNA profiles:", len(pdac_rna))

# 10. Check that PSN1 and HPAC are included
print("PSN1 included:", "ACH-000320" in pdac_rna["ModelID"].values)
print("HPAC included:", "ACH-000270" in pdac_rna["ModelID"].values)

# 11. Save the selected model IDs
results_folder = project / "results"
results_folder.mkdir(exist_ok=True)

pdac_rna[["ModelID"]].to_csv(
    results_folder / "pdac_model_ids.csv",
    index=False
)

print("Saved: results/pdac_model_ids.csv")

