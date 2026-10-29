
from pathlib import Path
import pandas as pd

# Find the folder where our project is stored
project = Path.home() / "Desktop" / "PDAC-FDG-PET-Molecular-Analysis"

# Find the RNA data file inside the data folder
rna_file = project / "data" / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"

# Read the CSV file
rna = pd.read_csv(rna_file)

# Print basic information
print("Number of rows:", rna.shape[0])
print("Number of columns:", rna.shape[1])

# Show the names of the first 10 columns
print("\nFirst 10 column names:")
print(rna.columns[:10].tolist())

# Show the first 5 cell-line IDs
print("\nFirst 5 ModelIDs:")
print(rna["ModelID"].head().tolist())

# Check whether PSN1 and HPAC are present
print("\nPSN1 found:", "ACH-000320" in rna["ModelID"].values)
print("HPAC found:", "ACH-000270" in rna["ModelID"].values)


