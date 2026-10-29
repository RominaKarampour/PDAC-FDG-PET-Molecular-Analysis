
from pathlib import Path
import pandas as pd

# STEP 1: Find our project folder
project = Path(__file__).resolve().parent.parent

# STEP 2: Locate the Sanger proteomics file
protein_file = project / "data" / "harmonized_Sanger_MS_2022.csv"

# STEP 3: Open the file as a table
proteins = pd.read_csv(protein_file)

# STEP 4: Count the rows and columns
print("Number of rows:", proteins.shape[0])
print("Number of columns:", proteins.shape[1])

# STEP 5: Display the first 15 column names
print("\nFirst 15 column names:")
print(proteins.columns[:15].tolist())

# STEP 6: Display the first 3 rows of the first 6 columns
print("\nPreview of the data:")
print(proteins.iloc[:3, :6].to_string(index=False))


