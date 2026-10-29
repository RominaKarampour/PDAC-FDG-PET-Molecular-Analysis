
from pathlib import Path
import pandas as pd
from scipy.stats import linregress

# Locate our project files
project = Path(__file__).resolve().parent.parent
data = project / "data"
results = project / "results"

# Read the E2F and G2/M pathway scores
scores = pd.read_csv(
    results / "e2f_g2m_scores.csv"
)

# Read the original Sanger proteomics data
protein = pd.read_csv(
    data / "harmonized_Sanger_MS_2022.csv"
)

protein = protein.rename(
    columns={"Unnamed: 0": "ModelID"}
)

# Keep only the two proteins of interest
protein = protein[
    ["ModelID", "P11166", "P52789"]
].rename(
    columns={
        "P11166": "GLUT1",
        "P52789": "HK2"
    }
)

# Combine RNA-derived pathway scores with protein data
df = scores.merge(
    protein,
    on="ModelID",
    how="inner",
    validate="one_to_one"
)

# Identify PSN1 and HPAC
targets = {
    "PSN1": "ACH-000320",
    "HPAC": "ACH-000270"
}

# Exclude PSN1 and HPAC from the regression training data
reference = df[
    ~df["ModelID"].isin(targets.values())
].copy()

# Define the two pathway scores
pathways = {
    "E2F": "HALLMARK_E2F_TARGETS",
    "G2M": "HALLMARK_G2M_CHECKPOINT"
}

output = []

# Fit one regression for each pathway/protein combination
for pathway_name, pathway_column in pathways.items():

    for protein_name in ["GLUT1", "HK2"]:

        # Remove missing measurements for this specific protein
        usable = reference[
            [pathway_column, protein_name]
        ].dropna()

        # Fit protein abundance against the RNA pathway score
        fit = linregress(
            usable[pathway_column],
            usable[protein_name]
        )

        print(
            f"{pathway_name} vs {protein_name}: "
            f"n={len(usable)}, "
            f"r={fit.rvalue:.3f}, "
            f"P={fit.pvalue:.4g}"
        )

        # Predict protein abundance for PSN1 and HPAC
        for cell_name, model_id in targets.items():

            row = df[df["ModelID"] == model_id]

            if len(row) != 1:
                raise ValueError(
                    f"Expected one row for {cell_name}"
                )

            pathway_score = row[pathway_column].iloc[0]
            observed = row[protein_name].iloc[0]

            predicted = (
                fit.intercept
                + fit.slope * pathway_score
            )

            residual = observed - predicted

            output.append({
                "Pathway": pathway_name,
                "Protein": protein_name,
                "CellLine": cell_name,
                "ReferenceN": len(usable),
                "ObservedProtein": observed,
                "PredictedProtein": predicted,
                "Residual": residual
            })

# Save all eight residuals
output_df = pd.DataFrame(output)

print("\nADJUSTED PROTEIN RESIDUALS")
print(output_df.to_string(index=False))

output_df.to_csv(
    results / "e2f_g2m_adjusted_protein_residuals.csv",
    index=False
)

print("\nSaved: results/e2f_g2m_adjusted_protein_residuals.csv")

