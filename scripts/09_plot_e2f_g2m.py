
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# STEP 1: Locate our files
project = Path(__file__).resolve().parent.parent
data = project / "data"
results = project / "results"

# STEP 2: Load the pathway scores
scores = pd.read_csv(results / "e2f_g2m_scores.csv")

# STEP 3: Load RNA data
rna = pd.read_csv(
    data / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
)

rna = rna[
    rna["IsDefaultEntryForModel"]
    .astype(str).str.strip().str.lower()
    .isin(["yes", "true", "1"])
].copy()

# STEP 4: Select the two RNA genes
rna_genes = {
    "SLC2A1 RNA": "SLC2A1 (6513)",
    "HK1 RNA": "HK1 (3098)"
}

for label, column in rna_genes.items():
    if column not in rna.columns:
        raise ValueError(f"Missing RNA column: {column}")

rna = rna[
    ["ModelID"] + list(rna_genes.values())
].rename(
    columns={value: key for key, value in rna_genes.items()}
)

# STEP 5: Load protein data
proteins = pd.read_csv(
    data / "harmonized_Sanger_MS_2022.csv"
).rename(columns={"Unnamed: 0": "ModelID"})

protein_columns = {
    "GLUT1 protein": "P11166",
    "HK1 protein": "P19367"
}

proteins = proteins[
    ["ModelID"] + list(protein_columns.values())
].rename(
    columns={value: key for key, value in protein_columns.items()}
)

# STEP 6: Combine measurements using cell-line IDs
df = scores.merge(rna, on="ModelID", validate="one_to_one")
df = df.merge(proteins, on="ModelID", how="left", validate="one_to_one")

print("PDAC cell lines:", len(df))

# STEP 7: Identify our two highlighted cell lines
targets = {
    "ACH-000320": "PSN1",
    "ACH-000270": "HPAC"
}

reference = df[~df["ModelID"].isin(targets)]

# STEP 8: Function to create one two-panel figure
def make_figure(pathway, measurements, filename, ylabel):

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    for ax, measurement in zip(axes, measurements):

        # Use only reference lines with complete measurements
        clean = reference[
            [pathway, measurement]
        ].dropna()

        x = clean[pathway].astype(float)
        y = clean[measurement].astype(float)

        # Draw reference cell lines
        ax.scatter(
            x, y,
            color="lightgray",
            edgecolor="gray",
            s=35,
            label="Reference PDAC lines"
        )

        # Fit regression using reference lines only
        if len(clean) >= 3 and x.nunique() > 1:
            fit = linregress(x, y)

            x_line = np.linspace(x.min(), x.max(), 100)
            y_line = fit.intercept + fit.slope * x_line

            ax.plot(x_line, y_line, color="black", linewidth=1.5)

            ax.text(
                0.04, 0.96,
                f"Reference n = {len(clean)}\n"
                f"r = {fit.rvalue:.3f}\n"
                f"P = {fit.pvalue:.3g}",
                transform=ax.transAxes,
                va="top",
                fontsize=9
            )

        # Highlight PSN1 and HPAC
        for model_id, name in targets.items():
            row = df[df["ModelID"] == model_id]

            if len(row) == 1:
                x_value = row[pathway].iloc[0]
                y_value = row[measurement].iloc[0]

                if pd.notna(x_value) and pd.notna(y_value):
                    color = "crimson" if name == "PSN1" else "royalblue"

                    ax.scatter(
                        x_value, y_value,
                        color=color,
                        edgecolor="black",
                        s=100,
                        zorder=5,
                        label=name
                    )

                    ax.annotate(
                        name,
                        (x_value, y_value),
                        xytext=(6, 6),
                        textcoords="offset points",
                        fontsize=9
                    )

        ax.set_xlabel(pathway.replace("HALLMARK_", "").replace("_", " "))
        ax.set_ylabel(ylabel)
        ax.set_title(measurement)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center",
               ncol=3, bbox_to_anchor=(0.5, 1.06))

    fig.tight_layout()

    fig.savefig(results / f"{filename}.pdf", bbox_inches="tight")
    fig.savefig(results / f"{filename}.png",
                dpi=300, bbox_inches="tight")

    plt.close(fig)

    print(f"Saved: {filename}.pdf and {filename}.png")

# STEP 9: Create the four figures
make_figure(
    "HALLMARK_E2F_TARGETS",
    ["SLC2A1 RNA", "HK1 RNA"],
    "e2f_vs_rna",
    "RNA expression, log2(TPM + 1)"
)

make_figure(
    "HALLMARK_G2M_CHECKPOINT",
    ["SLC2A1 RNA", "HK1 RNA"],
    "g2m_vs_rna",
    "RNA expression, log2(TPM + 1)"
)

make_figure(
    "HALLMARK_E2F_TARGETS",
    ["GLUT1 protein", "HK1 protein"],
    "e2f_vs_protein",
    "Sanger relative protein abundance"
)

make_figure(
    "HALLMARK_G2M_CHECKPOINT",
    ["GLUT1 protein", "HK1 protein"],
    "g2m_vs_protein",
    "Sanger relative protein abundance"
)

print("\nAll four figures created.")

