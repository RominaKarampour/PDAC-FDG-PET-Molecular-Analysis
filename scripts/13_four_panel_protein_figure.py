
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------
# 1. Locate our project and load the results
# --------------------------------------------------

project = Path.home() / "Desktop" / "PDAC-FDG-PET-Molecular-Analysis"
results = project / "results"

rna = pd.read_csv(
    results / "psn1_hpac_rna_expression.csv"
)

protein = pd.read_csv(
    results / "psn1_hpac_protein_abundance.csv"
)

residuals = pd.read_csv(
    results / "e2f_g2m_adjusted_protein_residuals.csv"
)

# --------------------------------------------------
# 2. Choose colors
# --------------------------------------------------

RED = "#C84E52"
BLUE = "#4D76B5"

# --------------------------------------------------
# 3. Helper: retrieve exactly one measurement
# --------------------------------------------------

def get_value(table, filters, column):

    selected = table.copy()

    for name, value in filters.items():
        selected = selected[selected[name] == value]

    if len(selected) != 1:
        raise ValueError(
            f"Expected one measurement for {filters}, "
            f"but found {len(selected)}"
        )

    return float(selected[column].iloc[0])

# --------------------------------------------------
# 4. Helper: draw grouped bars
# --------------------------------------------------

def draw_bars(ax, labels, psn1, hpac, title, ylabel):

    x = np.arange(len(labels))
    width = 0.34

    ax.bar(
        x - width / 2,
        psn1,
        width,
        color=RED,
        label="QM-PSN1"
    )

    ax.bar(
        x + width / 2,
        hpac,
        width,
        color=BLUE,
        label="C-HPAC"
    )

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=10)

    ax.axhline(0, color="black", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# --------------------------------------------------
# 5. Create four panels
# --------------------------------------------------

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

# --------------------------------------------------
# PANEL A: RNA
# --------------------------------------------------

rna_genes = ["SLC2A1", "SLC2A3", "HK1", "HK2"]

rna_psn1 = [
    get_value(
        rna,
        {"CellLine": "PSN1", "Gene": gene},
        "RNA_log2_TPM_plus_1"
    )
    for gene in rna_genes
]

rna_hpac = [
    get_value(
        rna,
        {"CellLine": "HPAC", "Gene": gene},
        "RNA_log2_TPM_plus_1"
    )
    for gene in rna_genes
]

draw_bars(
    axes[0, 0],
    rna_genes,
    rna_psn1,
    rna_hpac,
    "Transcriptomics",
    "RNA expression [log2(TPM + 1)]"
)

# --------------------------------------------------
# PANEL B: PROTEOMICS
# --------------------------------------------------

proteins = ["GLUT1", "HK1", "HK2", "LDHA", "MCT4"]

protein_psn1 = [
    get_value(
        protein,
        {"CellLine": "PSN1", "Protein": name},
        "ProteinAbundance"
    )
    for name in proteins
]

protein_hpac = [
    get_value(
        protein,
        {"CellLine": "HPAC", "Protein": name},
        "ProteinAbundance"
    )
    for name in proteins
]

draw_bars(
    axes[0, 1],
    proteins,
    protein_psn1,
    protein_hpac,
    "Proteomics",
    "Relative protein abundance"
)

# --------------------------------------------------
# PANELS C AND D: PROTEIN RESIDUALS
# --------------------------------------------------

for ax, pathway, title in [
    (axes[1, 0], "E2F", "E2F-adjusted proteomics"),
    (axes[1, 1], "G2M", "G2/M-adjusted proteomics")
]:

    names = ["GLUT1", "HK2"]

    psn1_values = [
        get_value(
            residuals,
            {
                "Pathway": pathway,
                "Protein": name,
                "CellLine": "PSN1"
            },
            "Residual"
        )
        for name in names
    ]

    hpac_values = [
        get_value(
            residuals,
            {
                "Pathway": pathway,
                "Protein": name,
                "CellLine": "HPAC"
            },
            "Residual"
        )
        for name in names
    ]

    draw_bars(
        ax,
        names,
        psn1_values,
        hpac_values,
        title,
        "Protein abundance residual"
    )

    # Use identical scales for panels C and D
    ax.set_ylim(-3.6, 3.1)

# --------------------------------------------------
# 6. Italicize RNA gene symbols only
# --------------------------------------------------

for label in axes[0, 0].get_xticklabels():
    label.set_fontstyle("italic")

# --------------------------------------------------
# 7. Add panel letters
# --------------------------------------------------

for ax, letter in zip(
    axes.flat,
    ["A", "B", "C", "D"]
):
    ax.text(
        -0.13,
        1.08,
        letter,
        transform=ax.transAxes,
        fontsize=17,
        fontweight="bold",
        va="top"
    )

# --------------------------------------------------
# 8. Add shared legend
# --------------------------------------------------

handles, labels = axes[0, 0].get_legend_handles_labels()

fig.legend(
    handles,
    labels,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.99),
    ncol=2,
    frameon=False,
    fontsize=11
)

# --------------------------------------------------
# 9. Adjust spacing
# --------------------------------------------------

fig.subplots_adjust(
    top=0.88,
    bottom=0.08,
    left=0.10,
    right=0.97,
    hspace=0.38,
    wspace=0.30
)

# --------------------------------------------------
# 10. Save the new figure
# --------------------------------------------------

pdf = results / "PDAC_four_panel_protein_adjusted.pdf"
png = results / "PDAC_four_panel_protein_adjusted.png"

fig.savefig(pdf, bbox_inches="tight")
fig.savefig(png, dpi=300, bbox_inches="tight")

plt.close(fig)

print("SUCCESS: Protein-adjusted four-panel figure created.")
print("PDF:", pdf)
print("PNG:", png)


