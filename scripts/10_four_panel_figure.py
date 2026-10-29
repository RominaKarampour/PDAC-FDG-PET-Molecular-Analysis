

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------
# PATHS
# --------------------------------------------------

project = Path(__file__).resolve().parent.parent
results = project / "results"

# --------------------------------------------------
# LOAD RESULTS
# --------------------------------------------------

rna = pd.read_csv(results / "psn1_hpac_rna_expression.csv")
protein = pd.read_csv(results / "psn1_hpac_protein_abundance.csv")
residuals = pd.read_csv(results / "e2f_g2m_adjusted_residuals.csv")

# --------------------------------------------------
# HELPER: READ VALUES FROM CSV
# --------------------------------------------------

def get_value(df, gene_column, gene, cell_line, value_column):
    row = df[
        (df[gene_column] == gene)
        & (df["CellLine"] == cell_line)
    ]

    if len(row) != 1:
        raise ValueError(
            f"Expected one row for {gene}, {cell_line}; found {len(row)}"
        )

    return float(row.iloc[0][value_column])


def find_column(df, candidates):
    for name in candidates:
        if name in df.columns:
            return name

    raise ValueError(
        f"Could not find any of {candidates}. "
        f"Available columns: {list(df.columns)}"
    )


# Identify columns in the expression tables
rna_gene_col = find_column(rna, ["Gene", "Symbol"])
rna_value_col = find_column(
    rna,
    ["RNA_log2_TPM_plus_1"]
)
protein_gene_col = find_column(
    protein,
    ["Protein", "Gene", "Symbol"]
)
protein_value_col = find_column(
    protein,
    ["Abundance", "ProteinAbundance", "Value"]
)

# --------------------------------------------------
# FIGURE SETTINGS
# --------------------------------------------------

fig, axes = plt.subplots(2, 2, figsize=(13, 10))

colors = {
    "PSN1": "#D62728",
    "HPAC": "#1F77B4",
}

width = 0.34

# --------------------------------------------------
# PANEL A: RNA EXPRESSION
# --------------------------------------------------

ax = axes[0, 0]

rna_genes = ["SLC2A1", "SLC2A3", "HK1", "HK2"]
x = np.arange(len(rna_genes))

for i, cell_line in enumerate(["PSN1", "HPAC"]):
    values = [
        get_value(
            rna,
            rna_gene_col,
            gene,
            cell_line,
            rna_value_col
        )
        for gene in rna_genes
    ]

    positions = x + (i - 0.5) * width

    bars = ax.bar(
        positions,
        values,
        width,
        label=cell_line,
        color=colors[cell_line]
    )

    ax.bar_label(bars, fmt="%.2f", padding=3, fontsize=8)

ax.set_xticks(x)
ax.set_xticklabels(rna_genes, fontstyle="italic")
ax.set_ylabel("RNA expression [log₂(TPM + 1)]")
ax.set_title("A. RNA expression", fontweight='bold')
ax.set_ylim(bottom=0)
ax.margins(x=0.12, y=0.18)
ax.legend(frameon=False)

# --------------------------------------------------
# PANEL B: PROTEIN ABUNDANCE
# --------------------------------------------------

ax = axes[0, 1]

proteins = ["GLUT1", "HK1", "HK2", "LDHA", "MCT4"]
x = np.arange(len(proteins))

for i, cell_line in enumerate(["PSN1", "HPAC"]):
    values = [
        get_value(
            protein,
            protein_gene_col,
            gene,
            cell_line,
            protein_value_col
        )
        for gene in proteins
    ]

    positions = x + (i - 0.5) * width

    bars = ax.bar(
        positions,
        values,
        width,
        label=cell_line,
        color=colors[cell_line]
    )

    ax.bar_label(bars, fmt="%.2f", padding=3, fontsize=8)

ax.set_xticks(x)
ax.set_xticklabels(proteins)
ax.set_ylabel("Relative protein abundance")
ax.set_title("B. Protein abundance",fontweight="bold" )
ax.set_ylim(bottom=0)
ax.margins(x=0.12, y=0.18)
ax.legend(frameon=False)

# --------------------------------------------------
# PANELS C AND D: RNA RESIDUALS
# --------------------------------------------------

for ax, pathway, title in [
    (axes[1, 0], "E2F", "C. E2F-adjusted RNA residuals"),
    (axes[1, 1], "G2M", "D. G2/M-adjusted RNA residuals"),
]:

    display_genes = ["SLC2A1", "HK2"]
    x = np.arange(len(display_genes))

    for i, cell_line in enumerate(["PSN1", "HPAC"]):

        values = []

        for gene in display_genes:

            # The residual CSV calls SLC2A1 "GLUT1"
            csv_gene = "GLUT1" if gene == "SLC2A1" else gene

            row = residuals[
                (residuals["Pathway"] == pathway)
                & (residuals["Gene"] == csv_gene)
                & (residuals["CellLine"] == cell_line)
            ]

            if len(row) != 1:
                raise ValueError(
                    f"Expected one residual for "
                    f"{pathway}, {csv_gene}, {cell_line}; "
                    f"found {len(row)}"
                )

            values.append(float(row.iloc[0]["Residual"]))

        positions = x + (i - 0.5) * width

        bars = ax.bar(
            positions,
            values,
            width,
            label=cell_line,
            color=colors[cell_line]
        )

        for bar, value in zip(bars, values):
            ax.annotate(
                f"{value:+.3f}",
                (
                    bar.get_x() + bar.get_width() / 2,
                    value
                ),
                xytext=(0, 5 if value >= 0 else -5),
                textcoords="offset points",
                ha="center",
                va="bottom" if value >= 0 else "top",
                fontsize=9
            )

    ax.axhline(0, color="black", linewidth=0.9)

    ax.set_xticks(x)
    ax.set_xticklabels(display_genes, fontstyle="italic")
    ax.set_ylabel("RNA residual [log₂(TPM + 1)]")
    ax.set_title(title, fontweight="bold")

    # Identical Y-axis range makes C and D comparable
    ax.set_ylim(-2.0, 2.0)
    ax.set_yticks([-2, -1, 0, 1, 2])

    # Automatically leave space for all grouped bars
    ax.margins(x=0.15)

    ax.legend(frameon=False, loc="upper right")

# --------------------------------------------------
# SAVE FIGURE
# --------------------------------------------------

fig.tight_layout(pad=2.5)

pdf_path = results / "PDAC_four_panel_molecular_figure.pdf"
png_path = results / "PDAC_four_panel_molecular_figure.png"

fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
fig.savefig(png_path, dpi=300, bbox_inches="tight")

plt.close(fig)

print("Saved:", pdf_path)
print("Saved:", png_path)

