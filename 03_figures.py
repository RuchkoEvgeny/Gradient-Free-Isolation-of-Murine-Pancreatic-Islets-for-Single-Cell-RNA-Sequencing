Generate the data-derived figures for the pancreatic-islet scRNA-seq
manuscript comparing Ficoll and Gradient-free isolation workflows.

Input:
    data/processed/islets_annotated.h5ad
    data/mouse_level_qc.csv

Expected AnnData fields:
    obs["workflow"]   : Ficoll / Gradient-free
    obs["cell_type"]  : final cell annotation
    obsm["X_umap"]    : final UMAP
    raw               : log-normalized expression for all genes

mouse_level_qc.csv columns:
    mouse_id, workflow, cell_yield, viability

Figures generated:
    Figure 5  cell yield and viability
    Figure 6  scRNA-seq QC metrics
    Figure 7  UMAPs and canonical marker expression
    Figure 8  cell composition and marker dot plot
    Figure 9  stress-associated expression and stress score

Workflow schematics and microscopy images are not generated
computationally.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


DATA_DIR = Path("data")
FIG_DIR = Path("figures")
FIG_DIR.mkdir(exist_ok=True)

ADATA_FILE = DATA_DIR / "processed" / "islets_annotated.h5ad"
MOUSE_QC_FILE = DATA_DIR / "mouse_level_qc.csv"

WORKFLOW_KEY = "workflow"
CELLTYPE_KEY = "cell_type"

WORKFLOW_ORDER = ["Ficoll", "Gradient-free"]

CELLTYPE_ORDER = [
    "Beta cells",
    "Alpha cells",
    "Delta cells",
    "PP cells",
    "Beta/alpha-like cells",
    "Ductal cells",
    "Acinar cells",
]

sc.settings.set_figure_params(
    dpi=120,
    dpi_save=300,
    frameon=False,
    fontsize=10,
    vector_friendly=True,
)


IDENTITY_GENES = [
    "Ins1",
    "Gcg",
    "Sst",
    "Ppy",
    "Prss2",
    "Krt19",
]

MARKER_PANEL = {
    "Beta cells": ["Ins1", "Ins2", "Pdx1", "Mafa", "Nkx6-1"],
    "Alpha cells": ["Gcg", "Ttr", "Mafb"],
    "Delta cells": ["Sst", "Hhex", "Rbp4"],
    "PP cells": ["Ppy"],
    "Acinar cells": ["Prss2", "Cpa1"],
    "Ductal cells": ["Krt19", "Sox9"],
}

SELECTED_STRESS_GENES = [
    "Fos",
    "Jun",
    "Atf3",
    "Dusp1",
    "Hspa1a",
    "Ddit3",
]

STRESS_PANEL = {
    "Immediate early response": [
        "Fos", "Fosb", "Jun", "Junb",
        "Jund", "Atf3", "Egr1", "Dusp1",
    ],
    "Heat shock response": [
        "Hspa1a", "Hspa1b", "Hsp90aa1", "Hspb1",
    ],
    "ER stress": [
        "Ddit3", "Atf4", "Xbp1",
    ],
    "Oxidative stress": [
        "Hmox1", "Sod2",
    ],
}

STRESS_SCORE_GENES = [
    "Fos", "Fosb", "Jun", "Junb", "Jund",
    "Atf3", "Egr1", "Dusp1",
    "Hspa1a", "Hspa1b", "Hsp90aa1", "Hspb1",
    "Ddit3", "Atf4", "Xbp1",
    "Hmox1", "Sod2",
]


def save_figure(fig, name):
    fig.savefig(
        FIG_DIR / f"{name}.svg",
        bbox_inches="tight",
    )
    fig.savefig(
        FIG_DIR / f"{name}.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def present_genes(adata, genes):
    return [
        gene
        for gene in genes
        if gene in adata.raw.var_names
    ]


adata = sc.read_h5ad(ADATA_FILE)

for key in (WORKFLOW_KEY, CELLTYPE_KEY):
    if key not in adata.obs:
        raise KeyError(
            f"{key!r} is missing from adata.obs."
        )

if adata.raw is None:
    raise ValueError(
        "adata.raw is required for marker and stress-gene plots."
    )

if "X_umap" not in adata.obsm:
    raise ValueError(
        "Final UMAP coordinates are missing from adata.obsm."
    )

adata.obs[WORKFLOW_KEY] = pd.Categorical(
    adata.obs[WORKFLOW_KEY],
    categories=WORKFLOW_ORDER,
    ordered=True,
)

present_types = set(
    adata.obs[CELLTYPE_KEY].astype(str)
)

celltypes_present = [
    x for x in CELLTYPE_ORDER
    if x in present_types
]

adata.obs[CELLTYPE_KEY] = pd.Categorical(
    adata.obs[CELLTYPE_KEY],
    categories=celltypes_present,
    ordered=True,
)

adata.obs["workflow_celltype"] = (
    adata.obs[WORKFLOW_KEY].astype(str)
    + " | "
    + adata.obs[CELLTYPE_KEY].astype(str)
)


# =====================================================================
# Figure 5: cell yield and viability
# =====================================================================

if MOUSE_QC_FILE.exists():
    mouse_qc = pd.read_csv(MOUSE_QC_FILE)

    required = {
        "mouse_id",
        "workflow",
        "cell_yield",
        "viability",
    }

    missing = required.difference(mouse_qc.columns)

    if missing:
        raise KeyError(
            "mouse_level_qc.csv is missing: "
            + ", ".join(sorted(missing))
        )

    mouse_qc["workflow"] = pd.Categorical(
        mouse_qc["workflow"],
        categories=WORKFLOW_ORDER,
        ordered=True,
    )

    summary = (
        mouse_qc
        .groupby("workflow", observed=True)[
            ["cell_yield", "viability"]
        ]
        .agg(["mean", "std"])
        .reindex(WORKFLOW_ORDER)
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(7.2, 3.5),
    )

    for ax, variable, ylabel, panel in [
        (
            axes[0],
            "cell_yield",
            "Total cells after dissociation",
            "A",
        ),
        (
            axes[1],
            "viability",
            "Cell viability (%)",
            "B",
        ),
    ]:
        means = summary[variable]["mean"]
        sd = summary[variable]["std"]

        ax.bar(
            WORKFLOW_ORDER,
            means,
            yerr=sd,
            capsize=4,
        )

        for x, workflow in enumerate(WORKFLOW_ORDER):
            values = mouse_qc.loc[
                mouse_qc["workflow"] == workflow,
                variable,
            ].to_numpy()

            ax.scatter(
                np.full(values.size, x),
                values,
                zorder=3,
            )

        ax.set_ylabel(ylabel)
        ax.set_xlabel("")
        ax.set_title(panel)

    fig.tight_layout()
    save_figure(
        fig,
        "Figure5_cell_yield_viability",
    )

else:
    print(
        "Figure 5 skipped because data/mouse_level_qc.csv "
        "was not found."
    )


# =====================================================================
# Figure 6: scRNA-seq QC
# =====================================================================

qc_panels = [
    ("n_genes_by_counts", "Detected genes per cell", "A"),
    ("total_counts", "Total transcript counts per cell", "B"),
    ("pct_counts_mt", "Mitochondrial transcripts (%)", "C"),
    ("pct_counts_ribo", "Ribosomal transcripts (%)", "D"),
]

missing_qc = [
    key
    for key, _, _ in qc_panels
    if key not in adata.obs
]

if missing_qc:
    raise KeyError(
        "QC columns missing from adata.obs: "
        + ", ".join(missing_qc)
    )

fig, axes = plt.subplots(
    2,
    2,
    figsize=(7.6, 6.5),
)

for ax, (key, ylabel, panel) in zip(
    axes.flat,
    qc_panels,
):
    sc.pl.violin(
        adata,
        keys=key,
        groupby=WORKFLOW_KEY,
        order=WORKFLOW_ORDER,
        jitter=0.20,
        rotation=0,
        ax=ax,
        show=False,
    )

    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.set_title(panel)

fig.tight_layout()
save_figure(fig, "Figure6_scRNAseq_QC")


# =====================================================================
# Figure 7: UMAPs and canonical markers
# =====================================================================

for workflow in WORKFLOW_ORDER:
    sample = adata[
        adata.obs[WORKFLOW_KEY] == workflow
    ].copy()

    fig = plt.figure(figsize=(10, 9))

    grid = fig.add_gridspec(
        3,
        3,
        height_ratios=[1.35, 1, 1],
    )

    ax = fig.add_subplot(grid[0, :])

    sc.pl.umap(
        sample,
        color=CELLTYPE_KEY,
        ax=ax,
        show=False,
        legend_loc="right margin",
        size=12,
        title=workflow,
    )

    feature_axes = [
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
        fig.add_subplot(grid[1, 2]),
        fig.add_subplot(grid[2, 0]),
        fig.add_subplot(grid[2, 1]),
        fig.add_subplot(grid[2, 2]),
    ]

    for ax, gene in zip(
        feature_axes,
        IDENTITY_GENES,
    ):
        if gene not in adata.raw.var_names:
            ax.axis("off")
            ax.set_title(f"{gene} not detected")
            continue

        sc.pl.umap(
            sample,
            color=gene,
            use_raw=True,
            ax=ax,
            show=False,
            size=10,
            title=gene,
        )

    fig.tight_layout()

    safe_name = (
        workflow.lower()
        .replace("-", "_")
        .replace(" ", "_")
    )

    save_figure(
        fig,
        f"Figure7_{safe_name}_cell_identity",
    )


# =====================================================================
# Figure 8A: relative cell composition
# =====================================================================

composition = (
    adata.obs
    .groupby(
        [WORKFLOW_KEY, CELLTYPE_KEY],
        observed=False,
    )
    .size()
    .unstack(fill_value=0)
)

composition = (
    composition
    .div(composition.sum(axis=1), axis=0)
    * 100
)

composition = composition.reindex(WORKFLOW_ORDER)

fig, ax = plt.subplots(figsize=(7.2, 4.2))

composition.plot(
    kind="bar",
    stacked=True,
    ax=ax,
)

ax.set_ylabel("Cells (%)")
ax.set_xlabel("")
ax.set_ylim(0, 100)

ax.legend(
    title="Cell type",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False,
)

plt.xticks(rotation=0)
fig.tight_layout()

save_figure(fig, "Figure8A_cell_composition")

composition.to_csv(
    FIG_DIR / "Figure8A_cell_composition.csv"
)


# =====================================================================
# Figure 8B: canonical marker dot plot
# =====================================================================

marker_panel = {
    group: present_genes(adata, genes)
    for group, genes in MARKER_PANEL.items()
}

marker_panel = {
    group: genes
    for group, genes in marker_panel.items()
    if genes
}

observed_groups = set(
    adata.obs["workflow_celltype"].astype(str)
)

group_order = []

for workflow in WORKFLOW_ORDER:
    for cell_type in CELLTYPE_ORDER:
        label = f"{workflow} | {cell_type}"
        if label in observed_groups:
            group_order.append(label)

adata.obs["workflow_celltype"] = pd.Categorical(
    adata.obs["workflow_celltype"],
    categories=group_order,
    ordered=True,
)

dotplot = sc.pl.dotplot(
    adata,
    var_names=marker_panel,
    groupby="workflow_celltype",
    use_raw=True,
    standard_scale="var",
    return_fig=True,
    show=False,
)

dotplot.savefig(
    FIG_DIR / "Figure8B_marker_dotplot.svg"
)

dotplot.savefig(
    FIG_DIR / "Figure8B_marker_dotplot.png",
    dpi=300,
)


# =====================================================================
# Figure 9: stress-associated transcriptional features
# =====================================================================

stress_genes = present_genes(
    adata,
    STRESS_SCORE_GENES,
)

if not stress_genes:
    raise ValueError(
        "None of the stress-response genes were found in adata.raw."
    )

missing_stress = sorted(
    set(STRESS_SCORE_GENES) - set(stress_genes)
)

if missing_stress:
    print(
        "Stress genes not found:",
        ", ".join(missing_stress),
    )


# Calculate the score from log-normalized, unscaled expression.
score_adata = adata.raw.to_adata()
score_adata.obs = adata.obs.copy()

sc.tl.score_genes(
    score_adata,
    gene_list=stress_genes,
    score_name="stress_response_score",
    random_state=0,
)

adata.obs["stress_response_score"] = (
    score_adata.obs["stress_response_score"]
    .reindex(adata.obs_names)
)


# Figure 9A-F: selected stress genes
fig, axes = plt.subplots(
    2,
    3,
    figsize=(10, 6.2),
)

for ax, gene, panel in zip(
    axes.flat,
    SELECTED_STRESS_GENES,
    list("ABCDEF"),
):
    if gene not in adata.raw.var_names:
        ax.axis("off")
        ax.set_title(f"{panel}  {gene} not detected")
        continue

    sc.pl.violin(
        adata,
        keys=gene,
        groupby=WORKFLOW_KEY,
        order=WORKFLOW_ORDER,
        use_raw=True,
        stripplot=False,
        ax=ax,
        show=False,
    )

    ax.set_xlabel("")
    ax.set_ylabel("Normalized expression")
    ax.set_title(f"{panel}  {gene}")

fig.tight_layout()
save_figure(fig, "Figure9A-F_stress_genes")


# Figure 9G: combined stress-response score
fig, ax = plt.subplots(figsize=(4.3, 4.0))

sc.pl.violin(
    adata,
    keys="stress_response_score",
    groupby=WORKFLOW_KEY,
    order=WORKFLOW_ORDER,
    stripplot=False,
    ax=ax,
    show=False,
)

ax.set_xlabel("")
ax.set_ylabel("Stress-response score")
ax.set_title("G")

fig.tight_layout()
save_figure(
    fig,
    "Figure9G_stress_response_score",
)

stress_summary = (
    adata.obs
    .groupby(WORKFLOW_KEY, observed=True)[
        "stress_response_score"
    ]
    .agg(["mean", "median"])
)

stress_summary.to_csv(
    FIG_DIR / "Figure9G_stress_score_summary.csv"
)

print("\nStress-response score:")
print(stress_summary)


# Figure 9H: stress genes by workflow and cell type
stress_panel = {
    group: present_genes(adata, genes)
    for group, genes in STRESS_PANEL.items()
}

stress_panel = {
    group: genes
    for group, genes in stress_panel.items()
    if genes
}

dotplot = sc.pl.dotplot(
    adata,
    var_names=stress_panel,
    groupby="workflow_celltype",
    use_raw=True,
    standard_scale="var",
    return_fig=True,
    show=False,
)

dotplot.savefig(
    FIG_DIR / "Figure9H_stress_dotplot.svg"
)

dotplot.savefig(
    FIG_DIR / "Figure9H_stress_dotplot.png",
    dpi=300,
)


# Export the number of cells used for each workflow/cell-type group.
cells_used = (
    adata.obs
    .groupby(
        [WORKFLOW_KEY, CELLTYPE_KEY],
        observed=True,
    )
    .size()
    .rename("n_cells")
    .reset_index()
)

cells_used.to_csv(
    FIG_DIR / "cells_used_for_figures.csv",
    index=False,
)

print(
    f"\nFinished. Figures saved to: "
    f"{FIG_DIR.resolve()}"
)
