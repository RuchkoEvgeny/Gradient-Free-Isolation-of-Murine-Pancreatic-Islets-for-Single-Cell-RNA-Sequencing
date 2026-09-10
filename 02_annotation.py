Marker-guided annotation of pancreatic endocrine and exocrine cells.

Input:
    data/processed/islets_preprocessed.h5ad

Output:
    data/processed/islets_annotated.h5ad
    results/annotation_marker_dotplot.svg
    results/cluster_markers.csv
    results/cluster_annotation.csv
    results/cell_annotations.tsv.gz

The exact final cluster-number mapping from the original notebook is not
hard-coded here because it was not preserved unambiguously. To keep the
script runnable, clusters are assigned automatically from canonical
marker-module scores. Optional manual overrides can be supplied in
MANUAL_CLUSTER_OVERRIDES without changing the rest of the pipeline.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


DATA_DIR = Path("data")
RESULT_DIR = Path("results")
RESULT_DIR.mkdir(exist_ok=True)

INPUT_FILE = DATA_DIR / "processed" / "islets_preprocessed.h5ad"
OUTPUT_FILE = DATA_DIR / "processed" / "islets_annotated.h5ad"

CLUSTER_KEY = "leiden_endo_exo_res1_0"
CELLTYPE_KEY = "cell_type"

RANDOM_STATE = 0


# Canonical markers used in the manuscript and in the original analysis.
MARKERS = {
    "Beta cells": [
        "Ins1", "Ins2", "Pdx1", "Mafa", "Nkx6-1"
    ],
    "Alpha cells": [
        "Gcg", "Ttr", "Mafb"
    ],
    "Delta cells": [
        "Sst", "Hhex", "Rbp4"
    ],
    "PP cells": [
        "Ppy"
    ],
    "Acinar cells": [
        "Prss2", "Cpa1"
    ],
    "Ductal cells": [
        "Krt19", "Sox9", "Spp1"
    ],
}

# Optional manual corrections. Leave empty for a fully automatic run.
# Example:
# MANUAL_CLUSTER_OVERRIDES = {"3": "Beta/alpha-like cells"}
MANUAL_CLUSTER_OVERRIDES = {}


def genes_present(adata, genes):
    return [
        gene
        for gene in genes
        if gene in adata.raw.var_names
    ]


def score_marker_modules(adata):
    """Calculate cell-level scores for the major pancreatic lineages."""
    for cell_type, genes in MARKERS.items():
        found = genes_present(adata, genes)

        if not found:
            raise ValueError(
                f"No marker genes found for {cell_type}: {genes}"
            )

        score_name = (
            "score_"
            + cell_type.lower()
            .replace(" ", "_")
            .replace("/", "_")
        )

        sc.tl.score_genes(
            adata,
            gene_list=found,
            score_name=score_name,
            use_raw=True,
            random_state=RANDOM_STATE,
        )


def automatic_cluster_annotation(adata):
    """
    Assign cluster identities from mean marker-module scores.

    A cluster is labelled Beta/alpha-like when Beta and Alpha are the
    two strongest endocrine programs and their standardized scores are
    close to one another. This provides a reproducible fallback when an
    exact historical manual cluster-number mapping is unavailable.
    """
    score_columns = {
        "Beta cells": "score_beta_cells",
        "Alpha cells": "score_alpha_cells",
        "Delta cells": "score_delta_cells",
        "PP cells": "score_pp_cells",
        "Acinar cells": "score_acinar_cells",
        "Ductal cells": "score_ductal_cells",
    }

    means = (
        adata.obs
        .groupby(CLUSTER_KEY, observed=True)[
            list(score_columns.values())
        ]
        .mean()
    )

    # Standardize each marker program across clusters so that programs
    # containing different numbers of genes are more comparable.
    standardized = means.copy()

    for column in standardized.columns:
        values = standardized[column]
        sd = values.std(ddof=0)

        if sd == 0 or np.isnan(sd):
            standardized[column] = 0.0
        else:
            standardized[column] = (
                values - values.mean()
            ) / sd

    reverse = {
        column: cell_type
        for cell_type, column in score_columns.items()
    }

    mapping = {}

    for cluster, row in standardized.iterrows():
        ranked = row.sort_values(ascending=False)
        top_col = ranked.index[0]
        top_type = reverse[top_col]

        beta = row["score_beta_cells"]
        alpha = row["score_alpha_cells"]

        top_two = set(ranked.index[:2])

        beta_alpha_are_top = top_two == {
            "score_beta_cells",
            "score_alpha_cells",
        }

        beta_alpha_close = abs(beta - alpha) <= 0.75

        if beta_alpha_are_top and beta_alpha_close:
            label = "Beta/alpha-like cells"
        else:
            label = top_type

        mapping[str(cluster)] = label

    mapping.update(
        {
            str(cluster): label
            for cluster, label
            in MANUAL_CLUSTER_OVERRIDES.items()
        }
    )

    return mapping, means, standardized


adata = sc.read_h5ad(INPUT_FILE)

if CLUSTER_KEY not in adata.obs:
    raise KeyError(
        f"{CLUSTER_KEY!r} is not present in adata.obs."
    )

if adata.raw is None:
    raise ValueError(
        "adata.raw is required for marker-based annotation."
    )


# ---------------------------------------------------------------------
# Marker visualization and cluster marker table
# ---------------------------------------------------------------------

marker_panel = {
    group: genes_present(adata, genes)
    for group, genes in MARKERS.items()
}

marker_panel = {
    group: genes
    for group, genes in marker_panel.items()
    if genes
}

dotplot = sc.pl.dotplot(
    adata,
    var_names=marker_panel,
    groupby=CLUSTER_KEY,
    use_raw=True,
    standard_scale="var",
    return_fig=True,
    show=False,
)

dotplot.savefig(
    RESULT_DIR / "annotation_marker_dotplot.svg",
)

dotplot.savefig(
    RESULT_DIR / "annotation_marker_dotplot.png",
    dpi=300,
)


sc.tl.rank_genes_groups(
    adata,
    groupby=CLUSTER_KEY,
    method="wilcoxon",
    use_raw=True,
)

cluster_markers = sc.get.rank_genes_groups_df(
    adata,
    group=None,
)

cluster_markers.to_csv(
    RESULT_DIR / "cluster_markers.csv",
    index=False,
)


# ---------------------------------------------------------------------
# Reproducible marker-based cluster annotation
# ---------------------------------------------------------------------

score_marker_modules(adata)

cluster_mapping, module_means, module_zscores = (
    automatic_cluster_annotation(adata)
)

adata.obs[CELLTYPE_KEY] = (
    adata.obs[CLUSTER_KEY]
    .astype(str)
    .map(cluster_mapping)
    .astype("category")
)


annotation_table = pd.DataFrame(
    {
        "cluster": list(cluster_mapping.keys()),
        "cell_type": list(cluster_mapping.values()),
    }
).sort_values("cluster")

annotation_table.to_csv(
    RESULT_DIR / "cluster_annotation.csv",
    index=False,
)

module_means.to_csv(
    RESULT_DIR / "cluster_marker_module_scores.csv"
)

module_zscores.to_csv(
    RESULT_DIR / "cluster_marker_module_zscores.csv"
)


# ---------------------------------------------------------------------
# Save final object and GEO-ready cell annotation table
# ---------------------------------------------------------------------

adata.write_h5ad(OUTPUT_FILE)

cell_annotations = adata.obs[
    ["workflow", CLUSTER_KEY, CELLTYPE_KEY]
].copy()

cell_annotations.index.name = "barcode"

cell_annotations.to_csv(
    RESULT_DIR / "cell_annotations.tsv.gz",
    sep="\t",
    compression="gzip",
)

print("\nCluster annotation:")
print(annotation_table.to_string(index=False))

print("\nFinal cell counts:")
print(
    pd.crosstab(
        adata.obs["workflow"],
        adata.obs[CELLTYPE_KEY],
    )
)

print(f"\nSaved: {OUTPUT_FILE}")
