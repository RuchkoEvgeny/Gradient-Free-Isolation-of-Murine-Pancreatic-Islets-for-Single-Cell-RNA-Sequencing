# Gradient-Free-Isolation-of-Murine-Pancreatic-Islets-for-Single-Cell-RNA-Sequencing
# Pancreatic islet scRNA-seq analysis

This repository contains the Python code used to reproduce the main
single-cell RNA-seq processing, cell-type annotation, and data-derived
figures for the manuscript:

**Gradient-Free Isolation of Murine Pancreatic Islets for Single-Cell RNA Sequencing**

The analysis compares two pancreatic islet isolation workflows:

- **Ficoll**
- **Gradient-free**

No historical internal treatment labels are used in the public scripts.

## Files

```text
01_preprocessing.py
02_annotation.py
03_figures.py
requirements.txt
```

## Analysis workflow

### 1. Preprocessing

`01_preprocessing.py` performs:

- sample-wise QC;
- removal of cells with fewer than 200 detected genes;
- removal of cells with fewer than 500 total counts;
- removal of cells with more than 10% mitochondrial transcripts;
- Scrublet doublet detection;
- normalization to 10,000 counts per cell;
- log1p transformation;
- initial analysis using 2,000 highly variable genes;
- selection of the endocrine/exocrine subset;
- repeated preprocessing of the subset using 3,000 highly variable genes;
- scaling with `max_value = 10`;
- PCA using 50 components;
- Harmony integration using workflow as the batch variable;
- nearest-neighbor graph construction with 30 neighbors and 30 Harmony-corrected PCs;
- UMAP with `min_dist = 0.3`;
- Leiden clustering at resolutions 0.5 and 1.0.

The script expects count-level AnnData files. If ambient-RNA correction
with CellBender was performed upstream, the corrected matrices should be
used as the input count matrices.

### 2. Cell-type annotation

`02_annotation.py` uses canonical pancreatic marker modules and cluster
marker genes to annotate the major populations:

- Beta cells
- Alpha cells
- Delta cells
- PP cells
- Beta/alpha-like cells
- Ductal cells
- Acinar cells

Because the exact historical cluster-number mapping was not preserved
unambiguously, the script contains a reproducible marker-score-based
fallback annotation. Optional manual cluster overrides can be entered in
`MANUAL_CLUSTER_OVERRIDES` if the original final mapping is available.

### 3. Figures

`03_figures.py` generates the data-derived manuscript figures:

- **Figure 5** — cell yield and viability;
- **Figure 6** — scRNA-seq QC metrics;
- **Figure 7** — UMAPs and canonical marker expression;
- **Figure 8** — cell composition and marker-gene dot plot;
- **Figure 9** — stress-associated gene expression, combined stress-response score, and stress-gene dot plot.

Workflow schematics and microscopy images are not generated
computationally.

## Expected input structure

```text
data/
├── ficoll_counts.h5ad
├── gradient_free_counts.h5ad
├── mouse_level_qc.csv
└── processed/
```

The `workflow` values used in metadata are:

```text
Ficoll
Gradient-free
```

`mouse_level_qc.csv` should contain:

```text
mouse_id,workflow,cell_yield,viability
```

Each row represents one mouse-derived preparation. Technical replicate
measurements should be averaged before entering this table.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The downstream analysis was performed with **Scanpy 1.12.1**.

## Running the analysis

Run the scripts from the repository root:

```bash
python 01_preprocessing.py
python 02_annotation.py
python 03_figures.py
```

Outputs are written to:

```text
data/processed/
results/
figures/
```

## Stress-response score

The combined stress-response score is calculated using
`scanpy.tl.score_genes` and the following genes:

```text
Fos, Fosb, Jun, Junb, Jund, Atf3, Egr1, Dusp1,
Hspa1a, Hspa1b, Hsp90aa1, Hspb1,
Ddit3, Atf4, Xbp1, Hmox1, Sod2
```

Comparisons between the two pooled scRNA-seq libraries are descriptive.

## Data availability

Raw sequencing data, processed expression matrices, cell annotations,
and metadata are intended to accompany the manuscript through the
associated public data record. The repository contains analysis code
rather than raw sequencing data.

## Citation

Please cite the associated manuscript when using this workflow or code.
Bibliographic details can be added to this section after publication.
