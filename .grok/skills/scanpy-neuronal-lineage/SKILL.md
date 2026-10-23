---
name: scanpy-neuronal-lineage
description: >
  Analyze the Allen et al. Cell Genomics 2026 Drosophila central-brain
  scRNA-seq matrices in this repo (GEO GSE296540, four 10x mtx samples)
  with scanpy. Use when the user says scanpy, QC, filter cells, doublets,
  scrublet, normalize, HVG, PCA, neighbors, UMAP, leiden, markers,
  rank_genes_groups, fru, dsx, neuronal lineage, GSE296540, FemaleBrain,
  or MaleBrain, or runs /scanpy-neuronal-lineage. Read references/dataset.md
  before any code. API source is the local scanpy checkout named in that file.
---

# Scanpy analysis — GSE296540 neuronal lineage

Run the scanpy basics tutorial on this repo's matrices. The dataset facts live in `references/dataset.md`. Do not restate them elsewhere.

## Before any code

1. Read `references/dataset.md` in this skill directory.
2. Read only the phase reference you are about to run:
   - Load, QC, or doublets: `references/qc.md`
   - Normalize through markers: `references/embed.md`
3. Confirm the import is the checkout named in the dataset reference (`uv run python -c "import scanpy,inspect; print(scanpy.__version__, inspect.getfile(scanpy))"`). If it is missing or a different tree, install that checkout as the dataset reference says. Do not code against a different scanpy.

## Order

Stop after the phase the user asked for. Write plots and `.h5ad` files under `processed/`, never back into `raw/`.

1. Load the four matrices and concatenate.
2. QC metrics and plots, grouped by `sample`.
3. Permissive cell and gene filters from the tutorial.
4. Scrublet on raw counts, with doublets flagged.
5. Store counts, then normalize and log1p.
6. HVG, PCA, neighbors, UMAP, Leiden.
7. Marker genes and the fru/dsx panel.

Review the QC and doublet plots before a tighter filter than the tutorial's `min_genes=100` / `min_cells=3`. Sex comparisons of single cells are exploratory; the rank-genes docstring in this checkout says cells are not independent replicates.

## What not to import from the Anthropic QC skill

`anthropics/life-sciences` `single-cell-rna-qc` expects `.h5ad` or 10x `.h5`, human `MT-` / mouse `mt-` / hemoglobin patterns, and its own MAD script. This dataset is gzipped mtx in one directory, Drosophila `mt:` genes, and no hemoglobin. Follow this checkout's `docs/tutorials/basics/clustering.ipynb` instead. That tutorial does not use MAD.
