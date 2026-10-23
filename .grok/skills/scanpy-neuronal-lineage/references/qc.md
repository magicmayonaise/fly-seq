# Load, QC, doublets

Source: `docs/tutorials/basics/clustering.ipynb` in the checkout named in `dataset.md`, adapted only for this matrix layout and FlyBase names.

Run from the fly-seq repo root with `uv run`. About 190,624 barcodes × 16,974 genes before filtering. Do the `min_genes` filter before neighbors.

```python
from pathlib import Path

import anndata as ad
import scanpy as sc

RAW = Path(
    "raw/scRNAseq_neuronal-lineage. Allen et al. Cell Genomics. 2026. /GSE296540_RAW"
)
OUT = Path("processed")
OUT.mkdir(exist_ok=True)
(OUT / "figures").mkdir(exist_ok=True)
sc.settings.figdir = OUT / "figures"

SAMPLES = {
    "FemaleBrainRep1": ("GSM8972280", "female", "1"),
    "FemaleBrainRep2": ("GSM8972281", "female", "2"),
    "MaleBrainRep1": ("GSM8972282", "male", "1"),
    "MaleBrainRep2": ("GSM8972283", "male", "2"),
}
PANEL = [
    "fru", "fru_plus", "fru_dP1", "fru_dP1_plus", "dsx",
    "EGFP", "Gal4", "dsxM", "dsxF", "dsxCOM",
    "fruA", "fruB", "fruC", "fruD1", "fruD2",
    "fruABint", "fruBCint", "fruCint", "fruPr1Ex1", "fruPr1Ex2", "fruCOM",
]

adatas = {}
for sample, (gsm, sex, rep) in SAMPLES.items():
    a = sc.read_10x_mtx(RAW, prefix=f"{gsm}_{sample}_")
    a.obs["sex"] = sex
    a.obs["replicate"] = rep
    a.obs["gsm"] = gsm
    adatas[sample] = a

adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()
adata.var["panel"] = adata.var_names.isin(PANEL)
adata.var["mt"] = adata.var_names.str.startswith("mt:")
adata.var["ribo"] = adata.var_names.str.startswith(("RpL", "RpS"))

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True)
adata.write_h5ad(OUT / "gse296540_with_qc.h5ad")

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    groupby="sample",
    jitter=0.4,
    multi_panel=True,
    save="_qc_by_sample.png",
)
sc.pl.scatter(
    adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
    save="_qc_scatter.png",
)
```

Male libraries have about twice the barcodes of the female libraries. Read the violins per sample before choosing a cutoff past the tutorial's.

```python
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
```

Scrublet wants unnormalised counts (`sc.pp.scrublet` docstring). Run it now, with `batch_key="sample"`. The tutorial default `expected_doublet_rate` is 0.05. Plot the score. Keep `predicted_doublet` on the object. Dropping doublets waits until they are visible on the UMAP.

```python
sc.pp.scrublet(adata, batch_key="sample")
sc.pl.scrublet_score_distribution(adata, save="_scrublet.png")
adata.write_h5ad(OUT / "gse296540_filtered.h5ad")
```

Report, per sample: cells in, cells after `min_genes=100`, median `total_counts`, median `n_genes_by_counts`, median `pct_counts_mt`, and the doublet fraction. A mito hard cutoff is a separate decision after those plots. Neurons often sit higher than other tissues, and this tutorial does not apply one.
