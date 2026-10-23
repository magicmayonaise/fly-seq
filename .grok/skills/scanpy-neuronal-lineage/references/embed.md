# Normalize, embed, markers

Continue from `processed/gse296540_filtered.h5ad`. Counts are still in `X`. Same tutorial as `qc.md`.

```python
from pathlib import Path

import scanpy as sc

OUT = Path("processed")
sc.settings.figdir = OUT / "figures"
adata = sc.read_h5ad(OUT / "gse296540_filtered.h5ad")

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)  # target_sum=None → median library size
sc.pp.log1p(adata)

# ScanpyV1 flavor is "seurat", which expects log data. n_top_genes is set.
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
sc.pl.highly_variable_genes(adata, save="_hvg.png")

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save="_pca_variance.png")
sc.pl.pca(
    adata,
    color=["sample", "sex", "pct_counts_mt"],
    dimensions=[(0, 1), (0, 1), (0, 1)],
    ncols=3,
    size=2,
    save="_pca_qc.png",
)

sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["sample", "sex"], size=2, save="_umap_sample.png")

# Tutorial overrides the ScanpyV1 leidenalg default.
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
for res in (0.2, 0.5, 1.0):
    sc.tl.leiden(
        adata,
        key_added=f"leiden_res_{res:.2f}",
        resolution=res,
        flavor="igraph",
        n_iterations=2,
    )

sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    wspace=0.5,
    size=3,
    save="_umap_doublets.png",
)
sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
    save="_umap_qc.png",
)
adata.write_h5ad(OUT / "gse296540_embedded.h5ad")
```

If a PC tracks `pct_counts_mt` or `sample` more than biology, say so before clustering and ask whether to drop that component. Do not drop it silently.

## Markers

`rank_genes_groups` expects log data. Pass `method="wilcoxon"` because the ScanpyV1 preset is `t-test`. The function docstring warns that cell-level p-values are inflated (cells are not independent) and points at pseudobulk for a real between-sample test. Use Wilcoxon here to propose cluster markers. For female vs male with two replicates each, aggregate with `sc.get.aggregate` before any claim that a gene differs by sex.

```python
PANEL = adata.var_names[adata.var["panel"]].tolist()
sc.pl.dotplot(adata, PANEL, groupby="leiden_res_0.50", standard_scale="var", save="_panel.png")

sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=5,
    save="_markers.png",
)
adata.write_h5ad(OUT / "gse296540_embedded.h5ad")
```

Name a cluster only from genes that are actually high in that cluster, including the panel when it is among them. Do not assign hemilineages from memory of the paper. The paper's labels live in the Seurat RDS and at flycns.com; use those as a lookup the user can request, not as labels to paste on.
