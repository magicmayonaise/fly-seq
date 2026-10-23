# GSE296540 dataset contract

Checked 2026-10-21 against the files on disk and scanpy checkout `1.13.0a2-31-gda45762d` at `/Users/dlandayan/00_Developer/tools/scanpy` (commit `da45762d`).

## Study

Folder name: `scRNAseq_neuronal-lineage. Allen et al. Cell Genomics. 2026.`

GEO series [GSE296540](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296540), *Drosophila melanogaster*, title "A High-Resolution Map Predicts Lineage and Birth Order Underlying Neuronal Identity". The same accession is the data for Allen et al., Cell Genomics 2026 (sexual dimorphism in the adult brain; DOI 10.1016/j.xgen.2025.101125). Deposited author objects are Seurat RDS. This skill reanalyzes the mtx matrices with scanpy. It does not reproduce the Seurat pipeline.

Author code and viewer, for later comparison: [aaron-allen/Dmel-adult-central-brain-atlas](https://github.com/aaron-allen/Dmel-adult-central-brain-atlas), [flycns.com](https://www.flycns.com/).

## Files

All four samples sit in one directory:

`raw/scRNAseq_neuronal-lineage. Allen et al. Cell Genomics. 2026. /GSE296540_RAW`

| sample | GSM | sex | rep | barcodes | prefix for `sc.read_10x_mtx` |
|---|---|---|---|---|---|
| FemaleBrainRep1 | GSM8972280 | female | 1 | 32709 | `GSM8972280_FemaleBrainRep1_` |
| FemaleBrainRep2 | GSM8972281 | female | 2 | 35879 | `GSM8972281_FemaleBrainRep2_` |
| MaleBrainRep1 | GSM8972282 | male | 1 | 62211 | `GSM8972282_MaleBrainRep1_` |
| MaleBrainRep2 | GSM8972283 | male | 2 | 59825 | `GSM8972283_MaleBrainRep2_` |

Each sample is Cell Ranger v3-style gzipped mtx: `{prefix}matrix.mtx.gz`, `{prefix}features.tsv.gz`, `{prefix}barcodes.tsv.gz`. `read_10x_mtx(path, prefix=...)` with `compressed=True` (the default) is the reader in `src/scanpy/io/_read.py`.

16,974 features, all `Gene Expression`. Column 1 and column 2 of `features.tsv` are the same symbol, not FlyBase `FBgn` ids. Default `var_names="gene_symbols"` is correct.

## Gene names that differ from the human tutorial

The tutorial sets `var_names.str.startswith("MT-")`, `("RPS","RPL")`, and `^HB[^(P)]`. For this matrix:

- Mitochondrial: FlyBase prefix `mt:` (`mt:CoI`, `mt:ND2`, `mt:Cyt-b`, mitochondrial tRNAs). Use `startswith("mt:")`. `Cyt-b5` and `Cyt-b5-r` are nuclear.
- Ribosomal: `RpL` (54) and `RpS` (38). Use `startswith(("RpL", "RpS"))`.
- Hemoglobin: none. Do not create an `hb` qc var.

Rows 1–21 are the sex-determination and transgene panel, then the normal transcriptome starts (`Nep3`, `CG…`). Keep the panel in the matrix. Mark it so plots can find it:

`fru`, `fru_plus`, `fru_dP1`, `fru_dP1_plus`, `dsx`, `EGFP`, `Gal4`, `dsxM`, `dsxF`, `dsxCOM`, `fruA`, `fruB`, `fruC`, `fruD1`, `fruD2`, `fruABint`, `fruBCint`, `fruCint`, `fruPr1Ex1`, `fruPr1Ex2`, `fruCOM`.

## Install that matches the checkout

From the fly-seq repo:

```bash
uv add --editable "/Users/dlandayan/00_Developer/tools/scanpy[leiden,scrublet]"
```

`leiden` pulls `igraph` and `leidenalg`. `scrublet` pulls `scikit-image`. Scrublet itself is in this checkout's `sc.pp.scrublet`. Default `sc.settings.preset` is `ScanpyV1` (`src/scanpy/_settings/__init__.py`): HVG flavor `seurat`, Leiden flavor `leidenalg`, rank-genes method `t-test`. Pass the tutorial's explicit arguments anyway.

## Outputs

`processed/` at the repo root. Figures in `processed/figures/`.
