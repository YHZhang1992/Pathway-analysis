# Environment: Pathway & enrichment analysis

## Standard profile

The canonical `workflow.py` portable profile uses Python 3.10+ and the standard library only, so its example and smoke test run without package installation. It establishes input validation, deterministic processing, output structure, figures, and provenance. Use the advanced implementation for production-grade domain methods.

Do not install the union of every historical dependency. Select the advanced canonical entry point, separate standard-library/local imports from external packages, and create a per-workflow lockfile.

Detected imports (static, may include local modules or command tokens):

`AnnotationDbi, Biobase, DESeq2, GEOparse, GEOquery, GSVA, KEGGREST, Matrix, ReactomePA, Rscript, Seurat, SeuratObject, __future__, apeglm, argparse, cd, clusterProfiler, collections, data.table, dplyr, echo, emmeans, export, fgsea, forcats, ggplot2, ggrepel, grid, gridExtra, gseapy, hashlib, json, jsonlite, limma, lmerTest, mashr, math, matplotlib, msigdbr, mygene, numpy, org.Hs.eg.db, pandas, patchwork, pathlib, pheatmap, pkg, plotly, prepare_slide_assets, pydeseq2, python, re, readr, readxl, report_pdf_customizations, scipy, seaborn, set, shutil, sklearn, statistics, statsmodels, sva, tidyverse, topGO, variancePartition, yaml`

- Python: declare dependencies and interpreter range in `pyproject.toml`; create and commit `uv.lock`.
- R: initialize `renv`, install approved CRAN/Bioconductor versions, then commit `renv.lock`.
- Hybrid/bioinformatics: use the supplied Conda manifest where available and record STAR, samtools, Cell Ranger, reference genome, annotation, and gene-set versions.

## Controlled/GxP profile

No constrained-package implementation is supplied for this topic.

Never install or update packages inside an analysis run. Qualification must pin exact versions on the target platform and retain installation evidence.
