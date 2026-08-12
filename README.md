# Single-cell-and-spatial-transcriptomic-analysis-of-SCLC
This repository contains the code used to process and visualize single-cell RNA sequencing (scRNA-seq) and Xenium spatial transcriptomic data for the study by Pan et al., "MYC drives neuroendocrine transformation of adenocarcinoma into RB1-proficient small cell lung cancer." In this work, scRNA-seq was performed using the 10x Genomics Chromium platform on 7 human small cell lung cancer (SCLC) samples and 2 primary lung tumors from KPM mice, and Xenium spatial transcriptomics was performed on 7 human combined SCLC/LUAD and 1 human transformed SCLC tumor sections. The code supports key analyses including scRNA-seq quality control, cell type annotation and proportion analysis, transcription factor activity analysis, and pseudotime trajectory reconstruction, as well as Xenium-based spatial cell-type annotation, epithelial lineage subtyping, and single-molecule-resolution transcript mapping.

# Contents of repository

**1. Human scRNA-seq data analysis**

**Fig1b-e & Extended data Fig1.R**

This R script contains the code used to generate Fig. 1b-e and Extended Data Fig. 1, including: preprocessing of scRNA-seq data from human SCLC clinical samples; copy number variation analysis with inferCNV; UMAP visualization of all cells and tumor subsets; bubble plot visualization of selected gene expression; proportion analysis of clusters and cell types; and preprocessing/visualization for pySCENIC transcription factor analysis. 

**pyscenic analysis--Fig1e.sh**

This script provides a standard pySCENIC workflow for transcription factor analysis, following the official tutorial (https://pyscenic.readthedocs.io/en/latest/installation.html). 

**2. Mouse scRNA-seq data analysis**

**Seurat analysis--Fig3a & Extended data Fig7a-e, h-i.R**

This R script reproduces Fig. 3a and Extended data Fig. 7a-e, h, including preprocessing of scRNA-seq data from 2 KPM tumors, copy number variation analysis using the inferCNV package, UMAP visualization of all and tumor cells, Bubble plot visualization of selected gene expression, cluster and cell-type proportion analysis, GSVA of hallmark signatures enriched in each tumor cell cluster, and preprocessing/visualization for pySCENIC transcription factor analysis. 

**Scanpy analysis--Fig3b-d, Extended data Fig7f-g & Extended data Fig13a.ipynb**

This Python script reproduces Fig. 3b-d, Extended data Fig. 7f-g, and Extended data Fig. 13a, including preprocessing of tumor cell scRNA-seq data defined by Seurat, and CellRank analysis to infer pseudotime ordering and expression dynamics of selected genes. 

**3. Human Xenium spatial transcriptomic analysis**

**Fig1f-k & Extended data Fig3.R**

This R script reproduces Fig. 1f-k and Extended data Fig. 3, including: loading and sample assignment of 10x Xenium data from human combined SCLC/LUAD and transformed SCLC tumor sections; data sketching and Harmony-based sample integration; cell-type annotation and epithelial subtype classification; spatial and UMAP visualization of cell types; and spatial mapping of marker gene expression at single-molecule resolution.

# Data availability
The raw and processed scRNA-seq data generated in this study is available in the Genome Sequence Archive (GSA) database under accession number HRA008577 and CRA019000. 

# Contact
For further inquiries regarding the analytical methods described in this study, please contact deshui.jia@shgh.cn. Additional relevant data can be obtained from the authors upon reasonable request.
