# Single-cell-RNA-sequencing-of-SCLC
This repository contains the code used to process and visualize single-cell RNA sequencing (scRNA-seq) data for the study by Pan et al., “MYC drives neuroendocrine transformation of adenocarcinoma into RB1-proficient small cell lung cancer.” In this work, scRNA-seq was performed using the 10x Genomics Chromium platform on both human small cell lung cancer (SCLC) samples and primary lung tumors from KPM mice. The code supports key analyses including scRNA-seq quality control, cell type annotation and proportion analysis, transcription factor activity analysis, and pseudotime trajectory reconstruction. 

**Contents of repository**
**1. Human scRNA-seq data analysis**
**Fig1 & Extended data Fig1.R**
This R script contains the code used to generate Fig. 1 and Extended Data Fig. 1, including: preprocessing of scRNA-seq data from human SCLC clinical samples; copy number variation analysis with inferCNV; UMAP visualization of all cells and tumor subsets; bubble plot visualization of selected gene expression; proportion analysis of clusters and cell types; and preprocessing/visualization for pySCENIC transcription factor analysis. 

**pyscenic analysis--Fig1c.sh**
This script provides a standard pySCENIC workflow for transcription factor analysis, following the official tutorial (https://pyscenic.readthedocs.io/en/latest/installation.html). 

**2. Mouse scRNA-seq data analysis**
**Seurat analysis--Fig3a & Extended data Fig5a-e, h-i.R**
This R script reproduces Fig. 3a and Extended data Fig. 5a-e, h, including preprocessing of scRNA-seq data from 2 KPM tumors, copy number variation analysis using the inferCNV package, UMAP visualization of all and tumor cells, Bubble plot visualization of selected gene expression, cluster and cell-type proportion analysis, GSVA of hallmark signatures enriched in each tumor cell cluster, and preprocessing/visualization for pySCENIC transcription factor analysis. 

**Scanpy analysis--Fig3b-d, Extended data Fig5f-g & Extended data Fig10a.ipynb**
This Python script reproduces Fig. 3b-d, Extended data Fig. 5f-g, and Extended data Fig. 10a, including preprocessing of tumor cell scRNA-seq data defined by Seurat, and CellRank analysis to infer pseudotime ordering and expression dynamics of selected genes. 


**Data availability**
The raw and processed scRNA-seq data generated in this study is available in the Genome Sequence Archive (GSA) database under accession number HRA008577 and CRA019000. 

**Contact**
For further inquiries regarding the analytical methods described in this study, please contact deshui.jia@shgh.cn. Additional relevant data can be obtained from the authors upon reasonable request.

