# Integrative Analysis of Single-Cell Transcriptomic and Multilayer Signaling Networks in Glioma Reveal Tumor Progression Stage

## Overview
this paper introduces a framework for integrative analysis of single-cell transcriptomic data and multilayer signaling networks in glioma, shedding light on tumor progression stages. The repository contains scripts and resources necessary for reproducing the results described in the corresponding paper.

## Data Resource
The single-cell RNA sequencing (scRNA-seq) raw count matrices used in this analysis can be downloaded from the following link: [GSE89567](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi).

## Workflow
1. **scRNAseq Data Analysis using Seurat Package (R):**
   - The scRNAseq data analysis was performed using the Seurat package in R. 
   - Scripts for this part of the analysis are available in the `seurat` folder in this repository.

2. **Extraction of Inter and Intra Signaling Networks (R):**
   - Inter and intra-signaling networks were extracted for both grade III and grade IV gliomas.
   - R code for this part of the analysis is available in the `seurat` folder in this repository.

3. **Quadruple Extraction and Identification (Python):**
   - Quadruples were extracted and the most important quadruples were identified using machine learning models.
   - Python code for this part of the analysis is available in the `python notebook` folder in this repository.

4. **Survival Analysis (R):**
   - Survival analysis related R code is available in the `seurat` folder in this repository.

## Usage
- Refer to individual folders for specific scripts and data files.
- Follow the instructions provided in each script to reproduce the corresponding analysis steps.

## Citation

## Contact
For any questions or inquiries, please contact [kkavousi@ut.ac.ir].

