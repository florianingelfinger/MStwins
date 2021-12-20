# MStwins
This repository contains code that was used in the analysis of the manuscript Ingelfinger, Gerdes et al. "Twin study reveals non-heritable immune perturbations in multiple sclerosis". The purpose of this repository is to provide additional details beyond the method section of the manuscript for the interested reader. 

Additional code used for the CITE-seq analysis can be accessed here:
https://github.com/beltranLab/twin_study_Nature_2021/blob/main/CITE-Seq_analysis.ipynb

The code in this repository is split into the following sections:

- **CyTOF_preprocessing**: 
Input for this script are the individual debarcoded mass cytometry .fcs files. In this script individual files are concatenated, transformed using a hyperbolic arcsin function and percentile normalized.

- **CyTOF_umap_FlowSOM_clustering**: 
The output of the previous script is used to perform dimensionality reduction using UMAP and FlowSOM clustering in order to assign each cell to a canonical immune subset thereby creating the cellular reference framework. After dissection of leukocytes into broad populations (such as CD4+ T cells, CD8+ T cells, myeloid cells, etc.) each of these subsets were further dissected in an iterative manner using the same computational pipeline.

- **Variance_component_analysis**: 
This script was used to estimate the contribution of genetics, early and late environment of the cell populations present in the reference framework by applying a structured equation model. The input was a table of immune cell population frequencies of the healthy monozygotic and dizygotic twin pairs determined by manual gating.

- **Diffcyt**: 
This script was used to retrieve data-driven immunophenotypes that distinguish MS twins from non-MS twins across all twin pairs and across only untreated twin pairs in order to retrieve immunne perturbations elicited by MS that are not elicited by therapy.

- **InfinityFlow**: 
This script was used to preprocess, cluster and map InfinityFlow data onto the mass cytometry data. Mapping of CITE-seq data was achieved in a similar manner.
