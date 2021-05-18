Data and code for the paper 

"Single-cell transcriptomics of the origin and molecular diversity of interneurons in the human fetal brain"

# Introduction
There are two main parts of analysis of this work, 
one is single cell sequencing, 
another is in situ sequencing.

# Single-cell sequencing
Following several scripts for the different parts of the analysis of single cell sequencing. 
Run them in this order:

1. Rscript R/Clustering of allsamples.r
2. Rscript R/progenitor_INPs_trajectory.r 
3. Rscript R/progenitor_analysis.r
4. Rscript R/MGE_GW9_GW11.r
5. Rscript R/MGE_allweek.r
6. Rscript R/LGE_GW9_GW11.r


# In situ sequencing 
In situ sequencing are performed in three barin sections, and analysis by a tools created by Zhou et al., 2021 (https://www.biorxiv.org/content/10.1101/2020.04.13.038901v1).

