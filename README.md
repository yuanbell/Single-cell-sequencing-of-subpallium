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
7. Rscript R/LGE_allweek.r
8. Rscript R/CGE_GW9_GW11.r
9. Rscript R/M1_analysis.r
10. Rscript R/Mapped_analysis.r


# In situ sequencing 
In situ sequencing are performed in three barin sections, and processed data are contained as follow:

1. ISS/GW9_CGE
2. ISS/GW12_MGE_LGE
3. ISS/GW12_CGE

We provide a web page for browse the distribution of single/multiple genes in sections performed by in situ sequencing (http://www.taolabissbrowser.site:3838/).

