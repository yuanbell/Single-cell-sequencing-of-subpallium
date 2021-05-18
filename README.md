Data and code for the paper 

"Single-cell transcriptomics of the origin and molecular diversity of interneurons in the human fetal brain"

# How to run
There are individual scripts for the different parts of the analysis. Run them in this order:

1.Rscript R/Clustering of allsamples.r

2.Rscript R/LGE_GW9_GW11.r

3.Rscript R/MGE_GW9_GW11.r

4.Rscript R/MGE_allweek.r

5.Rscript R/progenitor_INPs_trajectory.r 

6.Rscript R/progenitor_analysis.r

# R/Clustering of allsamples.r
this script performs the following steps:

1.read in 10X digital expression dataframe and 


