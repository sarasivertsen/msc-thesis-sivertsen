This repository contains R-scripts and datasets used to conduct statistical analyses for the MSc thesis "Exploratory sensory evaluation of added flavor compounds in beer through ranking and free-comment profiling" (NTNU, April 2026).


Contents
  - R-scripts
    1. Ranking analysis (rank_analysis.R): Used to perform Kendall's W, Friedman test, Nemenyi post-hoc test, exact  binomial test
    2. Count data analysis (count_table_analysis.R): Used to perform correspondence analysis, hierarchical clustering analysis, Fisher's exact test
    3. Long-form dataset (long_table_analysis.R): Used to perform PERMANOVA, GLMM
  - Datasets containing preprosessed sensory data from the main beer tasting 
    1. Ranking datasets (rank_v.csv, rank_aa.csv, rank_rk.csv, rank_ia.csv): Contains flavor intensity rankings and metadata
    2. Count dataset (fd_count.csv): Contains aggregated flavor descriptions in a contingency table
    3. Long-form dataset (fd_long.csv): Contains flavor descriptions and metadata 

Requirements: 
R (version ≥ 4.0), and the following packages: DescTools, PMCMRplus, FactoMineR, factoextra, vegan, lme4



