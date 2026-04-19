This repository contains R-scripts and datasets used to conduct statistical analyses for the MSc thesis "Exploratory sensory evaluation of added flavor compounds in beer through ranking and free-comment profiling" (NTNU, April 2026).

The thesis aimed to evaluate the perception of four different flavor compounds (vanillin, acetic acid, raspberry ketone, and isoamyl acetate) in two different base beers (a porter and an amber ale). Two flavor compounds were evaluated in the porter (vanillin and acetic acid), and two flavor compounds were evaluated in the amber ale (raspberry ketone and isoamyl acetate). An untrained panel consisting of 31 participants recruited from a university setting performed a ranking test to evaluate the flavor intensity, and the free-comment method was used to collect flavor descriptions. The resulting sensory data have been used to generate the six following datasets.  

<ins>Contents</ins>
 - **Datasets containing preprosessed sensory data from the main beer tasting** 
    1. Ranking datasets (_rank_v.csv, rank_aa.csv, rank_rk.csv, rank_ia.csv_): Contains flavor intensity rankings and metadata
    2. Count dataset (_fd_count.csv_): Contains aggregated flavor descriptions in a contingency table
    3. Long-form dataset (_fd_long.csv_): Contains flavor descriptions and metadata 

  - **R-scripts**
    1. Ranking analysis (_rank_analysis.R_): Used to perform Kendall's W, Friedman test, Nemenyi post-hoc test, exact binomial test
    3. Count data analysis (_count_table_analysis.R_): Used to perform correspondence analysis, hierarchical clustering analysis, Fisher's exact test
    5. Long-form dataset (_long_table_analysis.R_): Used to perform PERMANOVA, GLMM
 
Requirements: 
R (version ≥ 4.0), and the following packages: DescTools, PMCMRplus, FactoMineR, factoextra, vegan, lme4.



