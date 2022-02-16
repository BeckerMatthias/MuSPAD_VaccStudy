# MuSPAD_VaccStudy

R code related to the publication: "Comparative magnitude and persistence of humoral SARS-CoV-2 vaccination responses in the adult population in Germany", published in Frontiers in Immunology (https://doi.org/10.3389/fimmu.2022.828053)

The code in the file "main_figures.R" used to generate the graphs used in the paper. In general graphs are exported .svg format, from where they were further processed into the finalized figures using a vector graphics editing software such as InkScape.
For the paper figures Rstudio Version 1.2.5001 was used with R version 3.6.1
Required are the three input files containing measurement and meta data "data_longitudinal_shared.csv", "data_mix_and_match_shared.csv", "data_timepoints_shared.csv".

The code in "supplementary_tables_4_5_6.r" generates results reported in the Supplementary tables 4, 5, and 6. For the paper, the code was run in Rstudio Version 1.4.1717 with R version 4.1.2. Required is the input file "data_mix_and_match_shared.csv" containing measurement and meta data for the mix-and-match cohort.
