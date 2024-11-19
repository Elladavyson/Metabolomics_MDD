# Metabolomic-investigation-of-major-depressive-disorder-identifies-a-potentially-causal-association-w 
[![DOI](https://zenodo.org/badge/606032513.svg)](https://doi.org/10.5281/zenodo.7788315)

Code for analysis presented in the paper at : https://www.biologicalpsychiatryjournal.com/article/S0006-3223(23)00055-0/fulltext

## Analysis Scripts 

### MetWAS

extracting_QC.R - R script to extract spectrometer covariate from metabolite QC data
metabolite_processing_prelm.R - R script extracting the baseline metabolite data, normalising the metabolite data and merging with all covariates of interest 
Unadjusted_MetWAS_lm.R - R script for the linear model looking at associations between each metabolite and MDD with base covariates 
Adjusted_MetWAS_lm.R - R script for the linear model looking at associations between each metabolite and MDD with additional covariates 

### LD-score

LDscore_script - Bash script for running LD-score (ldsc) regression to find the heritability estimates for each metabolite 
rg_mdd_metabolite_script - Bash script for running LD-score (ldsc) regression to estimate the genetic correlation between MDD and metabolites
reading_rg_logfiles_mdd_metabolite - Python script for reading the LD-score regression output files

### MR 

Metabolite_MDD_MR.R - R script running the metabolite -> MDD MR analysis 
MR_loop_MDD_metabolite.R - R script running the MDD -> metabolite MR analysis 


### Colocalisation 

Select_top_SNP_per_region_per_metabolite_ellas_version.R - R script to go through the lead SNPs for each metabolite (in order of significance), and remove all SNPs within a 1Mb region to avoid running multiple of the same colocalisation tests (coloc kb distance = 1Mb)
colocalisation_beginning - R script running colocalisation on metabolites (significant in MR analysis) and MDD


### GTEx MR and colocalisation 

GTEx_parquet_docu - R script documenting the extraction of GTEx data 
GTEx_MR_inprog - R script for eQTLs (identified through GTEx data) -> MDD MR analysis 
colocalisation_GTEx_fullsumstats.txt - R script for colocalisation analysis between GTEx eQTLs and MDD

### Misc scripts (plots and functions)

functions_metabolite.R - R script with common functions used throughout other scripts (getting metabolite names etc)
locus_zoom_meta.R - R script for making locus-zoom plots
