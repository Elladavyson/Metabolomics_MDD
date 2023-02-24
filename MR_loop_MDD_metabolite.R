##MDD > metabolite MR script ##

library(remotes)
library(TwoSampleMR)
library(readr)
library(dplyr)
library(ggplot2)

setwd('/exports/igmm/eddie/GenScotDepression/users/edavyson')

#reading in the list of metabolites (three to start with )

list_all_metabolites <- read.table("/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/all_metabolites_list.txt")

print('List of all metabolites read in')

#the linker file 
meta_name_linker <- read.csv('/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/Nightingale_markers_desc_full.csv', header = TRUE)

meta_name_linker$linker_ID <- paste(meta_name_linker$Field.ID, ".0.0", sep = '')
meta_name_linker$linker_ID <- paste('f.',meta_name_linker$linker_ID, sep = '')

get_metabolitename <- function(ID_vector) {
  metabolite_names <- c()
  for (i in ID_vector){
    metabolite <- meta_name_linker$Description[meta_name_linker$linker_ID == i]
    metabolite_names <- append(metabolite_names, metabolite)
  }
  return(metabolite_names)
}
print('linker file read in')

#Setting up a table for the number of genetic instruments for each test to be saved 
print("Setting up the GI table")

num_exposures <- 1
GI_summary <- data.frame(f.eid = 1:1, exposyre_name = 1:1, total_SNPs = 1:1, total_GI_signif = 1:1,  GI_signif_maf_filter= 1:1, GI_final_postclumping = 1:1)

#read in the exposure data (same for every metabolite)
print("Reading in the exposure data ")
exposure_data <- read_exposure_data(filename = '/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/raw_GWAS/mdd_sumstats_noUKB_withN_PGC2_andBETA.gz',
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "FRQ_U_314566",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  ncase_col = "Nca",
  ncontrol_col = "Nco",
  samplesize_col = "N",
  min_pval = 1e-200,
  log_pval = FALSE,
  chr_col = "CHR",
  pos_col = "BP")

exposure_data$id.exposure <- 'MDD'
exposure_data$exposure <- 'MDD'

print('Selecting the significant SNPs')

exposure_data_signif <- exposure_data[exposure_data$pval.exposure < 5e-08,]

print("Filtering snps for MAF > 0.005")

exposure_data_signif_maf <- exposure_data_signif %>% filter(eaf.exposure < 0.995 & eaf.exposure > 0.005)

print("Clumping the significant SNPs")
exposure_data_clumped <- clump_data(exposure_data_signif_maf)

print("Appending number of genetic instruments to GI table")

GI_summary[1,1] <- 'MDD'
GI_summary[1,2] <- 'MDD'
GI_summary[1,3] <- dim(exposure_data)[1]
GI_summary[1,4] <- dim(exposure_data_signif)[1]
GI_summary[1,5] <- dim(exposure_data_signif_maf)[1]
GI_summary[1,6] <- dim(exposure_data_clumped)[1]

print('Beginning of loop with the outcome data')

for (i in 1:dim(list_all_metabolites)[1]) {
  metabolite <- list_all_metabolites[i,]
  print(metabolite)
  file_name <- paste("/exports/eddie/scratch/s2112198/col_num_loop","/norm_all_gwas.", metabolite, ".glm.linear.ldsc", sep = "")
  print(file_name)

  print("Reading in outcome data")

 metabolite_outcome_data <- read_outcome_data(filename = file_name,
   snps = exposure_data_clumped$SNP,
   sep = '\t',
   snp_col = 'SNP', 
   beta_col = 'BETA', 
   se_col ='SE',
   effect_allele_col = 'A1',
   other_allele_col = 'A2', 
   eaf_col = 'A1_FREQ', 
   pval_col = 'P', 
   samplesize_col= 'N')
print('Outcome data read in successfully')

metabolite_name <- get_metabolitename(metabolite)
print(metabolite_name)
metabolite_outcome_data$outcome <- metabolite_name
metabolite_outcome_data$id.outcome <- metabolite

print("Harmonising the data")
harmonised_data <- harmonise_data(exposure_dat=exposure_data_clumped, outcome_dat = metabolite_outcome_data)

print("Performing MR")
mr_results <- mr(harmonised_data)
print(head(mr_results))
mr_hetero <- mr_heterogeneity(harmonised_data)
mr_horiz_pleio <- mr_pleiotropy_test(harmonised_data)
res_single <- mr_singlesnp(harmonised_data)
res_loo <- mr_leaveoneout(harmonised_data)

#combining and saving the results 
print("Combining and saving the results")

output_dir <- '/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/MR_analysis/norm_GWAS_MR/MR_MDD_tometabolite/MR_results/'
mr_horiz_pleio <- mr_horiz_pleio %>% slice(rep(1:n(), each = 5)) 
mr_hetero <- rbind(mr_hetero, NA)
mr_hetero <- rbind(mr_hetero, NA)
mr_hetero <- rbind(mr_hetero, NA)
colnames(mr_hetero)[c(1,2,5)] <- c("het_exposure", "het_outcome", "het_method")
colnames(mr_horiz_pleio)[c(1,2,6,7)] <- c("egger_exposure", "egger_outcome", "egger_se", "egger_pval")
combined_results <- cbind(mr_results, mr_hetero[c(5,6,7,8)])
combined_results <- cbind(combined_results, mr_horiz_pleio[c(5,6,7)])

combinedresults_filename <- paste(output_dir, metabolite, ".mr_results.tsv", sep = "")
singlesnp_filename <- paste(output_dir,  metabolite, ".mr_singlesnp_results.tsv", sep = "")
leaveoneout_filename <- paste(output_dir, metabolite, ".mr_leaveoneout_results.tsv", sep = "")

print("Writing the results out to EDDIE")
write.table(combined_results, combinedresults_filename, sep = "\t", quote = F, row.names = F)
write.table(res_single, singlesnp_filename, sep = "\t", quote = F, row.names = F)
write.table(res_loo, leaveoneout_filename, sep = "\t", quote = F, row.names = F)

#plots 
print("Making and saving the plots")
plot_directory <- paste(output_dir, "plots/", sep  = "")

#create a new directory for each report for the figures? 
report_directory <- paste(output_dir, "reports/", metabolite, sep = "")
dir.create(report_directory)

p1 <- mr_scatter_plot(mr_results, harmonised_data)
ggsave(p1[[1]], file = paste(plot_directory, metabolite, ".mr_scatter_plot.pdf", sep = ""), width = 8, height = 8)

p2 <- mr_forest_plot(res_single)
ggsave(p2[[1]], file = paste(plot_directory, metabolite, ".mr_forest_plot.pdf", sep = ""), width = 8, height = 8)

p3 <- mr_leaveoneout_plot(res_loo)
ggsave(p3[[1]], file = paste(plot_directory,metabolite,".mr_leaveoneout.pdf", sep = ""), width = 8, height = 8)

p4 <- mr_funnel_plot(res_single)
ggsave(p4[[1]], file = paste(plot_directory, metabolite, ".mr_funnel_plot.pdf", sep = ""), width = 8 , height = 8)

mr_report(harmonised_data, output_path = report_directory)

}

write.table(GI_summary, '/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/MR_analysis/norm_GWAS_MR/MR_MDD_tometabolite/MR_results/MDD_tometabolite_GI_summary.tsv', sep = '\t', quote = F, row.names = F)



