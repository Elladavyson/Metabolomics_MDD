#Performing two sample MR

library(remotes)
library(TwoSampleMR)
library(readr)
library(dplyr)
library(ggplot2)
library(data.table)

############## MR LOOP  ##############
#reading in the list of metabolites (three to start with )
  
print("Reading in the metabolite list")
list_all_metabolites <- read.table("/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/all_metabolites_list.txt")

#the linker file and function to get meaningful metabolite_name 
print("Reading in the metabolite linker file ")
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

#Setting up a table for the number of genetic instruments for each test to be saved 
print("Setting up the GI table")

num_metabolites = 249
GI_summary <- data.frame(f.eid = 1:num_metabolites, metabolite_name = 1:num_metabolites, total_SNPs = 1:num_metabolites, total_GI_signif = 1:num_metabolites,  GI_signif_maf_filter= 1:num_metabolites, GI_final_postclumping = 1:num_metabolites)

#The MR loop (to loop through all the metabolites)
print("Beginning the loop")
for (i in 1:dim(list_all_metabolites)[1]) {
  metabolite <- list_all_metabolites[i,]
  print(metabolite)
  file_name <- paste("/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/norm_GWAS/all_part_GWAS/", metabolite, ".gwas/norm_all_gwas.", metabolite, ".glm.linear.ldsc.gz", sep = "")
  print(file_name)
  exposure_data <- read_exposure_data(filename = file_name,
   sep = '\t',
   snp_col = 'SNP', 
   beta_col = 'BETA', 
   se_col ='SE',
   effect_allele_col = 'A1',
   other_allele_col = 'A2', 
   eaf_col = 'A1_FREQ', 
   pval_col = 'P', 
   samplesize_col= 'N')
metabolite_name <- get_metabolitename(metabolite)
print(metabolite_name)
exposure_data$id.exposure <- metabolite
exposure_data$exposure <- metabolite_name

print("Selecting the significant SNPs:")

exposure_data_signif <- exposure_data[exposure_data$pval.exposure < 5e-08/42,]

print("Filtering snps for MAF > 0.005")

exposure_data_signif_maf <- exposure_data_signif %>% filter(eaf.exposure < 0.995 & eaf.exposure > 0.005)

print("Clumping the exposure data:")

exposure_data_clumped <- clump_data(exposure_data_signif_maf)

print("Appending number of genetic instruments to GI table")

GI_summary[i,1] <- metabolite
GI_summary[i,2] <- metabolite_name
GI_summary[i,3] <- dim(exposure_data)[1]
GI_summary[i,4] <- dim(exposure_data_signif)[1]
GI_summary[i,5] <- dim(exposure_data_signif_maf)[1]
GI_summary[i,6] <- dim(exposure_data_clumped)[1]

print(paste("Number of genetic instruments for MR", dim(exposure_data_clumped)[1]))
print("reading in the outcome data")

mdd_outcome_data <- read_outcome_data(filename = '/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/raw_GWAS/mdd_sumstats_noUKB_withN_PGC2_andBETA.gz',
  snps = exposure_data_clumped$SNP,
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
  pos_col = "BP"
)

mdd_outcome_data$outcome <- 'MDD'
mdd_outcome_data$id.outcome <- 'MDD'

print("Harmonising the data")
harmonised_data <- harmonise_data(exposure_dat=exposure_data_clumped, outcome_dat = mdd_outcome_data)

if (dim(harmonised_data )[1] >= 15) {
print("Performing MR")
mr_results <- mr(harmonised_data)
print(head(mr_results))
mr_hetero <- mr_heterogeneity(harmonised_data)
mr_horiz_pleio <- mr_pleiotropy_test(harmonised_data)
res_single <- mr_singlesnp(harmonised_data)
res_loo <- mr_leaveoneout(harmonised_data)

#combining and saving the results 
print("Combining and saving the results")
mr_horiz_pleio <- mr_horiz_pleio %>% slice(rep(1:n(), each = 5)) 
mr_hetero <- rbind(mr_hetero, NA)
mr_hetero <- rbind(mr_hetero, NA)
mr_hetero <- rbind(mr_hetero, NA)
colnames(mr_hetero)[c(1,2,5)] <- c("het_exposure", "het_outcome", "het_method")
colnames(mr_horiz_pleio)[c(1,2,6,7)] <- c("egger_exposure", "egger_outcome", "egger_se", "egger_pval")
combined_results <- cbind(mr_results, mr_hetero[c(5,6,7,8)])
combined_results <- cbind(combined_results, mr_horiz_pleio[c(5,6,7)])

results_directory <- '/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/MR_analysis/norm_GWAS_MR/MR_metabolite_toMDD/MR_results/'

combinedresults_filename <- paste(results_directory, metabolite, ".mr_results.tsv", sep = "")
singlesnp_filename <- paste(results_directory, metabolite, ".mr_singlesnp_results.tsv", sep = "")
leaveoneout_filename <- paste(results_directory,metabolite, ".mr_leaveoneout_results.tsv", sep = "")

print("Writing the results out to EDDIE")

write.table(combined_results, combinedresults_filename, sep = "\t", quote = F, row.names = F)
write.table(res_single, singlesnp_filename, sep = "\t", quote = F, row.names = F)
write.table(res_loo, leaveoneout_filename, sep = "\t", quote = F, row.names = F)

#plots 
print("Making and saving the plots")

plot_directory <- paste(results_directory, "plots/", sep = "")

report_directory <- paste(results_directory, "reports/", metabolite, sep = "")
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

else print("Not enough genetic instruments for MR analysis")

write.table(GI_summary, '/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/MR_analysis/norm_GWAS_MR/MR_metabolite_toMDD/MR_results/met_MDD_GI_summary.tsv', sep = '\t', quote = F, row.names = F)

}