#libraries
library(tidyverse)

metabolite_cidi_bmi <- read.table('/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/norm_metabolite_covariate_revisions.tsv', sep = '\t', header = TRUE)

#running the models 
length = 249
results <- data.frame(metabolite = 1:length, mdd_beta = 1:length, mdd_SE = 1:length, mdd_P = 1:length, n = 1:length)

for(i in 3:251){
  metabolite <- colnames(metabolite_cidi_bmi)[i]
  print(i)
  print(metabolite)
  mod <- lm(metabolite_cidi_bmi[,metabolite] ~ as.factor(Depressed.Ever) + scale(Age) + as.factor(Sex) + as.factor(AS) + scale(BMI) + as.factor(spectrometer) + as.factor(smoking_stat)+ as.factor(ethnicity_collapsed)+ as.factor(uni_nonuni) + scale(townsend_index) , data=metabolite_cidi_bmi)
  results[i-1,1] <- metabolite
  results[i-1,2] <- coef(summary(mod))[2]
  results[i-1,3] <- coef(summary(mod))[39]
  results[i-1,4] <- coef(summary(mod))[113]
  results[i-1,5] <- nobs(mod)
}


results$FDR_P <- p.adjust(results[,"mdd_P"], method = 'BH')
results <- results[order(results$FDR_P),]

write.csv(results ,'/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/cidi_pheno/metabol_mdd_results/adjusted_results_all_covariates_revisions.csv', row.names = F) #change to appropriate file

