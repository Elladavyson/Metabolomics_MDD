#libraries
library(tidyverse)

metabolite_cidi <- read.csv('/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/cidi_pheno/norm_metabol_allmerged_wispectro.csv', header = TRUE)


#running the models 
length = 249
in
for(i in 2:250){
  metabolite <- colnames(metabolite_cidi)[i]
  print(i)
  print(metabolite)
  mod <- lm(metabolite_cidi[,metabolite] ~ as.factor(Depressed.Ever) + scale(Age) + as.factor(Sex) + as.factor(AS) + as.factor(spectrometer), data=metabolite_cidi)
  results[i-1,1] <- metabolite
  results[i-1,2] <- coef(summary(mod))[2]
  results[i-1,3] <- coef(summary(mod))[31]
  results[i-1,4] <- coef(summary(mod))[89]
  results[i-1,5] <- nobs(mod)
}


results$FDR_P <- p.adjust(results[,"mdd_P"], method = 'BH')

results <- results[order(results$FDR_P),]

write.csv(results ,'/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/cidi_pheno/metabol_mdd_results/unadjusted_withAS_spectr_results.csv', row.names = F) #change to appropriate file

