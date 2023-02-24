#for documenting how I extracted the spectrometer 
library(stringr)
library(dplyr)
nmr_withqc <- read.table('/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/NMRMetabolomics_withqc.tsv', sep = '\t', header = TRUE)

nmr_baseline_withqc <- select(nmr_withqc, -ends_with("1.0"))
nmr_withqc_noallNA <- nmr_baseline_withqc[rowSums(is.na(nmr_baseline_withqc)) != ncol(nmr_baseline_withqc)-1,]

#edit the spectrometer field to just be numeric 
nmr_withqc_noallNA$spectrometer <- str_sub(nmr_withqc_noallNA$f.23650.0.0, -1,-1)
#merge with the NMR spectrometer field 

spectrometer <- nmr_withqc_noallNA[, c('f.eid', 'spectrometer')]
write.table(spectrometer, '/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/spectrometer_covariate.tsv', sep = '\t', row.names = F, quote = F)
write.table(nmr_withqc_noallNA, '/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/NMR_baseline_withqc.tsv', sep = '\t', row.names = F, quote = F)



#merge with the baseline file currently used in the analysis for the CIDI phenotype 

cidi_pheno <- read.csv('/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/cidi_pheno/norm_metabol_allmerged.csv', header = TRUE)

cidi_spectro_merged <- merge(cidi_pheno, nmr_withqc_noallNA[,c('f.eid', 'spectrometer')], by = 'f.eid')
write.csv(cidi_spectro_merged, '/exports/igmm/eddie/GenScotDepression/users/edavyson/Metabolomics_Aut2021/norm_metabol_allmerged_wispectro.csv', row.names = F)