


################ Extracting covariates and metabolite data for lm ################
################ Be on a STAGING node in EDDIE ################

##### phenotype : CIDI depression status 
##### covariates : Age, Sex, AS, BMI, Spectrometer, Ethnic background, Qualifications, Smoking, SES
##### metabolite data : 249 normalised baseline measures from the NMR UKB data 

### libraries 
### make sure the library with all packages in is in the library paths

.libPaths('') 
library(data.table)
library(dplyr)
library(readr)
library(bestNormalize)
library(tidyverse)


### files 
#copied from /exports/igmm/datastore/igmm/datastore/GenScotDepression/data/ukb/phenotypes/fields/2021-04-phenotypes-ukb44797/
recruitment_file <- readRDS("Recruitment.rds")
physicalmeasures <- readRDS("PhysicalMeasures.rds")
baseline_charac <- readRDS("BaselineCharacteristics.rds")
nmr_withqc <- read.table('Metabolomics_Aut2021/NMRMetabolomics_withqc.tsv', sep = '\t', header = TRUE)
touchscreen <- readRDS('Touchscreen.rds')
mdd_cidi_phenotype <- readRDS("MHQ.1907.ukb24262.Process_MH_Questionnaire_Output.rds")

## getting the MDD phenotype and removing the missing values 
mdd_cidi_phenotype <- mdd_cidi_phenotype[,c('f.eid', 'Depressed.Ever')]
mdd_cidi_phenotype_noNA <- na.omit(mdd_cidi_phenotype)

###### GETTING ALL VARIABLES READY #######

### AGE, SEX, AS, BMI ###

## AGE AND AS##
# Recruitment file, age at baseline is the field '21003.0.0. AS is field f.54.0.0'

recruitment_age_AS <- recruitment_file[,c('f.eid', 'f.21003.0.0', 'f.54.0.0')]
colnames(recruitment_age_AS) <- c("f.eid", "Age", "AS")

##or if staging node not available read in using the RDS function
#recruitment_age_AS <- readRDS('Metabolomics_Aut2021/age_AS_recruit.rds')

## SEX ##
#baseline characteristic file - sex is field 31-0.0 

baseline_charac <- readRDS("2021-04-phenotypes-ukb44797/BaselineCharacteristics.rds")
baseline_sex <- baseline_charac[, c("f.eid", "f.31.0.0")]

#recoding the Male and Female strings as 1 and 0 respectively

baseline_sex_coded <- baseline_sex %>% mutate(y = ifelse(f.31.0.0 == "Female", 0, 1))

baseline_sex_coded <- baseline_sex_coded[,c(1,3)] #leaving out the "Male" and "Female" column 
colnames(baseline_sex_coded) <- c("f.eid", "Sex")

##or if staging node not available read in using the RDS function
#baseline_sex_coded <- readRDS('Metabolomics_Aut2021/baseline_sex.rds')

## Townsend deprivation index 
# baseline characteristic file, field is f.189.0.0

baseline_townsend <- baseline_charac[, c("f.eid", "f.189.0.0")]
colnames(baseline_townsend) <- c('f.eid', 'townsend_index')

## BMI ## 

BMI <- physicalmeasures[,c('f.eid', 'f.21001.0.0')]
colnames(BMI)[2] <- 'BMI'

##or if staging node not available read in using the RDS function
#BMI <- readRDS('Metabolomics_Aut2021/BMI_covariate.rds')

##TOUCHSCREEN MEASURES- ethnicity, smoking status and education ## 

touchscreen_variables <- touchscreen[, c('f.eid', 'f.20116.0.0', 'f.6138.0.0', 'f.21000.0.0')]
colnames(touchscreen_variables) <- c('f.eid', 'smoking_stat', 'qualifications', 'ethnicity')

# make a variable for if a participant went to university or not (binary)
# if the participant went to uni was coded as 1, if not, it was coded as 0
touchscreen_variables$uni_nonuni <- ifelse(touchscreen_variables$qualifications == 'College or University degree', 1, 0)

#collapse the ethnicity variable to white/mixed and other 

touchscreen_variables <- touchscreen_variables %>% mutate(ethnicity_collapsed = case_when(
  ethnicity %in% c('White', 'British','Irish','Any other white background') ~ 0, 
  ethnicity %in% c('Asian or Asian British', 'White and Black African', 'Any other mixed background', 'Mixed' , 'Black or Black British' , 'White and Black Caribbean', 'White and Asian') ~ 1,
  ethnicity %in% c('Chinese', 'Pakistani' , 'African' , 'Do not know' , 'Other ethnic group' , 'Indian' , 'Bangladeshi' , 'Caribbean' , 'Any other black background') ~ 2
))

### SPECTROMETER ### 

nmr_baseline_withqc <- select(nmr_withqc, -ends_with("1.0"))
nmr_withqc_noallNA <- nmr_baseline_withqc[rowSums(is.na(nmr_baseline_withqc)) != ncol(nmr_baseline_withqc)-1,]

#edit the spectrometer field to just be numeric 
nmr_withqc_noallNA$spectrometer <- str_sub(nmr_withqc_noallNA$f.23650.0.0, -1,-1)
#merge with the NMR spectrometer field 

spectrometer <- nmr_withqc_noallNA[, c('f.eid', 'spectrometer')]

### METABOLITE DATA ###

raw_metabolite <- nmr_withqc_noallNA[,c(1:250)]

## NORMALISE ##

# replace any 0 entries with an NA #

for (i in colnames(raw_metabolite)[2:250]){
  raw_metabolite[which(raw_metabolite[,i]==0),i] = NA
}

# log transform #
log_metabolite <- raw_metabolite 

for(i in colnames(raw_metabolite)[2:250]){ 
  log_metabolite[,i]<- log(raw_metabolite[,i])
}

# inverse rank normalise transform # 
norm_metabolite <- log_metabolite

for(i in colnames(raw_metabolite)[2:250]){ 
  norm_metabolite[,i]<- orderNorm(log_metabolite[,i])$x.t
}


###### MERGING ALL VARIABLES ######

baseline_variable_list <- list(recruitment_age_AS, baseline_sex_coded, BMI, baseline_townsend, touchscreen_variables)
baseline_df <- baseline_variable_list %>% reduce(inner_join, by = c('f.eid'))

metabolite_variable_list <- list(norm_metabolite,recruitment_age_AS, baseline_sex_coded, BMI, spectrometer, baseline_townsend, touchscreen_variables)
metabolite_df <- metabolite_variable_list %>% reduce(inner_join, by = c('f.eid'))


mdd_variable_list <- list(mdd_cidi_phenotype_noNA, recruitment_age_AS, baseline_sex_coded, BMI, baseline_townsend, touchscreen_variables)
mdd_df <- mdd_variable_list %>% reduce(inner_join, by = 'f.eid')

metabolite_mdd_variable_list <- list(mdd_cidi_phenotype_noNA, norm_metabolite,recruitment_age_AS, baseline_sex_coded, BMI, spectrometer, baseline_townsend, touchscreen_variables)
all_combined_df <- metabolite_mdd_variable_list %>% reduce(inner_join, by = c('f.eid'))

cases_mdd_df <- all_combined_df %>% filter(Depressed.Ever == 1)
controls_mdd_df <- all_combined_df %>% filter(Depressed.Ever == 0)

get_characteristics_from_data <- function(dataset, datasetname, mdd_cases = FALSE) {
  name <- datasetname
  N <- dim(dataset)[1]
  mean_age <- mean(dataset[,'Age'], na.rm = T)
  percentage_female <- (dim(dataset[dataset[,'Sex']==0,])[1]/N)*100
  mean_BMI <- mean(dataset[,'BMI'], na.rm = T)
  mean_townsend <- mean(dataset[, 'townsend_index'], na.rm = T)
  no_NA_smok_status <- as.numeric(table(is.na(dataset[, 'smoking_stat']))[1])
  smok_status <- (length(which(dataset[,'smoking_stat']=='Never'))/no_NA_smok_status)*100
  no_NA_ethnic <- as.numeric(table(is.na(dataset[,'ethnicity_collapsed']))[1])
  ethnic_back <- paste(signif(length(which(dataset[,'ethnicity_collapsed']==0))/no_NA_ethnic,2)*100, ":",signif(length(which(dataset[,'ethnicity_collapsed']==1))/no_NA_ethnic,2)*100, ":", signif(length(which(dataset[,'ethnicity_collapsed']==2))/no_NA_ethnic,2)*100)
  no_NA_uni <- as.numeric(table(is.na(dataset[, 'uni_nonuni']))[1])
  percentage_uni <- (dim(dataset[dataset[,'uni_nonuni']==1,])[1]/no_NA_uni)*100
  if (mdd_cases == TRUE) {
    percent_mdd <- (dim(dataset[dataset[,'Depressed.Ever']==1,])[1]/N)*100
    characteristic_df <- data.frame(Name=name, Number=N, Age = mean_age,female = percentage_female,BMI =mean_BMI, SES= mean_townsend, smoking=smok_status, N_ethnic_back=ethnic_back,uni= percentage_uni, percent_cases = percent_mdd)
  } 
  else{
  characteristic_df <- data.frame(Name=name, Number=N, Age = mean_age,female = percentage_female,BMI =mean_BMI, SES= mean_townsend, smoking=smok_status, N_ethnic_back=ethnic_back,uni= percentage_uni, percent_cases = NA)
  }
  return(characteristic_df)
}

baseline_characteristics <- get_characteristics_from_data(baseline_df, 'Total')
metabolite_characteristics <- get_characteristics_from_data(metabolite_df, 'Baseline Metabolic Data')
mdd_characteristics <- get_characteristics_from_data(mdd_df, 'MDD status', mdd_cases = TRUE)
all_characteristics <- get_characteristics_from_data(all_combined_df, 'MDD status and Metabolic data ', mdd_cases = TRUE)
mdd_cases_characteristics<- get_characteristics_from_data(cases_mdd_df, 'MDD_cases')
mdd_controls_characteristics <- get_characteristics_from_data(controls_mdd_df, 'MDD_controls')

characteristic_table <- rbind(baseline_characteristics, metabolite_characteristics, mdd_characteristics, all_characteristics, mdd_cases_characteristics, mdd_controls_characteristics)
characteristic_table_round <- characteristic_table %>% mutate_if(is.numeric, round, digits = 2)

write.table(all_combined_df, 'Metabolomics_Aut2021/norm_metabolite_covariate_revisions.tsv', sep = '\t', row.names = F, quote = F)
write.table(characteristic_table_round,'Metabolomics_Aut2021/characteristic_paper_table1.tsv', sep = '\t', row.names = F, quote = F) 



##for pippa

write.table(metabolite_df, 'metabolite_UKB.tsv', sep = '\t', row.names = F, quote = F)
write.table(all_combined_df, 'metabolite_MDD_UKB.tsv', sep = '\t', row.names = F, quote = F)




