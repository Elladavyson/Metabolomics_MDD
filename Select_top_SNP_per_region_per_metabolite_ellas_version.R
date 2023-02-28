## Load requisite libraries 
library(dplyr)
library(zoo)
library(data.table)

## Load in file with independent SNPs and their associated metabolite (here protein for me)
snps = fread("Metabolomics_Aut2021/FUMA_fullsumstats/leadSNPs/stacked_fumasnps_coloc.tsv")
snps <- as.data.frame(snps)
## keep relevant columns-- what is this for?
snps_cols <- snps[,c(1,2,3)]

## rename columns to match Ellas (hey Ella!) 
names(snps_cols)[2] <- "SNP"
names(snps_cols)[1] <- "#CHROM"
names(snps_cols)[3] <- "POS"

## set desired distance threshold (kb) to separate SNPs (here 1Mb at first instance)
kb <- 1e6 #can be changed if needed based on discussions with Andrew, Riccardo etc. 


##### For Rob's file ########################################
## You can ignore these but run if trying to see the example #### 
## make fake metabolite column with two metabolites 
snps_cols$Metabolite = rep(c("Metab.1", "Metab.2"), each = nrow(snps)/2, length.out  = nrow(snps))
## make fake P value column 
snps_cols$P <- runif(nrow(snps), 5e-15, 5e-8)
#######################################################

##########################################################################
###### Loop to choose top SNP per region, and for each metabolite  #######
##########################################################################

## We will create a column called 'keep'
#### keep==TRUE will tell us that it is the top SNP in a region 
#### keep==FALSE will tell us that there is a SNP nearby with a lower P value (so this one goes)
#### we will bring those with keep==TRUE forward to coloc 

## This is only relevant for coloc, the other downstream analyses don't need this approach
## We are just taking the top SNP per region as we'd basically be re-running the same 
## or a very similar coloc analysis using alternative SNPs in the same region (as the 400kb window would move slightly)


## Create column with unique SNP_metabolite pair
## This helps for later as one SNP can be annotated to >1 metabolite (so this prevents it being mixed up)
snps$unique <- paste(snps$SNP, snps$metabolite_name, sep = "_")

## Create column which will tell us whether snp gets kept for coloc 
snps$keep <- NA

##### Now we are ready to run the loop #####

####### STEP ONE ###########
## Subset to a given metabolite (e.g. Metabolite 1)
for(i in unique(snps$metabolite_name)){ 
 ## Subset to Metabolite 
  tmp.met = snps[snps$metabolite_name %in% i,]    
  ## Print  Metabolite 
  print(i)
  
  ######## STEP TWO #########
  ## For that metabolite, we will next loop through each available chromosome  
  for(j in tmp.met$'#CHROM'){
    ## subset to chromosome 
    tmp.chr = tmp.met[tmp.met$'#CHROM' %in% j, ]  
  
    ##### STEP THREE ###########
    
    ## now that we have a unique chromosome + for a unique metabolite
    ## we will find all significant, clumped SNPs on that chromosome 
    ## we will then split the SNPs into different groups (i.e. bins) based on how far apart they are 
    ## the distance criterion relates to the 'kb' variable at the top of the script 

    ## order by POS
    tmp.chr = tmp.chr[order(tmp.chr$POS), ]
    ## Split into different groups that are all separated by a given threshold 
    idx <- c(0, cumsum(abs(diff(tmp.chr$POS)) > kb))
    tmp.chr1 <- dplyr::bind_rows(split(tmp.chr, idx), .id = "Distance.Group")
   
   ##### STEP FOUR ######## 
    ## We have now (i) selected a unique metabolite, (ii) selected a unique chromosome for that metabolite 
    ## and (iii) and zoned in one region on that chromosome (based on a 1Mb window)
   
    ## Next we will find the SNP with the minimum P value (most significant) in that region
    ## this SNP will be kept (set keep=="TRUE")
    ## other SNPs in that region will be discarded (keep=="FALSE")
   
    ## loop through SNPs in a given region  
   for(snp in tmp.chr1$Distance.Group){ 
     tmp.chr2 = tmp.chr1[tmp.chr1$Distance.Group %in% snp, ]
     
     ### if there is just one SNP remaining, we can just keep it (there are no others in region to complete)
      if(nrow(tmp.chr2) <= 1){ 
       snps[which(snps$unique == tmp.chr2$unique),"keep"] <- "TRUE"
       }  else { 
    
     ### if there is more than one, we can only keep the one with the minimum P value 
     ## determine the minimum P value and set keep to true 
     tmp.min <- tmp.chr2[which.min(tmp.chr2$P),]
     tmp.min$keep <- "TRUE"
     ## if there are others, set them to false as they can be discarded 
     tmp.not.min <- tmp.chr2[-which.min(tmp.chr2$P),]
     tmp.not.min$keep <- "FALSE"
     ## combine these SNPs 
     tmp.store = rbind(tmp.min, tmp.not.min)
     ## store the output - we will match it to the other file based on unique SNP+Metabolite combo
     ## make sure order matches 
     snps.tmp = snps[which(snps$unique %in% tmp.store$unique),]
     ids = snps.tmp$unique
     ## match the files 
     tmp.store = tmp.store[match(ids,tmp.store$unique),]
     ## store the keep variable in original dataframe 
     snps[which(snps$unique %in% tmp.store$unique),"keep"] <- tmp.store$keep
       }
     } 
  }
}



## QUALITY CONTROL #####
## Tidy up file 
snps1 = snps[order(snps$metabolite_name, snps$'#CHROM', snps$POS), ]
snps.final = snps1[snps1$keep == "TRUE",]

## Make sure you have all your metabolites 

table(unique(snps.final$metabolite_name) == unique(snps1$metabolite_name)) # i only had 2 fake ones 

## Make sure all original chromosomes are present for all the metabolites 

## Loop through metabolites 
for(orig in unique(snps$metabolite_name)){
  print(orig)
  ## Extract metabolite
  tar = snps[snps$metabolite_name %in% orig,]
  tar1 = snps.final[snps.final$metabolite_name %in% orig,]
  ## Extract chromosomes 
  tar.chrs = unique(tar$'#CHROM')
  tar1.chrs = unique(tar1$'#CHROM')
  ## Are they the same? TRUE will be printed for each metabolite if all is good 
  print(table(tar.chrs == tar1.chrs))
} 


### FINAL STEPS #### 
## if all looks sensible, subset to those being kept for coloc 
#### NOTE: because mine has fake data, there are duplicate SNPs for the same metabolite (and therefore some are not sensible, but please ignore as it is fictitious)


## Save out snps file with keep information (for record keeping + eyeballing)
write.csv(snps1, "filtered_SNPS_annotated_with_top_SNP_info.csv",row.names = F)

## Save out final snps.final file - you can bring this to coloc 
write.csv(snps.final, "filtered_SNPs_for_coloc.csv", row.names = F)