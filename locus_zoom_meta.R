###LD structure around the colocalised SNP ###

library(data.table)
library(plyr)
library(dplyr)
library(readr)
library(ggplot2)

source('/Users/ellad/UniversityEdinburgh/PhD/CodePhD/Metabolomics/functions_metabolite.R')
#scale function for the x axis 
scaleFUN <- function(x) {
  tmp <- x/1000000
  return(tmp)
}

#LD colours 
# 1.0--0.8 - #FF0000
# 0.8 --0.6 - #FFA500
# 0.6--0.4- #00FF00
# 0.4--0.2 -#87CEFA
#0.2--0.0 = #000080
#NA = #7F7F7F
#top hit = #7D26CD
setwd('/Users/ellad')
LD.colours <- data.frame(LD = as.character(seq(from=0, to=1, by = 0.1)), Colour = c('#000080', rep(c('#000080', '#87CEFA', '#00FF00', '#FFA500', '#FF0000'), each = 2)), stringsAsFactors =  FALSE)

locus_plot <- function(sumstats, LD_file, lowerlim, upperlim, colocalised_SNP, field_id) {
  colnames(sumstats)[1:2] <- c('CHR', 'BP')
  subset_gwas <- sumstats[which(sumstats$BP >= lowerlim),]
  subset_gwas <- subset_gwas[which(subset_gwas$BP <= upperlim),]
  SNPs_in_LD_matrix_df <- merge(subset_gwas, LD_file[,c('SNP_B', 'R2')], by.x = 'SNP', by.y = 'SNP_B')
  #manually change the 0 P values to 1e-320
  SNPs_in_LD_matrix_df[which(SNPs_in_LD_matrix_df$P < 1e-320), 'P'] <- 1e-320
  coloc_snp <- SNPs_in_LD_matrix_df %>% filter(SNP == colocalised_SNP)
  plot <- ggplot(SNPs_in_LD_matrix_df, aes(x = BP, y = -log10(P))) +
    geom_point(aes(colour = cut(R2, c(0.0,0.2, 0.4, 0.6, 0.8, 1.0)), fill = cut(R2, c(0.0,0.2, 0.4, 0.6, 0.8, 1.0))), shape = 21, size = 4, color = 'black') +
    scale_color_manual(name = 'R2', values = c('#000080','#87CEFA', '#00FF00', '#FFA500',  '#FF0000'), na.translate = F)+ scale_fill_manual(name = expression(R^{2}), values = c('#000080','#87CEFA', '#00FF00', '#FFA500',  '#FF0000'), na.translate = F, labels = c('0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0'))+
    theme_minimal() + theme(plot.title= element_text(face = "bold", hjust = 0.5), axis.line = element_line(color = "black"), text = element_text(size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
    labs(x= 'Chromosome 11 (Mb)', y = '-log10(P)', fill = expression(R^{2}), title = get_metabolitename(field_id)) + 
    scale_x_continuous(labels =scaleFUN) +ylim(0,350)+geom_text(x = 61700000, y = 340, label = colocalised_SNP)
  return(plot)
}

#label argument I have taken out for now - perhaps put back in? 
#+geom_label(data = SNPs_in_LD_matrix_df %>% filter(SNP == colocalised_SNP), aes(label = SNP))

##LA to TFA plot ## 

la_sumstats_chr11 <- read.table('/Users/ellad/UniversityEdinburgh/PhD/Data/UKB/Metabolomics/chr11_signif_metabolites/f.23456.0.0chr_11.tsv', sep = '\t', header = T)
rs174578_ld <- read.table('/Users/ellad/UniversityEdinburgh/PhD/Data/UKB/Metabolomics/SNP_ld/rs174578.ld', header = T)

LA_TFA_locus_plot <- locus_plot(la_sumstats_chr11, rs174578_ld, 61450000, 61800000, 'rs174578', 'f.23456.0.0')
## DHA plot ##

DHA_sumstats_chr11 <- read.table('/Users/ellad/UniversityEdinburgh/PhD/Data/UKB/Metabolomics/chr11_signif_metabolites/f.23450.0.0chr_11.tsv', sep = '\t', header = T)
rs2727271_ld <- read.table('/Users/ellad/UniversityEdinburgh/PhD/Data/UKB/Metabolomics/SNP_ld/rs2727271.ld', header = T)

DHA_locus_plot <- locus_plot(DHA_sumstats_chr11, rs2727271_ld, 61450000, 61800000, 'rs2727271', 'f.23450.0.0')
DHA_locus_plot_sameLD <- locus_plot(DHA_sumstats_chr11, rs174578_ld, 61450000, 61800000, 'rs174578', 'f.23450.0.0')

##DHA to Total Fatty Acids ## 

DHA_TFA_sumstats_chr11 <- read.table('/Users/ellad/UniversityEdinburgh/PhD/Data/UKB/Metabolomics/chr11_signif_metabolites/f.23457.0.0chr_11.tsv', sep = '\t', header = T)
rs174575_ld <- read.table('/Users/ellad/UniversityEdinburgh/PhD/Data/UKB/Metabolomics/SNP_ld/rs174575.ld', header = T)

DHA_TFA_locus_plot <- locus_plot(DHA_TFA_sumstats_chr11, rs174575_ld, 61450000, 61800000, 'rs174575', 'f.23457.0.0')
  
## MDD locus plot - slightly different scales etc ##

mdd_chr11 <- read.table('/Users/ellad/UniversityEdinburgh/PhD/Data/UKB/Metabolomics/mdd_sumstats_chr11_noUKB_withNandBETA.tsv', sep = '\t', header = T)
subset_mdd_chr11 <- mdd_chr11[which(mdd_chr11$BP >= 61450000),]
subset_mdd_chr11 <- subset_mdd_chr11[which(subset_mdd_chr11$BP <= 61800000),]
mdd_snps_in_LD_matrix_df <- merge(subset_mdd_chr11, rs174578_ld[,c('SNP_B', 'R2')], by.x = 'SNP', by.y = 'SNP_B')

mdd_locus_plot <- ggplot(mdd_snps_in_LD_matrix_df, aes(x = BP, y = -log10(P))) +
  geom_point(aes(colour = cut(R2, c(0.0,0.2, 0.4, 0.6, 0.8, 1.0)), fill = cut(R2, c(0.0,0.2, 0.4, 0.6, 0.8, 1.0))), shape = 21, size = 4, color = 'black') +
  scale_color_manual(name = 'R2', values = c('#000080','#87CEFA', '#00FF00', '#FFA500',  '#FF0000'), na.translate = F)+ scale_fill_manual(name = expression(R^{2}), values = c('#000080','#87CEFA', '#00FF00', '#FFA500',  '#FF0000'), na.translate = F, labels = c('0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0'))+
  theme_minimal() + theme(plot.title= element_text(face = "bold", hjust = 0.5), axis.line = element_line(color = "black"), text = element_text(size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  labs(x= 'Chromosome 11 (Mb)', y = '-log10(P)', fill = expression(R^{2}), title = 'MDD') + 
  scale_x_continuous(labels =scaleFUN)+geom_text(x = 61700000, y = 5, label = 'rs174578')+ylim(0,8)

## the gene tracks of the region ##

#### not sure what the above is , but doing it with ggplot ####
more_UCSC_genes <- read.table('/Users/ellad/UniversityEdinburgh/PhD/Data/UKB/Metabolomics/6145-618_UCSC_genes', sep = '\t', header = F)
colnames(more_UCSC_genes) <- c('Chrom', 'Start', 'End', 'Gene', 'Description')

more_UCSC_genes$Colour <- 'black' #change this depending on the genes which we tested 
genes_in_GTEx <- c('FEN1', 'MYRF', 'FADS1', 'FADS2', 'TMEM258', 'FADS3')

more_UCSC_genes$Colour <- ifelse(more_UCSC_genes$Gene %in% genes_in_GTEx, 1, 0)
#setting the Y variable for stagering the genes 

if(length(more_UCSC_genes[, "Gene"]) > 15) {
  y = rep(c(5,4,3,2,1), times = length(more_UCSC_genes[,'Gene']))
  more_UCSC_genes$Y = y[1:length(more_UCSC_genes$Gene)]
} else {
  y = rep(c(1.2, 1.1, 1.0, 0.9, 0.8), times = length(more_UCSC_genes[,"Gene"]))
  more_UCSC_genes$Y = y[1:length(more_UCSC_genes$Gene)]
}

gene_track_plot <- ggplot(more_UCSC_genes, aes(x = BP, y = Y)) + geom_segment(aes(x = Start, xend = End, y = Y, yend = Y, color = as.factor(Colour)))+
  geom_text(data = more_UCSC_genes, aes(x= (Start+End)/2, y = Y+0.02, label = Gene, fontface = 'italic', color = as.factor(Colour))) + 
  theme_minimal() + theme(plot.title= element_text(face = "bold", hjust = 0.5), axis.line = element_line(color = "black"), text = element_text(size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  labs(y = "", x = 'Chromosome 11 (Mb)') +scale_color_manual(values = c('black', 'red')) + scale_x_continuous(labels =scaleFUN)

gene_track_plot <- gene_track_plot+guides(color="none")                                                                                                                                                                                                                                                                                                                                                                     

all_plots <- ggarrange(LA_TFA_locus_plot,DHA_locus_plot, DHA_TFA_locus_plot, mdd_locus_plot, gene_track_plot, common.legend = T, nrow = 5, ncol = 1)
ggsave('/Users/ellad/UniversityEdinburgh/PhD/Data/UKB/Metabolomics/plots/locus_zoom_plots.tiff', device = 'tiff', all_plots, scale = 1, width = 184, height = 285, units ="mm")

#DHA_plots <- ggarrange(DHA_locus_plot, DHA_TFA_locus_plot, mdd_locus_plot, gene_track_plot, common.legend = T, labels = c('A', 'B', 'C', 'D'), nrow = 4, ncol =1)
#ggsave('Metabolomics/plots/DHA_all_locus_zoom.tiff', device = 'tiff', DHA_plots, scale = 1, width = 190, height = 243.2, units ="mm", bg = "white")
DHA_plot_single <- ggarrange(mdd_locus_plot,DHA_locus_plot, gene_track_plot, common.legend = T, labels = c('A', 'B', 'C'), nrow = 3, ncol = 1)
#la_tfa_plots <- ggarrange(LA_TFA_locus_plot, mdd_locus_plot, gene_track_plot, common.legend = T, labels = c('A', 'B', 'C'), nrow = 3, ncol = 1)
ggsave('/Users/ellad/UniversityEdinburgh/PhD/Year_1/Metabolomics_Submission/Acceptance/Proofd/Figure_3.tiff', device = 'tiff', DHA_plot_single, scale = 1, width = 190, height = 243.2, units ="mm", bg = "white")
