
#library#

library(data.table)
library(rlang)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(zoo)



#reading in the NMR data table (useful for annotating)

nmr_info <- read.table('Metabolomics/nmr_info.tsv', sep = '\t', header = TRUE)
nmr_info$matchingmet <- paste('f.', nmr_info$UKB.Field.ID, '.0.0', sep = "")

#linker file 

meta_name_linker <- read.csv("Nightingale_markers_desc_full.csv", header = TRUE)
meta_name_linker$linker_ID <- paste(meta_name_linker$Field.ID, ".0.0", sep = '')
meta_name_linker$linker_ID <- paste('f.',meta_name_linker$linker_ID, sep = '')

#### functions to get the metabolite name or the short version of the metabolite name ####

## the full name ##
get_metabolitename <- function(ID_vector) {
  metabolite_names <- c()
  for (i in ID_vector){
    metabolite <- meta_name_linker$Description[meta_name_linker$linker_ID == i]
    metabolite_names <- append(metabolite_names, metabolite)
  }
  return(metabolite_names)
}

## the short name ## 
get_metaboliteshort <- function(ID_vector) {
  metabolite_shorts <- c()
  for(i in ID_vector) {
    metabolite <- nmr_info$Biomarker[nmr_info$matchingmet == i]
    metabolite_shorts <- append(metabolite_shorts, metabolite)
  }
  return(metabolite_shorts)
}

## get the symbol back from the long name ## 

get_metabolitefieldID <- function(ID_vector) {
  metabolite_fieldIDs <- c()
  for(i in ID_vector) {
  metabolite_fieldID <- meta_name_linker$linker_ID[meta_name_linker$Description == i]
  metabolite_fieldIDs <- append(metabolite_fieldIDs, metabolite_fieldID)
  }
  return(metabolite_fieldIDs)
}

#### summarise the findings from an analysis ####
summary_analysis <- function(results, resultstype, beta_colname, p_colname) {
  positive_effect <- results[results[,beta_colname] > 0,]
  negative_effect <- results[results[,beta_colname] < 0,]
  num_neg <- dim(negative_effect)[1]
  num_pos <- dim(positive_effect)[1]
  significant_pos <- positive_effect[positive_effect[,p_colname] == min(positive_effect[,p_colname]),]
  significant_neg <- negative_effect[negative_effect[,p_colname] == min(negative_effect[,p_colname]),]
  beta_pos <- positive_effect[positive_effect[,beta_colname] == max(positive_effect[,beta_colname]),]
  beta_neg <- negative_effect[negative_effect[,beta_colname] == min(negative_effect[,beta_colname]),]
  max_P_value <- max(results[,p_colname])
  min_P_value <- min(results[,p_colname])
  percent_signif <- (sum(results[,p_colname] < 0.05)/length(results[,p_colname])*100)
  sum_signif_results <- sum(results[,p_colname] < 0.05)
  max_beta <- max(results[,beta_colname])
  min_beta <- min(results[,beta_colname])
  summary_df <- data.frame(Analysis= resultstype, Min_beta = min_beta, Max_beta = max_beta, Num_Signif = sum_signif_results,  percent_Signif= percent_signif, Min_FDR_P = min_P_value, Max_FDR_P=max_P_value, Num_pos_Effect= num_pos, No_neg_Effect=num_neg, Most_sig_pos_met = significant_pos$Description ,Most_sig_neg_met= significant_neg$Description, Biggest_beta_pos_met = beta_pos$Description, Biggest_beta_neg_met=beta_neg$Description)
  summary <- c(resultstype, min_beta, max_beta, sum_signif_results, percent_signif, min_P_value, max_P_value, num_pos, num_neg, significant_pos$Description, significant_neg$Description, beta_pos$Description, beta_neg$Description) 
  return(summary_df)
}

#### plotting #### 

#color palette for plot 
palette1_names <- setNames(object=scales::hue_pal()(21), nm = unique(nmr_info$Group))

### plotting MetWAS results ### 

# getting the ID for the metabolite groups for the x axis labels #
get_labelids <- function(results) {
  last_id_positions <- c(0)
  ordered_results <- results[order(results[,'Group']),]  ##order the results by group 
  ordered_results$id <- seq(1, nrow(ordered_results)) ## add an id index 
  for(i in 1:length(unique(ordered_results[,'Group']))) {
    group_id <- unique(ordered_results[,'Group'])[i]  ## select each group 
    group_subset <- ordered_results %>% filter(Group== group_id)
    last_ind <- tail(group_subset$id, n = 1) ##find the last index of each group and append to the vector 
    last_id_positions[i+1] <- last_ind 
    label_id_positions <- rollmean(last_id_positions, k =2) ## now find the middle index for each group, for the plotting of the label on the x axis 
  }
  return(label_id_positions)
}

MetWASplot <- function(results, beta_colname, p_colname) {
  label_ids <- get_labelids(results)
  results <- results[order(results[,'Group']),]
  results[,'id'] <- seq(1, nrow(results)) #adding an identifier 
  annotation_results <- as.data.frame(results %>% filter(!!parse_expr(p_colname) < 0.05) %>% group_by(Group) %>% filter(!!parse_expr(beta_colname) == max(!!parse_expr(beta_colname))))
  plot <- ggplot(results, aes(x = 1:249, y = !!parse_expr(beta_colname), color = as.factor(Group))) + 
    geom_point(alpha = ifelse(results[,p_colname] < 0.05, 1, 0.2)) + 
    ylim(min(results[,beta_colname])- 0.01, max(results[,beta_colname] + 0.01)) +
    scale_x_continuous(breaks = label_ids, labels = unique(results[,'Group']))+
    theme_minimal() + 
    theme(axis.text.x = element_text(face = "bold", size = 5, angle = 45, vjust = 0.6), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black"), legend.position = "bottom") +
    guides(color = guide_legend(title = "Metabolite Group")) + 
    labs(x = 'Metabolite', y = 'Beta') +
    geom_label_repel(annotation_results, mapping = aes(x = id, y = !!parse_expr(beta_colname), label = Biomarker), size = (8*0.352777778), fontface = "bold", show.legend = FALSE, force = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed") 
  return(plot)
}

## function for aligning legend ##

align_legend <- function(p, hjust = 0.5)
{
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
  
  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  
  # there can be multiple guides within one legend box  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    
    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    
    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }
  
  # reconstruct legend and write back
  g$grobs[[legend_index]] <- legend
  g
}


### MR ###

#get the legend for a plot 

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

MR_effectsizes <- function(MR_method, results, col, ymin, ymax, label_size, beta_colname) {
  MR_results_method <- results %>% filter(method == MR_method)
  label_ids <-get_labelids(MR_results_method)
  MR_results_method <- MR_results_method[order(MR_results_method[,'Group']),]
  MR_results_method[,'id'] <- seq(1, nrow(MR_results_method)) #adding an identifier 
  MR_results_method$signif <- ifelse(MR_results_method$FDR_P < 0.05, 'Significant', 'Not Significant')
  
  MR_effect_plot <- ggplot(MR_results_method, aes(x = id, y = !!parse_expr(beta_colname), color = as.factor(Group))) + 
    geom_point(alpha = ifelse(MR_results_method$signif == 'Not Significant', 0.3, 1)) + 
    geom_label_repel(data = subset(MR_results_method, signif == "Significant"), aes(label=Biomarker), size = label_size, show.legend = FALSE) +
    labs(title = MR_method, color = 'Metabolite Group', y = 'Beta', x = 'Metabolite') + geom_hline(yintercept =0, linetype = 'dashed') +
    ylim(ymin, ymax) + 
    scale_x_continuous(breaks = label_ids, labels = unique(MR_results_method[,'Group']))+ 
    theme_minimal() + theme(plot.title= element_text(face = "bold", hjust = 0.5), axis.line = element_line(color = "black"), axis.text.x = element_text(face = "bold", size = 5, angle = 45, vjust = 0.6), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 12))
  return(MR_effect_plot)
}
### MR Odds Ratios ### 

MR_oddsratios <- function(MR_method, results, col, ymin, ymax, label_size, OR_colname) {
  
  #same as the beta graph 
  
  MR_results_method <- results %>% filter(method == MR_method)
  label_ids <-get_labelids(MR_results_method)
  MR_results_method <- MR_results_method[order(MR_results_method[,'Group']),]
  MR_results_method[,'id'] <- seq(1, nrow(MR_results_method)) #adding an identifier 
  MR_results_method$signif <- ifelse(MR_results_method$FDR_P < 0.05, 'Significant', 'Not Significant')

  MR_OR_plot <- ggplot(MR_results_method, aes(x = id, y = !!parse_expr(OR_colname), color = as.factor(Group))) + 
    geom_point(alpha = ifelse(MR_results_method$signif == 'Not Significant', 0.3, 1)) +
    geom_errorbar(aes(ymin = LOW_CI, ymax = HIGH_CI), alpha = ifelse(MR_results_method$signif == 'Not Significant', 0.3, 1))+ 
    geom_label_repel(data = subset(MR_results_method, signif =='Significant'), aes(label = Biomarker), show.legend = F) + 
    geom_hline(yintercept = 1, linetype = 'dashed') + 
    theme_minimal() + theme(plot.title= element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(face = "bold", size = 5, angle = 45, vjust = 0.6), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black"), legend.position = "bottom") + 
    labs(y = 'Odds Ratio (95% CI)', title = MR_method, color = 'Metabolite Group')+ easy_remove_x_axis() + ylim(ymin, ymax)
  return(MR_OR_plot)
}




### colocalisation hypothesis probability plots ### 

coloc_prob_plots <- function(coloc_results_df_chrformat, chraxis_dataframe, hypothesis) {
  
  ##coloc_results_df_chrformat == the coloc results formatted with the chromosome and positions relative (so will look semi like a sparse GWAS plot)
  ##chraxis_dataframe == the axis dataframe which just tells R where to label the chromosome positions 
  ##hypothesis == MUST BE 'PP.H0/1/2/3/4', which probability do you want to plot
  coloc_plot <- ggplot(coloc_results_formatted, aes(x = POScum, y = !!parse_expr(hypothesis)))+geom_point(aes(color = as.factor(Metabolite_short)), size = 1.3, alpha = 0.7) + 
    #custom x axis 
    
    scale_x_continuous(label = coloc_axis_df$CHR, breaks= coloc_axis_df$center)+ geom_hline(yintercept = 0.80, linetype = 'dashed') +
    scale_y_continuous(limits=c(0,1),breaks = seq(0,1.0, by = 0.1)) + 
    theme_minimal() + theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.title= element_text(face = "bold", hjust = 0.5), axis.line.x = element_line(color='black'), axis.line.y = element_line('black')) + 
    labs(x = 'Chromosome', title = hypothesis, color = 'Metabolite')
  
  return(coloc_plot)
}






