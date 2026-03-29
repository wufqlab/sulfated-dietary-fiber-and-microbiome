### fig 4, use sleuth to do the differential expression analysis 

## load the library
library("sleuth") 
library("files")  # 
library("optparse")
library(biomaRt) 
library(pheatmap)
library(dplyr)
library(tidyverse)
library(ggplot2)

## load the kallisto results into sleuth
print('Loading projects')

setwd("/projects/Prebiotics-Antibiotics_updated/RNA seq/RNA seq/RNAseq-Rockhopper/src/RNA-seq Analysis/")

base_dir <- '../../Data/Clean/kallisto'    # directory where quantified abundances are
sample_id <- dir(file.path(base_dir))
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))

s2c <- read.table("../../Data/Clean/FucKan.RNAseq.DesignMatrix.copy.txt", header = TRUE, stringsAsFactors= FALSE)
s2c <- dplyr::select(s2c, sample=Name, Kanamycin, Fucoidan, group)
s2c <- dplyr::mutate(s2c, path=kal_dirs)
print(s2c)

## 'Fitting model':
so <- sleuth_prep(s2c, ~ Kanamycin + Fucoidan, extra_bootstrap_summary=TRUE)
so <- sleuth_fit(so,~ Kanamycin + Fucoidan + Kanamycin*Fucoidan, fit_name='full') #(full model)
so <- sleuth_fit(so,~ Kanamycin + Fucoidan, fit_name='minimal')  #('reduced')
so <- sleuth_wt(so, which_beta ='Kanamycinb_y', which_model='full') #sleuth_wt

## examine the results of the tests and save the data
# to check the models that have been fit
models(so)

# to examine the results of the test:
results_table <- sleuth_results(so, 'Kanamycinb_y', 'full', test_type='wt')
fname = paste('KanamycinEffects', '.csv', sep='')

# Run the test on beta coefficient on Fucoidan factor
so <- sleuth_wt(so, which_beta ='Fucoidanb_y', which_model='full')
results_table <- sleuth_results(so, 'Fucoidanb_y', 'full', test_type='wt')
fname = paste('FucoidanEffects', '.csv', sep='')

# Run the test on beta coefficient on Kanamycin:Fucoidan factor
so <- sleuth_wt(so, which_beta ='Kanamycinb_y:Fucoidanb_y', which_model='full')
results_table <- sleuth_results(so, 'Kanamycinb_y:Fucoidanb_y', 'full', test_type='wt')
fname = paste('EpistasisEffects', '.csv', sep='')

## go live on the first model built
sleuth_live(so)


## pca plot
plot_pca(so, pc_x = 1L, pc_y = 2L, units = "est_counts", text_labels = FALSE,
         point_size = 2, point_alpha = 0.8, color_by = "group") + theme_light()

# variance explained by PC1, PC2
aa <- plot_pc_variance(so, use_filtered = TRUE, units = "est_counts",
                       pca_number = NULL, scale = FALSE, PC_relative = NULL)
aa
aa$data ### PC1: 68.19%; PC2: 11.97%


## organize the data using the results table from (Kanamycinb_y or 'Kanamycinb_y:Fucoidanb_y') 
# Run the test on beta coefficient on Kanamycin:Fucoidan factor
so <- sleuth_wt(so, which_beta ='Kanamycinb_y:Fucoidanb_y', which_model='full')
results_table <- sleuth_results(so, 'Kanamycinb_y:Fucoidanb_y', 'full', test_type='wt')


so <- sleuth_wt(so, which_beta ='Kanamycinb_y', which_model='full')
results_table <- sleuth_results(so, 'Kanamycinb_y', 'full', test_type='wt')


sig_transcripts <- results_table %>%
  filter(qval < 0.05) # 0.05 0.1
sig_transcripts <- sig_transcripts[order(sig_transcripts$qval), ]

plot_transcript_heatmap(so, units = "tpm",transcripts = sig_transcripts$target_id[1:100], 
                        cluster_transcripts = TRUE, # FALSE
                        trans = "log",color_high = "red",
                        color_mid = "#f7f7f7", color_low = "blue",
                        show_rownames = FALSE)

dev.off()
