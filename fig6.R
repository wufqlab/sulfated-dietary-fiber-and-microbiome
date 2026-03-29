## for figure 6

library("phyloseq")
library("ggplot2")
library("ggpubr")
library("dunn.test")
library("ape")
library("gridExtra")
library("DESeq2")
library("dplyr")
library("vegan")

#### data and build phyloseq object ####
# directory
setwd("/projects/miseq/")

otutb=as.data.frame(readRDS('seqtab_rm_chimera_all.rds'))  # OTU table
taxtb = read.csv("taxa_species_all.csv", header = TRUE)   # taxonomy table 
sampletb = read.csv("metatable.all.csv", row.names = 1, header = TRUE)    # sample table

SAP = sample_data(sampletb)

otutb_rowname = substring(rownames(otutb), 1, nchar(rownames(otutb))-9) # change row names, deleting the ".fastq.gz"
rownames(otutb) = otutb_rowname

# reformating the tax table
rownames(taxtb)=taxtb$Row.names     # get rid of the row names
taxtb$Row.names=NULL
taxtb_1=as.matrix(taxtb)            # format from data.frame to matrix
OTU = otu_table(otutb[, ], taxa_are_rows = F) # all data
TAX = tax_table(taxtb_1)

# combine the OTU, Tax table, and sample data into phyloseq object
ms = phyloseq(OTU, TAX, SAP) 
random_tree = rtree(ntaxa(ms), rooted = TRUE, tip.label = taxa_names(ms))
ms = merge_phyloseq(ms, random_tree)

# Show available ranks in the dataset, and remove NA at phylum level
rank_names(ms)
table(tax_table(ms)[, "Phylum"], exclude = NULL)
# ensures that features with ambiguous phylum annotation are also removed.
ms <- subset_taxa(ms, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# to remove the tata (column sum or row sum) samller than 2 from the original OTU table and Taxo table
ms.pt = prune_taxa(taxa_sums(ms) > 1, ms) # some using taxa_sums(ms) >0, maybe using >1
ms.pts = prune_samples(sample_sums(ms.pt)>=1000, ms.pt)   # Remove samples with less than 1000 reads (sum)
ms.pts.nor = transform_sample_counts(ms.pts, function(otu) otu/sum(otu))
# ms.pts.nor.log = transform_sample_counts(ms.pts.nor, function(x) log(1 + x)) # log: not improve much

#### end ####


#### plot NMDS for time points
ms.Aug29 = subset_samples(ms.pts.nor, Timepoint == "Aug.29") # t0, right before Treatment
ms.Aug31 = subset_samples(ms.pts.nor, Timepoint == "Aug.31")
ms.Sep2 = subset_samples(ms.pts.nor, Timepoint == "Sep.02")
ms.Sep5 = subset_samples(ms.pts.nor, Timepoint == "Sep.05")
ms.Sep8 = subset_samples(ms.pts.nor, Timepoint == "Sep.08")
ms.Sep10 = subset_samples(ms.pts.nor, Timepoint == "Sep.10")
ms.Sep14 = subset_samples(ms.pts.nor, Timepoint == "Sep.14")
ms.Sep18 = subset_samples(ms.pts.nor, Timepoint == "Sep.18")
ms.Sep24 = subset_samples(ms.pts.nor, Timepoint == "Sep.24")
ms.Sep30 = subset_samples(ms.pts.nor, Timepoint == "Sep.30")

ord.nmds.bray = ordinate(ms.treat.fk , method = "NMDS", distance = "bray") 
plot_ordination(ms.treat.fk, ord.nmds.bray, color="Treatment", title="Bray MDS") +
  theme_light() + theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank())

metadf <- as(sample_data(ms.treat.fk), "data.frame")
adonis2(phyloseq::distance(ms.treat.fk, method = "bray", perm=999) ~ Treatment,
       data = metadf)
 
#### plot relative abundance with time (baseline, treatment, recovery)
ms.period = merge_phyloseq(ms.Aug29, ms.Sep14, ms.Sep30) 
ord.nmds.bray = ordinate(ms.period  , method = "PCoA", distance = "bray") # PCoA, NDMS
plot_ordination(ms.period , ord.nmds.bray, color="Timepoint", title="Bray NMDS") +
  theme_light() + theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank())


# plot the top 100 taxa
top100otu <- names(sort(taxa_sums(ms.period), TRUE)[1:100])  
dat.dataframe = psmelt(ms.period)
dat.dataframe$Genus.x[which(!(dat.dataframe$OTU %in% top100otu))]=NA 

getPalette = colorRampPalette(brewer.pal(8, "Paired")) 
PhylaPalette = getPalette(7)

ggplot(dat.dataframe, aes(x=Sample, y=Abundance, fill=Phylum)) +
      geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90))+ 
      scale_fill_manual(values = PhylaPalette)


 ############################################################################## 
 ##### to identify significant taxas after Kan treatmnet, and then test / compare their recovery in FK group vs K groups
 # 1) filter microbes with abundance decreased after Kan; 2) filter microbes that can not recover to original levels (~0) even after Kan stopped
 
 #### Upload the deseq2 results bw CTL vs K, ctl vs FK
 ctl.k.deseq  <- read.csv("Aug31Sep02051014/OTUs_CtlvsK_Aug29toSep30_pval=1.csv", header = TRUE)  
 ctl.fk.deseq <- read.csv("Deseq2/Aug31Sep02051014/OTUs_CtlvsFK_Aug29toSep30_pval=1.csv", header = TRUE) 
 
 sig.ctl.k = read.csv("Deseq2/Aug31Sep02051014/sigOTUs_CtlvsK_Aug31Sep02.05.10.14.csv", header = TRUE) # padj <= 0.01
 sig.ctl.k.otu <- sig.ctl.k %>% select(1)   
 colnames(sig.ctl.k.otu)[1] <- "otu"
 
 sig.ctl.fk = read.csv("Deseq2/Aug31Sep02051014/sigOTUs_CtlvsFK_Aug31Sep02.05.10.14.csv", header = TRUE)    # padj <= 0.01
 sig.ctl.fk.otu <- sig.ctl.fk %>% select(1)     
 colnames(sig.ctl.fk.otu)[1] <- "otu"
 
 sig.ctl.k.deseq <- merge(ctl.k.deseq, sig.ctl.k.otu, "otu")
 sig.ctl.fk.deseq <- merge(ctl.fk.deseq, sig.ctl.fk.otu, "otu")
 
 sig.ctl.k.deseq.dec.fc0 <- sig.ctl.k.deseq[which((sig.ctl.k.deseq$date=="Aug31" | 
                                               sig.ctl.k.deseq$date=="Sep02" ) &
                                               sig.ctl.k.deseq$log2FoldChange <= 0), ]
 
 ## 1: setup the threshold, using the size of the interquartile range 
 sig.ctl.k.deseq.dec.t0 <- sig.ctl.k.deseq[which(sig.ctl.k.deseq$date=="Aug29"), ] 
 sum.t0.fc <- summary(sig.ctl.k.deseq.dec.t0$log2FoldChange)
 thres <- abs(sum.t0.fc[2]) + abs(sum.t0.fc[5])
 thres
 
 sig.ctl.k.deseq.rec <- sig.ctl.k.deseq[which((sig.ctl.k.deseq$date=="Sep24" | 
                                                 sig.ctl.k.deseq$date=="Sep30") &
                                                 sig.ctl.k.deseq$log2FoldChange > -(thres)), ]
 sig.ctl.k.deseq.dec <-  sig.ctl.k.deseq.dec.fc0[which( !sig.ctl.k.deseq.dec.fc0$otu %in% sig.ctl.k.deseq.rec$otu), ]
 
 # seperate the data into t0, t1, t2, t3, t4 for Ctl vs K
 sig.ctl.k.deseq.dec.t1 <- sig.ctl.k.deseq.dec
 sig.ctl.k.deseq.dec.t1$period <- "t1"  

 sig.ctl.k.deseq.dec.tm <- sig.ctl.k.deseq[which((sig.ctl.k.deseq$date=="Sep14" | sig.ctl.k.deseq$date=="Sep10")), ]  
 sig.ctl.k.deseq.dec.tm2 <- sig.ctl.k.deseq.dec.tm[which(sig.ctl.k.deseq.dec.tm$otu %in% unique(sig.ctl.k.deseq.dec$otu)), ]
 sig.ctl.k.deseq.dec.tm2$period <- "t2" 
 
 sig.ctl.k.deseq.dec.te  <- sig.ctl.k.deseq[which(( sig.ctl.k.deseq$date=="Sep24" | sig.ctl.k.deseq$date=="Sep30")), ]   
 sig.ctl.k.deseq.dec.te2 <- sig.ctl.k.deseq.dec.te[which(sig.ctl.k.deseq.dec.te$otu %in% unique(sig.ctl.k.deseq.dec$otu)), ]
 sig.ctl.k.deseq.dec.te2$period <- "t3"   
 
 sig.ctl.k.deseq.dec.t02 <- sig.ctl.k.deseq.dec.t0[which(sig.ctl.k.deseq.dec.t0$otu %in% unique(sig.ctl.k.deseq.dec$otu)), ]
 sig.ctl.k.deseq.dec.t02$period <- "t0" 
 sig.ctl.k.deseq.dec.all <- rbind(sig.ctl.k.deseq.dec.t1, sig.ctl.k.deseq.dec.tm2, sig.ctl.k.deseq.dec.te2,sig.ctl.k.deseq.dec.t02)  
 
 ggplot(sig.ctl.k.deseq.dec.all, aes(x=period, y=log2FoldChange, color = period)) +
   geom_boxplot(outlier.shape = NA, lwd=0.3, alpha=1)  + theme_light() + labs(title = "control vs Kan-fill0s") +
   geom_point(size = 1,position = position_jitterdodge(0.2), alpha = 0.6) #+ ylim(-12, 10)


# seperate the data into t0, t1, t2, t3, t4 for Ctl vs FK
 sig.ctl.fk.deseq.dec <- sig.ctl.fk.deseq[which((sig.ctl.fk.deseq$date=="Aug31" | sig.ctl.fk.deseq$date=="Sep02" ) ), ]   
 sig.ctl.fk.deseq.dec.t1 <- sig.ctl.fk.deseq.dec[which(sig.ctl.fk.deseq.dec$otu %in% unique(sig.ctl.k.deseq.dec$otu)), ]  
 sig.ctl.fk.deseq.dec.t1$period <- "t1"  
 
 sig.ctl.fk.deseq.dec.tm <- sig.ctl.fk.deseq[which(( sig.ctl.fk.deseq$date=="Sep14" | sig.ctl.fk.deseq$date=="Sep10" )), ]   
 sig.ctl.fk.deseq.dec.tm2 <- sig.ctl.fk.deseq.dec.tm[which(sig.ctl.fk.deseq.dec.tm$otu %in% unique(sig.ctl.k.deseq.dec$otu)), ]  
 sig.ctl.fk.deseq.dec.tm2$period <- "t2" 
 
 sig.ctl.fk.deseq.dec.te  <- sig.ctl.fk.deseq[which(( sig.ctl.fk.deseq$date=="Sep24"  |  sig.ctl.fk.deseq$date=="Sep30")), ]
 sig.ctl.fk.deseq.dec.te2 <- sig.ctl.fk.deseq.dec.te[which(sig.ctl.fk.deseq.dec.te$otu %in% unique(sig.ctl.k.deseq.dec$otu)), ]
 sig.ctl.fk.deseq.dec.te2$period <- "t3"  
 
 sig.ctl.fk.deseq.dec.t0 <- sig.ctl.fk.deseq[which(sig.ctl.fk.deseq$date=="Aug29"), ] 
 sig.ctl.fk.deseq.dec.t02 <- sig.ctl.fk.deseq.dec.t0[which(sig.ctl.fk.deseq.dec.t0$otu %in% unique(sig.ctl.k.deseq.dec$otu)), ]
 sig.ctl.fk.deseq.dec.t02$period <- "t0"  
 
 sig.ctl.fk.deseq.dec.all <- rbind(sig.ctl.fk.deseq.dec.t1, sig.ctl.fk.deseq.dec.tm2, sig.ctl.fk.deseq.dec.te2, sig.ctl.fk.deseq.dec.t02)  

 ggplot(sig.ctl.fk.deseq.dec.all, aes(x=period, y=log2FoldChange, color = period)) +
   geom_boxplot(outlier.shape = NA, lwd=0.3, alpha=1)  + theme_light() +labs(title = "control vs Fucoidan/Kan-filled0s") +
   geom_point(size = 1,position = position_jitterdodge(0.2), alpha = 0.6)  
 
 sig.ctl.k.deseq.dec.all$comparison <- "Ctl vs K"
 sig.ctl.fk.deseq.dec.all$comparison <- "Ctl vs FK"
 sig.ctl.k.fk.deseq.dec.all.comb <- rbind(sig.ctl.k.deseq.dec.all, sig.ctl.fk.deseq.dec.all)
 ggplot(sig.ctl.k.fk.deseq.dec.all.comb, aes(x=period, y=log2FoldChange, color = comparison)) +
   geom_boxplot(outlier.shape = NA, lwd=0.3, alpha=1)  + theme_light() +labs(title = "Ctl-K vs Ctl-FK (sigCtl-K)")+
   geom_point(size = 1.2,position = position_jitterdodge(0.2), alpha = 0.5)  +   ylim(-11.5, 4.5) + 
   theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank()) # +
 
