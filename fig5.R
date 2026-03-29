# fig 5, ex vivo 

library("phyloseq")
library("ggplot2")
library("ggpubr")
library("dunn.test")
library("vegan")
library("devtools")
library("pairwiseAdonis")
library("ape")
########################################################################################
# data read and preliminary Group
setwd("/project/exvivo/Miseq data/")
otutb=as.data.frame(readRDS('seqtab_rm_chimera.rds'))  # OTU table
taxtb = read.csv("taxa_species.csv", header = TRUE)   # taxonomy table 
sampletb = read.csv("metadata.csv", row.names = 1, header = TRUE)    # sample table

otutb_rowname = substring(rownames(otutb), 1, nchar(rownames(otutb))-9)  
rownames(otutb) = otutb_rowname

# reformating the tax table
rownames(taxtb)=taxtb$Row.names     
taxtb$Row.names=NULL
taxtb_1=as.matrix(taxtb)             

OTU = otu_table(otutb[97:288, ], taxa_are_rows = F) # 1-96 is data of Plate 3
TAX = tax_table(taxtb_1)
SAP = sample_data(sampletb)
exvivo = phyloseq(OTU, TAX, SAP) 

random_tree = rtree(ntaxa(exvivo), rooted = TRUE, tip.label = taxa_names(exvivo))
exvivo.t = merge_phyloseq(exvivo, random_tree)

# to remove the tata (column sum or row sum) samller than 2 from the original OTU table and Taxo table
exvivo.pt = prune_taxa(taxa_sums(exvivo.t) > 2, exvivo.t)
# Remove samples with less than 1000 reads (sum)
exvivo.pts = prune_samples(sample_sums(exvivo.pt)>=1000, exvivo.pt)  

exvivo.pts.nor = transform_sample_counts(exvivo.pts, function(otu) otu/sum(otu))

#####################################################################################
######################################################################## host 1, ampicillin
ps.nor.P = subset_samples(exvivo.pts.nor, Type == "Pt")
ps.nor.P.24h = subset_samples(ps.nor.P, Time == "24h")

p00 = plot_bar(ps.nor.P.24h, fill = "Family")
plot_bar(ps.nor.P.24h, fill = "Genus.x")
p00 <- p00 + facet_wrap(~Group)  
p00
ggsave("ps.nor.P.24h-family.pdf", p00, device = "pdf", dpi = 480, width = 8.27, height = 5.87)

# select the desired groups
ps.nor.P.24h.Amp.1 = subset_samples(ps.nor.P.24h, Group == "Ctl" )
ps.nor.P.24h.Amp.2 = subset_samples(ps.nor.P.24h, Group == "Fcdn")
ps.nor.P.24h.Amp.3 = subset_samples(ps.nor.P.24h, Group == "H.Amp.Fcdn")
ps.nor.P.24h.Amp.4 = subset_samples(ps.nor.P.24h, Group == "H.Amp")
ps.nor.P.24h.Amp = merge_phyloseq(ps.nor.P.24h.Amp.1, ps.nor.P.24h.Amp.2, ps.nor.P.24h.Amp.3, ps.nor.P.24h.Amp.4)

top30otu <- names(sort(taxa_sums(ps.nor.P.24h.Amp), TRUE)[1:30])

dat.dataframe = psmelt(ps.nor.P.24h.Amp)
dat.dataframe$Genus.x[which(!(dat.dataframe$OTU %in% top30otu))]=NA 

p <- ggplot(dat.dataframe, aes(x=Sample, y=Abundance, fill=Genus.x)) + 
  geom_bar(stat="identity") +  theme(axis.text.x = element_text(angle = 90))

# plot and reduce the legend size (key size and text size)
p + guides(fill=guide_legend(ncol = 1,byrow = FALSE, 
                             keywidth =unit(0.4,"cm"),
                             keyheight=unit(0.4,"cm")))+
  theme(legend.text = element_text(size=9))

# ndms plot
ord.nmds.bray = ordinate(ps.nor.P.24h.Amp, method = "NMDS", distance = "bray")
p24A <- plot_ordination(ps.nor.P.24h.Amp, ord.nmds.bray, color="Group", title="Bray NMDS") + theme_bw() 
p24A$layers <- p24A$layers[-1]
p24A <- p24A + geom_point( shape = 1, size = 3, stroke = 1.2)  
p24A
# weighted unifrac plot
ordu = ordinate(ps.nor.P.24h.Amp, "PCoA", "unifrac", weighted=TRUE)
p24A.f <- plot_ordination(ps.nor.P.24h.Amp, ordu, color="Group", title="Weighted Unifrac")+ theme_bw() 
p24A.f$layers <- p24A.f$layers[-1]
p24A.f <- p24A.f + geom_point( shape = 1, size = 3, stroke = 1.2)  
p24A.f

# PERMANOVA
metada <- as(sample_data(ps.nor.P.24h.Amp), "data.frame")
adonis(distance(ps.nor.P.24h.Amp, method = "bray", perm=999) ~ Group,
       data = metada)

#### Pairwise PERMANOVA
z <- pairwise.adonis(distance(ps.nor.P.24h.Amp, method = "bray", perm=999), metada$Group, 
                sim.method = "bray", p.adjust.m = "bonferroni")
write.csv(z, 'pairwise.adonis.P.24h.Amp.csv')

#####################################################################################
######################################################################## host 1, kanamycin
ps.nor.P = subset_samples(exvivo.pts.nor, Type == "Pt")
ps.nor.P.24h = subset_samples(ps.nor.P, Time == "24h")

# select the desired groups
ps.nor.P.24h.Kan.1 = subset_samples(ps.nor.P.24h, Group == "Ctl2" )
ps.nor.P.24h.Kan.2 = subset_samples(ps.nor.P.24h, Group == "Fcdn")
ps.nor.P.24h.Kan.3 = subset_samples(ps.nor.P.24h, Group == "H.Kan.L.Fcdn")
ps.nor.P.24h.Kan.4 = subset_samples(ps.nor.P.24h, Group == "H.Kan")
ps.nor.P.24h.Kan = merge_phyloseq(ps.nor.P.24h.Kan.1, ps.nor.P.24h.Kan.2, ps.nor.P.24h.Kan.3, ps.nor.P.24h.Kan.4)

# plot the top 30 taxa
top30otu <- names(sort(taxa_sums(ps.nor.P.24h.Kan), TRUE)[1:30])
dat.dataframe = psmelt(ps.nor.P.24h.Kan)
dat.dataframe$Genus.x[which(!(dat.dataframe$OTU %in% top30otu))]=NA 

p <- ggplot(dat.dataframe, aes(x=Sample, y=Abundance, fill=Genus.x)) + 
  geom_bar(stat="identity") +  theme(axis.text.x = element_text(angle = 90))

p + guides(fill=guide_legend(ncol = 1,byrow = FALSE, 
                             keywidth =unit(0.4,"cm"),
                             keyheight=unit(0.4,"cm"))) +
  theme(legend.text = element_text(size=12))
p
# ndms plot
ord.nmds.bray = ordinate(ps.nor.P.24h.Kan, method = "NMDS", distance = "bray")
p24K <- plot_ordination(ps.nor.P.24h.Kan, ord.nmds.bray, color="Group", title="Bray NMDS") + theme_bw() 
p24K$layers <- p24K$layers[-1]
p24K <- p24K + geom_point( shape = 1, size = 3, stroke = 1.2) # + xlim(-0.65, 0.81) + ylim(-1.1, 0.43)
p24K
# weighted unifrac plot
ordu = ordinate(ps.nor.P.24h.Kan, "PCoA", "unifrac", weighted=TRUE)
p24K.f <- plot_ordination(ps.nor.P.24h.Kan, ordu, color="Group", title="Weighted Unifrac")+ theme_bw() 
p24K.f$layers <- p24K.f$layers[-1]
p24K.f <- p24K.f + geom_point( shape = 1, size = 3, stroke = 1.2) # + xlim(-0.65, 0.81) + ylim(-1.1, 0.43)
p24K.f

# PERMANOVA
metada <- as(sample_data(ps.nor.P.24h.Kan), "data.frame")
adonis(distance(ps.nor.P.24h.Kan, method = "bray", perm=999) ~ Group,
       data = metada)

#### Pairwise PERMANOVA
z <- pairwise.adonis(distance(ps.nor.P.24h.Kan, method = "bray", perm=999), metada$Group, 
                     sim.method = "bray", p.adjust.m = "bonferroni")
write.csv(z, 'pairwise.adonis.P.24h.Kan.csv')

#####################################################################################
######################################################################## host 2, kanamycin
ps.nor.W = subset_samples(exvivo.pts.nor, Type == "Wf")
ps.nor.W.24h = subset_samples(ps.nor.W, Time == "24h")

ps.nor.W.24h.Kan.1 = subset_samples(ps.nor.W.24h, Group == "Ctl" )
ps.nor.W.24h.Kan.2 = subset_samples(ps.nor.W.24h, Group == "Fcdn")
ps.nor.W.24h.Kan.3 = subset_samples(ps.nor.W.24h, Group == "H.Kan.L.Fcdn")
ps.nor.W.24h.Kan.4 = subset_samples(ps.nor.W.24h, Group == "H.Kan")
ps.nor.W.24h.Kan = merge_phyloseq(ps.nor.W.24h.Kan.1, ps.nor.W.24h.Kan.2, ps.nor.W.24h.Kan.3, ps.nor.W.24h.Kan.4)

top30otu <- names(sort(taxa_sums(ps.nor.W.24h.Kan), TRUE)[1:30])
dat.dataframe = psmelt(ps.nor.W.24h.Kan)
dat.dataframe$Genus.x[which(!(dat.dataframe$OTU %in% top30otu))]=NA 

p <- ggplot(dat.dataframe, aes(x=Sample, y=Abundance, fill=Genus.x)) + 
  geom_bar(stat="identity") +  theme(axis.text.x = element_text(angle = 90))

p + guides(fill=guide_legend(ncol = 2,byrow = FALSE, 
                             keywidth =unit(0.3,"cm"),
                             keyheight=unit(0.3,"cm")))+
  theme(legend.text = element_text(size=9))

# ndms plot
ord.nmds.bray = ordinate(ps.nor.W.24h.Kan, method = "NMDS", distance = "bray")
w24K <- plot_ordination(ps.nor.W.24h.Kan, ord.nmds.bray, color="Group", title="Bray NMDS") + theme_bw() 
w24K$layers <- w24K$layers[-1]
w24K <- w24K + geom_point( shape = 1, size = 3, stroke = 1.2) # + xlim(-0.65, 0.81) + ylim(-1.1, 0.43)
w24K
# weighted unifrac plot
ordu = ordinate(ps.nor.W.24h.Kan, "PCoA", "unifrac", weighted=TRUE)
w24K.f <- plot_ordination(ps.nor.W.24h.Kan, ordu, color="Group", title="Weighted Unifrac")+ theme_bw() 
w24K.f$layers <- w24K.f$layers[-1]
w24K.f <- w24K.f + geom_point( shape = 1, size = 3, stroke = 1.2) # + xlim(-0.65, 0.81) + ylim(-1.1, 0.43)
w24K.f

# PERMANOVA
metada <- as(sample_data(ps.nor.W.24h.Kan), "data.frame")
adonis(distance(ps.nor.W.24h.Kan, method = "bray", perm=999) ~ Group,
       data = metada)

#### Pairwise PERMANOVA
z <- pairwise.adonis(distance(ps.nor.W.24h.Kan, method = "bray", perm=999), metada$Group, 
                     sim.method = "bray", p.adjust.m = "bonferroni")
write.csv(z, 'pairwise.adonis.W.24h.Kan.csv')

#####################################################################################
######################################################################## host 2, ampicillin
ps.nor.W = subset_samples(exvivo.pts.nor, Type == "Wf")
ps.nor.W.24h = subset_samples(ps.nor.W, Time == "24h")

ps.nor.W.24h.Amp.1 = subset_samples(ps.nor.W.24h, Group == "Ctl2" )
ps.nor.W.24h.Amp.2 = subset_samples(ps.nor.W.24h, Group == "Fcdn")
ps.nor.W.24h.Amp.3 = subset_samples(ps.nor.W.24h, Group == "H.Amp.Fcdn")
ps.nor.W.24h.Amp.4 = subset_samples(ps.nor.W.24h, Group == "H.Amp")
ps.nor.W.24h.Amp = merge_phyloseq(ps.nor.W.24h.Amp.1, ps.nor.W.24h.Amp.2, ps.nor.W.24h.Amp.3, ps.nor.W.24h.Amp.4)

top30otu <- names(sort(taxa_sums(ps.nor.W.24h.Amp), TRUE)[1:30])
dat.dataframe = psmelt(ps.nor.W.24h.Amp)
dat.dataframe$Genus.x[which(!(dat.dataframe$OTU %in% top30otu))]=NA 

p <- ggplot(dat.dataframe, aes(x=Sample, y=Abundance, fill=Genus.x)) + 
  geom_bar(stat="identity") +  theme(axis.text.x = element_text(angle = 90))

p + guides(fill=guide_legend(ncol = 2,byrow = FALSE, 
                             keywidth =unit(0.3,"cm"),
                             keyheight=unit(0.3,"cm")))+
  theme(legend.text = element_text(size=9))

# ndms plot
ord.nmds.bray = ordinate(ps.nor.W.24h.Amp, method = "NMDS", distance = "bray")
w24K <- plot_ordination(ps.nor.W.24h.Amp, ord.nmds.bray, color="Group", title="Bray NMDS") + theme_bw() 
w24K$layers <- w24K$layers[-1]
w24K <- w24K + geom_point( shape = 1, size = 3, stroke = 1.2) # + xlim(-0.65, 0.81) + ylim(-1.1, 0.43)
w24K
# weighted unifrac plot
ordu = ordinate(ps.nor.W.24h.Amp, "PCoA", "unifrac", weighted=TRUE)
w24K.f <- plot_ordination(ps.nor.W.24h.Amp, ordu, color="Group", title="Weighted Unifrac")+ theme_bw() 
w24K.f$layers <- w24K.f$layers[-1]
w24K.f <- w24K.f + geom_point( shape = 1, size = 3, stroke = 1.2) # + xlim(-0.65, 0.81) + ylim(-1.1, 0.43)
w24K.f

# PERMANOVA
metada <- as(sample_data(ps.nor.W.24h.Amp), "data.frame")
adonis(distance(ps.nor.W.24h.Amp, method = "bray", perm=999) ~ Group,
       data = metada)

#### Pairwise PERMANOVA
z <- pairwise.adonis(distance(ps.nor.W.24h.Amp, method = "bray", perm=999), metada$Group, 
                     sim.method = "bray", p.adjust.m = "bonferroni")
write.csv(z, 'pairwise.adonis.W.24h.Amp.csv')