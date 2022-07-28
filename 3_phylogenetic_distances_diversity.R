# Heli Juottonen 2021 heli.juottonen@alumni.helsinki.fi
# R scripts for PacBio sequence analysis of the woodtiger moth frass transplant experiment
#  https://github.com/helijuottonen/mothtransplant

# For the article: xxx
# Juottonen H 1, Moghadam NN 1, Murphy L 1, Mappes J 1,2, Galarza JA 1,2 
# 1 Department of Biological and Environmental Sciences, University of Jyväskylä, Finland
# 2 Organismal and Evolutionary Biology Research Program, University of Helsinki, Finland

# 3: Phylogenetic diversity and NMDS with phylogenetic distances (MPD, MNTD)

# Run before this: 1_basic_processing.R, 2_subsampling.R

# Disclaimer: I'm a self-taught R user. I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices. Use at your own risk!

library(picante) # v. 1.8.2
library(Biostrings)
library(ggplot2) # v. 3.3.5
library(cowplot) # v. 1.1.1
library(dplyr) #v. 1.0.8

# extracting sequences for making a tree
mothseq <- taxa_names(mothps3s_seq)
mothseq2 <- Biostrings::DNAStringSet(mothseq)
taxa_names(mothps3s_seq) <- paste0("ASV", seq(ntaxa(mothps3s_seq)))
names(mothseq2) <- taxa_names(mothps3s_seq)
Biostrings::writeXStringSet(mothseq2,  file = "data/mothseq3s.fasta", format = "fasta")

# RAxML tree built in Silva ACT with the sequences in 'mothseq3s.fasta'
# -> mothps3s.tree

#### preparations for phylogenetic distances & betadiversity

mothtree <- read.tree("data/mothps3s.tree")
mothtree
mothtree2 = multi2di(mothtree) #tip taken from Researchgate, tree not read correctly if this is left out
#plot(mothtree2)

# making sure same samples in OTU table (from 2_subsampling.R) and tree
mothtree3 <-prune.sample(mothps3s.asv.st, mothtree2)

# making sure the order of samples is the same in tree and OTU table
mothps3s.asv.st.ord <- mothps3s.asv.st[, mothtree3$tip.label]
mothps3s.asv.st.ord2 <- mothps3s.asv.st.ord[order(row.names(mothps3s.asv.st.ord)),]

####### phylogenetic diversity

pd.result <- pd(mothps3s.asv.st.ord2, mothtree3, include.root = TRUE)
mothmeta.ord.pd <- mothmeta.ord[order(row.names(mothmeta.ord)),]
pd.result$tissue2 <- mothmeta.ord.pd$tissue2
pd.result$tissue <- mothmeta.ord.pd$tissue
pd.result$genotype <- mothmeta.ord.pd$genotype
pd.result$transpl <- mothmeta.ord.pd$transpl

write.csv2(pd.result, "results/pd_result.txt")
pd.result <- read.csv2("results/pd_result.txt", header=T, row.names=1)


# plotting phylogenetic diversity

pd.result$tissue <- factor(pd.result$tissue, c("Whole_Larva", "Poo_before", "Poo_after", "Gut", "AF", "Diet"))

levels(pd.result$tissue) 
levels(pd.result$tissue) <- c("whole larva", "frass before", "frass after", "gut", "abdominal fluid", "diet")

pd.result$tissue[is.na(pd.result$tissue)] <-"diet"

# colours with genotype
moth.gen.pal2 <- c("WW"="#BBBBBB", "yy"="#CCBB44", "Diet" = "white")
pdlabs <- c("ctrl","ctrl","transpl","transpl")

# plot of phylogenetic diversity per sample type
pdplot <- ggplot(pd.result, aes(x=tissue2, y=PD, fill=genotype)) +
  geom_boxplot(outlier.shape=NA, width=0.4)  +
  geom_point(position=position_jitterdodge(jitter.width=0.2), shape=21, aes(fill=genotype)) +
  scale_fill_manual(values=moth.gen.pal2) +
  scale_color_manual(values=moth.gen.pal2) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_x_discrete(labels=pdlabs) +
  facet_grid(cols = vars(tissue), scales="free", space="free") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

pdplot

# plot of ASV numbers per sample type
srplot <- ggplot(pd.result, aes(x=tissue2, y=SR, fill=genotype)) +
  geom_boxplot(outlier.shape=NA, width=0.4)  +
  geom_point(position=position_jitterdodge(jitter.width=0.2), shape=21, aes(fill=genotype)) +
  scale_fill_manual(values=moth.gen.pal2) +
  scale_color_manual(values=moth.gen.pal2) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_x_discrete(labels=pdlabs) +
  facet_grid(cols = vars(tissue), scales="free", space="free") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

srplot

# combined plot
pd_sr_plot <- plot_grid(pdplot, srplot, labels="AUTO", ncol=1)

pd_sr_plot
pdf("Figures/pd_sr_plot.pdf") 
print(pd_sr_plot)
dev.off()

####### phylogenetic distances
# MPD & MNTD

###  calculating MPD and MNTD
mothtree.dist <- cophenetic(mothtree3)

# MPD
comdist.result <- comdist(mothps3s.asv.st.ord2, mothtree.dist, abundance.weighted=TRUE)

# MNTD
comdistnt.result <- comdistnt(mothps3s.asv.st.ord2, mothtree.dist, abundance.weighted=TRUE)


#### nmds with phylogenetic distances

# NMDS with MPD
mdp.nmds <- metaMDS(comdist.result, trymax=200, trace=FALSE, autotransform=FALSE)
mdp.nmds

# NMDS with MNTD
mntd.nmds <- metaMDS(comdistnt.result, trymax=200, trace=FALSE, autotransform=FALSE)
mntd.nmds

# extracting sample scores
mdp.sc.df <- as.data.frame(scores(mdp.nmds))
mntd.sc.df <- as.data.frame(scores(mntd.nmds))

mothmeta.ord <- mothmeta.ord[order(row.names(mothmeta.ord)),]
mdp.sc.df <- mdp.sc.df[order(row.names(mdp.sc.df)),]
mntd.sc.df <- mntd.sc.df[order(row.names(mntd.sc.df)),]

colnames(mdp.sc.df) <- c("mdp_NMDS1", "mdp_NMDS2")
colnames(mntd.sc.df) <- c("mntd_NMDS1", "mntd_NMDS2")
mothps3s.nmds.sc <- merge(mdp.sc.df, mntd.sc.df, by="row.names")
row.names(mothps3s.nmds.sc) <- mothps3s.nmds.sc$Row.names
mothps3s.nmds.sc <- mothps3s.nmds.sc[,2:5]
mothps3s.nmds.sc <- merge(mothps3s.nmds.sc, mothmeta.ord, by="row.names")
row.names(mothps3s.nmds.sc) <- mothps3s.nmds.sc$Row.names

mothps3s.nmds.sc$tissue_transpl <- paste(mothps3s.nmds.sc$tissue, mothps3s.nmds.sc$transpl, sep="_")

mothps3s.nmds.nd <- subset(mothps3s.nmds.sc, transpl != "diet")

mothps3s.nmds.nd$tissue_transpl <- paste(mothps3s.nmds.nd$tissue, mothps3s.nmds.nd$transpl, sep="_")

moth.gen.pal <- c("WW"="#BBBBBB", "yy"="#CCBB44")

# plotting mpd separated by tissue types / experiment stage

wl_pb <- c("Whole_Larva", "Poo_before")
wl_pb.nmds <- subset(mothps3s.nmds.nd, tissue %in% wl_pb)
pb_pa <- c("Poo_before", "Poo_after")
pb_pa.nmds <- subset(mothps3s.nmds.nd, tissue %in% pb_pa)
gut_af <- c("Gut", "AF")
gut_af.nmds <- subset(mothps3s.nmds.nd, tissue %in% gut_af)

wl_pb.plot <- ggplot(wl_pb.nmds, aes(mdp_NMDS1, mdp_NMDS2, shape=tissue, color=genotype)) + 
  geom_point(size=4) +
  theme_cowplot(11) + 
  coord_fixed() +
  xlim(-0.7, 0.3) +
  ylim(-0.3, 0.3) +
  scale_colour_manual(values=moth.gen.pal, name=NULL) + 
  scale_shape_manual(values=c(17,19), name=NULL, labels=c("frass before", "whole larvae")) +
  xlab("NMDS1") +
  ylab("NMDS2")

wl_pb.plot

pb_pa.plot <- ggplot(pb_pa.nmds, aes(mdp_NMDS1, mdp_NMDS2, shape=tissue_transpl, color=genotype)) + 
  geom_point(size=4) +
  theme_cowplot(11) + 
  coord_fixed() +
  xlim(-0.7, 0.3) +
  ylim(-0.3, 0.3) +
  scale_colour_manual(values=moth.gen.pal, name=NULL) +
  scale_shape_manual(values=c(2,6,17), name=NULL, labels=c("frass after ctrl", "frass after transpl", "frass before")) +
  xlab("NMDS1") +
  ylab("NMDS2")

pb_pa.plot

gut_af.plot <- ggplot(gut_af.nmds, aes(mdp_NMDS1, mdp_NMDS2, shape=tissue_transpl, color=genotype)) + 
  geom_point(size=4) +
  theme_cowplot(11) + 
  coord_fixed() +
  xlim(-0.7, 0.3) +
  ylim(-0.3, 0.3) +
  scale_colour_manual(values=moth.gen.pal, name=NULL) + 
  scale_shape_manual(values=c(18,5,20,1), name=NULL, labels=c("abd. fluid ctrl", "abd. fluid transpl", "gut ctrl", "gut transpl")) + 
  xlab("NMDS1") +
  ylab("NMDS2")

gut_af.plot

#  combined plot of MPD NMDS
plot_grid(wl_pb.plot, pb_pa.plot, gut_af.plot, labels="AUTO", ncol=1, align="v")

# mntd separated by tissue type / experiment stage

wl_pb.mntd <- subset(mothps3s.nmds.nd, tissue %in% wl_pb)
pb_pa.mntd <- subset(mothps3s.nmds.nd, tissue %in% pb_pa)
gut_af.mntd <- subset(mothps3s.nmds.nd, tissue %in% gut_af)

wl_pb.plot <- ggplot(wl_pb.mntd, aes(mntd_NMDS1, mntd_NMDS2, shape=tissue, color=genotype)) + 
  geom_point(size=4) +
  theme_cowplot(11) + 
  #coord_fixed() +
  xlim(-0.45, 0.3) +
  ylim(-0.35, 0.35) +
  scale_colour_manual(values=moth.gen.pal, name=NULL) + 
  scale_shape_manual(values=c(17,19), name=NULL, labels=c("frass before", "whole larvae")) +
  xlab("NMDS1") +
  ylab("NMDS2")

wl_pb.plot

pb_pa.plot <- ggplot(pb_pa.mntd, aes(mntd_NMDS1, mntd_NMDS2, shape=tissue_transpl, color=genotype)) + 
  geom_point(size=4) +
  theme_cowplot(11) + 
  #coord_fixed() +
  xlim(-0.45, 0.3) +
  ylim(-0.35, 0.35) +
  scale_colour_manual(values=moth.gen.pal, name=NULL) + 
  scale_shape_manual(values=c(2,6,17), name=NULL, labels=c("frass after ctrl", "frass after transpl", "frass before")) +
  xlab("NMDS1") +
  ylab("NMDS2")

pb_pa.plot

gut_af.plot <- ggplot(gut_af.mntd, aes(mntd_NMDS1, mntd_NMDS2, shape=tissue_transpl, color=genotype)) + 
  geom_point(size=4) +
  theme_cowplot(11) + 
  #coord_fixed() +
  xlim(-0.45, 0.3) +
  ylim(-0.35, 0.35) +
  scale_colour_manual(values=moth.gen.pal, name=NULL) + 
  scale_shape_manual(values=c(18,5,20,1), name=NULL, labels=c("abd. fluid ctrl", "abd. fluid transpl", "gut ctrl", "gut transpl")) + 
  xlab("NMDS1") +
  ylab("NMDS2")

gut_af.plot

# combined plot of MNTD NMDS
plot_grid(wl_pb.plot, pb_pa.plot, gut_af.plot, labels="AUTO", ncol=1, align="v")

