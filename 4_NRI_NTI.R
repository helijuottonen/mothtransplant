# Heli Juottonen 2021 heli.juottonen@alumni.helsinki.fi
# R scripts for PacBio sequence analysis of the woodtiger moth frass transplant experiment
#  https://github.com/helijuottonen/mothtransplant

# For the article: xxx
# Juottonen H 1, Moghadam NN 1, Murphy L 1, Mappes J 1,2, Galarza JA 1,2 
# 1 Department of Biological and Environmental Sciences, University of Jyväskylä, Finland
# 2 Organismal and Evolutionary Biology Research Program, University of Helsinki, Finland

# 4: NRI, NTI

# run first:
#   1_basic_processing.R
#   2_subsampling.R
#   3_phylogenetic_distances_diversity.R

# Disclaimer: I'm a self-taught R user. I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices. Use at your own risk!

# library(picante)
# library(ggplot2)
# library(cowplot)
# library(dplyr)

library(reshape2) # v. 1.4.4

# calculate NRI, NTI with picante

nri.res <- ses.mpd(mothps3s.asv.st.ord2, mothtree.dist, abundance.weighted=TRUE, null.model ="taxa.labels", runs = 999)

nti.res <- ses.mntd(mothps3s.asv.st.ord2, mothtree.dist, abundance.weighted=TRUE, null.model ="taxa.labels", runs = 999)

# combine the results with metadata and save in a file

mothmeta.ord <- mothmeta.ord[order(row.names(mothmeta.ord)),]
nri.res.ord <- nri.res[order(row.names(nri.res)),]
nti.res.ord <- nti.res[order(row.names(nti.res)),]
nti.res.ord <- subset(nti.res.ord, select = -c(ntaxa, runs))

nri.nti.ord <- cbind(nri.res.ord, nti.res.ord)
nri.nti.ord$tissue <- mothmeta.ord$tissue
nri.nti.ord$tissue_gen <- mothmeta.ord$tissue_gen
nri.nti.ord$tissue2 <- mothmeta.ord$tissue2
nri.nti.ord$tissue_tr_gen <- mothmeta.ord$tissue_tr_gen
nri.nti.ord$genotype <- mothmeta.ord$genotype

write.csv2(nri.nti.ord, "results/nri_nti.txt")
nri.nti.ord <- read.csv2("results/nri_nti.txt")


# nri, nti plots 

# define treatment names and order
nri.nti.ord.nd <- na.omit(nri.nti.ord)
nri.nti.ord.nd$tissue <- factor(nri.nti.ord.nd$tissue, levels= c("Whole_Larva", "Poo_before", "Poo_after", "Gut", "AF"), labels=c("whole larva", "frass before", "frass after", "gut", "abdominal fluid"))

# define colours
moth.gen.pal <- c("WW"="#BBBBBB", "yy"="#CCBB44")

# plotting
nriplot <- ggplot(nri.nti.ord.nd, aes(x=tissue, y=(mpd.obs.z)*-1, fill=genotype)) +
  geom_hline(yintercept=0, color="black", linetype="dotted") +
  geom_boxplot(outlier.shape=NA, width=0.4)  +
  geom_point(position=position_jitterdodge(jitter.width=0.1), shape=21, aes(fill=genotype)) +
  scale_fill_manual(values=moth.gen.pal) +
  scale_color_manual(values=moth.gen.pal) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  facet_grid(cols = vars(tissue), scales="free", space="free") +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(y = "NRI") 

nriplot

ntiplot <- ggplot(nri.nti.ord.nd, aes(x=tissue, y=(mntd.obs.z)*-1, fill=genotype)) +
  geom_hline(yintercept=0, color="black", linetype="dotted") +
  geom_boxplot(outlier.shape=NA, width=0.4)  +
  geom_point(position=position_jitterdodge(jitter.width=0.1), shape=21, aes(fill=genotype)) +
  scale_fill_manual(values=moth.gen.pal) +
  scale_color_manual(values=moth.gen.pal) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  facet_grid(cols = vars(tissue), scales="free", space="free") +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(y = "NTI") 

ntiplot

# combine plots with cowplot
nri_nti_plot <- plot_grid(nriplot, ntiplot, labels="AUTO", ncol=1)


# Do the NTI and NRI values differ statistically from 0?
# differs - phylogenetic clustering more than by chance/ with neutral processes

wl.t <- nri.nti.ord %>% filter(tissue == "Whole_Larva")
t.test(wl.t$mntd.obs.z) 
t.test(wl.t$mpd.obs.z) 

pb.t <- nri.nti.ord %>% filter(tissue == "Poo_before")
t.test(pb.t$mntd.obs.z) 
t.test(pb.t$mpd.obs.z) 

pa.t <- nri.nti.ord %>% filter(tissue == "Poo_after")
t.test(pa.t$mntd.obs.z) 
t.test(pa.t$mpd.obs.z) 

gut.t <- nri.nti.ord %>% filter(tissue == "Gut")
t.test(gut.t$mntd.obs.z) 
t.test(gut.t$mpd.obs.z) 

af.t <- nri.nti.ord %>% filter(tissue == "AF")
t.test(af.t$mntd.obs.z) 
t.test(af.t$mpd.obs.z)


# betaNTI

# Stegen et al. 2013 ISMEJ https://www.nature.com/articles/ismej201393
# script from: https://github.com/stegen/Stegen_etal_ISME_2013
# library(picante)

## read in OTU table
otu = mothps3s.asv.st;
dim(otu); # this gives the dimensions
otu[1:5,1:5]; # this gives a look at the first 5 rows and columns

## read in the phylogeny
phylo = mothtree3
phylo; # a summary of the phylogeny
plot.phylo(phylo,typ="fan"); # a quick plot

## make sure the names on the phylogeny are ordered the same as the names in otu table

match.phylo.otu = match.phylo.data(phylo, t(otu));
str(match.phylo.otu);

## calculate empirical betaMNTD

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
write.csv2(beta.mntd.weighted,'results/betaMNTD_weighted.csv',quote=F);

identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

# calculate randomized betaMNTD

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv2(weighted.bNTI,"results/weighted.bNTI.txt",quote=F);

pdf("Figures/weighted_bNTI_Histogram.pdf")
hist(weighted.bNTI)
dev.off()

# stegen 2013 code ends

weighted.bNTI <- read.csv2("results/weighted.bNTI.txt", row.names=1, header=T)

mothps3s.bNTI <- melt(as.matrix(weighted.bNTI), varnames = c("row", "col"))
mothps3s.bNTI <- na.omit(mothps3s.bNTI)
mothps3s.bNTI <- mothps3s.bNTI %>% rename(bNTI = value)


# betaNTI for tissue pairs along the experiment

# whole larva vs. poo before
wl_pb <- mothps3s.bNTI %>% filter(grepl('Whole_Larva', col) & grepl('Poo_before', row))

wl_pb_ww <- wl_pb %>% filter(grepl('WW', col) & grepl('WW', row)) %>% mutate(genotype ="WW")
wl_pb_yy <- wl_pb %>% filter(grepl('yy', col) & grepl('yy', row)) %>% mutate(genotype ="yy")
wl_pb.bnti <- dplyr::bind_rows(list(wl_pb_ww=wl_pb_ww, wl_pb_yy=wl_pb_yy), .id = 'transition')

# poo before vs. poo after
pb_pa <- mothps3s.bNTI %>% filter(grepl('Poo_before', col) & grepl('Poo_after', row))

pb_pa_ww <- pb_pa %>% filter(grepl('WW', col) & grepl('WW', row)) %>% mutate(genotype ="WW")
pb_pa_yy <- pb_pa %>% filter(grepl('yy', col) & grepl('yy', row)) %>% mutate(genotype ="yy")

pb_pa.bnti <- dplyr::bind_rows(list(pb_pa_ww=pb_pa_ww, pb_pa_yy=pb_pa_yy), .id = 'transition')

# poo after vs. gut
pa_gut <- mothps3s.bNTI %>% filter(grepl('Poo_after', col) & grepl('Gut', row))

pa_gut_ww <- pa_gut %>% filter(grepl('WW', col) & grepl('WW', row)) %>% mutate(genotype ="WW")
pa_gut_yy <- pa_gut %>% filter(grepl('yy', col) & grepl('yy', row)) %>% mutate(genotype ="yy")

pa_gut.bnti <- dplyr::bind_rows(list(pa_gut_ww=pa_gut_ww, pa_gut_yy=pa_gut_yy), .id = 'transition')

# gut vs. AF
gut_AF <- mothps3s.bNTI %>% filter(grepl('Gut', col) & grepl('AF', row))

gut_AF_ww <- gut_AF %>% filter(grepl('WW', col) & grepl('WW', row)) %>% mutate(genotype ="WW")
gut_AF_yy <- gut_AF %>% filter(grepl('yy', col) & grepl('yy', row)) %>% mutate(genotype ="yy")

gut_AF.bnti <- dplyr::bind_rows(list(gut_AF_ww=gut_AF_ww, gut_AF_yy=gut_AF_yy), .id = 'transition')


all.bnti.cross <- dplyr::bind_rows(list(wl_pb=wl_pb.bnti, pb_pa=pb_pa.bnti, pa_gut=pa_gut.bnti, gut_AF=gut_AF.bnti), .id = 'transition')

all.bnti.cross$transition <- factor(all.bnti.cross$transition, c("wl_pb", "pb_pa", "pa_gut", "gut_AF"))
moth.gen.pal <- c("WW"="#BBBBBB", "yy"="#CCBB44")

levels(all.bnti.cross$transition) <- c("whole larva-frass", "frass before-after", "frass after-gut", "gut-abdominal fluid")

# plotting betaNTI for tissue pairs along the experiment

bnti.cross.plot <- ggplot(all.bnti.cross, aes(x=transition, y=bNTI)) +
  geom_hline(yintercept=-2, color="black", linetype="dotted") +
  geom_hline(yintercept=2, color="black", linetype="dotted") +
  geom_boxplot(aes(fill=genotype), outlier.shape=NA, width=0.4)  +
  geom_point(position=position_jitterdodge(jitter.width=0.1), shape=21, aes(fill=genotype)) +
  scale_fill_manual(values=moth.gen.pal) +
  scale_color_manual(values=moth.gen.pal) +
  facet_grid(cols = vars(transition), scales="free", space="free") +
  ylim(-3, 5) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("betaNTI")

bnti.cross.plot

# NRI, NTI, betaNTI plots together
nri_nti_bnti_plot <- plot_grid(nriplot, ntiplot, bnti.cross.plot, labels="AUTO", ncol=1)

nri_nti_bnti_plot
pdf("Figures/nri_nti_bnti_plot.pdf") 
print(nri_nti_bnti_plot)
dev.off()

