# Heli Juottonen 2021 heli.juottonen@alumni.helsinki.fi
# R scripts for PacBio sequence analysis of the woodtiger moth frass transplant experiment
#  https://github.com/helijuottonen/mothtransplant

# For the article: xxx
# Juottonen H 1, Moghadam NN 1, Murphy L 1, Mappes J 1,2, Galarza JA 1,2 
# 1 Department of Biological and Environmental Sciences, University of Jyväskylä, Finland
# 2 Organismal and Evolutionary Biology Research Program, University of Helsinki, Finland


# 6: PERMANOVAs based on phylogenetic distances

# run first:
#   1_basic_processing.R
#   2_subsampling.R
#   3_phylogenetic_distances_diversity.R

# Disclaimer: I'm a self-taught R user. I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices. Use at your own risk!


# library(vegan)
# library(picante)

### whole larvae, effect of genotype
# selecting only whole larvae in metadata, subsampled ASV table and tree
mothmeta.ord.pd.wl <- subset(mothmeta.ord.pd, tissue == "Whole_Larva")
mothps3s.asv.st.wl =  mothps3s.asv.st.ord2[row.names(mothps3s.asv.st.ord2) %in% row.names(mothmeta.ord.pd.wl),]
mothtree_wl <-prune.sample(mothps3s.asv.st.wl, mothtree2) 

# calculating phylogenetic distances (MNTD)
mothtree_wl.dist <- cophenetic(mothtree_wl)
comdistnt_wl.result <- comdistnt(mothps3s.asv.st.wl, mothtree_wl.dist, abundance.weighted=TRUE)

# PERMDISP to check for multivariate heterogeneity
wl.bedisp <- betadisper(comdistnt_wl.result, mothmeta.ord.pd.wl$genotype)
anova(wl.bedisp)

# PERMANOVA with adonis2 in vegan
adonis2(comdistnt_wl.result ~ genotype, permutations=999, data=mothmeta.ord.pd.wl)

# frass before+after, effect of genotype
mothmeta.ord.pd.pba <- subset(mothmeta.ord.pd, tissue == "Poo_before" | tissue == "Poo_after")
mothps3s.asv.st.pba =  mothps3s.asv.st.ord2[row.names(mothps3s.asv.st.ord2) %in% row.names(mothmeta.ord.pd.pba),]
mothtree_pba <-prune.sample(mothps3s.asv.st.pba, mothtree2) 
mothtree_pba.dist <- cophenetic(mothtree_pba)
comdistnt_pba.result <- comdistnt(mothps3s.asv.st.pba, mothtree_pba.dist, abundance.weighted=TRUE)

pba.bedisp <- betadisper(comdistnt_pba.result, mothmeta.ord.pd.pba$genotype)
anova(pba.bedisp)

adonis2(comdistnt_pba.result ~ genotype + tissue, by="margin", permutations=999, data=mothmeta.ord.pd.pba)

# frass after, effect of genotype and transplantation
mothmeta.ord.pd.pa <- subset(mothmeta.ord.pd, tissue == "Poo_after")
mothps3s.asv.st.pa =  mothps3s.asv.st.ord2[row.names(mothps3s.asv.st.ord2) %in% row.names(mothmeta.ord.pd.pa),]
mothtree_pa <-prune.sample(mothps3s.asv.st.pa, mothtree2) 
mothtree_pa.dist <- cophenetic(mothtree_pa)

comdist_pa.result <- comdist(mothps3s.asv.st.pa, mothtree_pa.dist, abundance.weighted=TRUE)
adonis2(comdist_pa.result ~ transpl + genotype, by="margin", permutations=999, data=mothmeta.ord.pd.pa)

pa.bedisp <- betadisper(comdist_pa.result, mothmeta.ord.pd.pa$genotype)
anova(pa.bedisp)


