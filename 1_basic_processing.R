# Heli Juottonen 2021 heli.juottonen@alumni.helsinki.fi
# R scripts for PacBio sequence analysis of the woodtiger moth frass transplant experiment
#  https://github.com/helijuottonen/mothtransplant

# For the article: xxx
# Juottonen H 1, Moghadam NN 1, Murphy L 1, Mappes J 1,2, Galarza JA 1,2 
# 1 Department of Biological and Environmental Sciences, University of Jyväskylä, Finland
# 2 Organismal and Evolutionary Biology Research Program, University of Helsinki, Finland

# 1: Basic processing of data:
# reading in data, modifying names, removing unwanted taxa/samples,
# making a phyloseq object
# starts from ASV & taxonomy tables (in the folder results) created with 0_dada2.R

# Disclaimer: I'm a self-taught R user. I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices. Use at your own risk!

library(phyloseq)
# v. 1.32.0

# reading in data
# ASV table from dada2
st.nochim <- readRDS("results/RDS/moth_transpl_st_nochim.rds")
# taxonomy from dada2
tax <- readRDS("results/RDS/moth_transpl_tax_138.1.rds")
# metadata
mothmeta <- read.csv2("data/moth_transpl_metadata.csv", header=T, row.names=1)

# adding more informative samples names to metadata
mothmeta$sample3 <- paste(mothmeta$tissue, mothmeta$transpl, mothmeta$m.family, mothmeta$sample_ID, sep="_")
mothmeta$sample4 <- paste(mothmeta$m.family, mothmeta$sample_ID, sep="_")
mothmeta$tissue_gen <- paste(mothmeta$tissue, mothmeta$genotype, sep="_")
mothmeta$tissue_tr_gen <- paste(mothmeta$tissue, mothmeta$transpl, mothmeta$genotype, sep="_")

# replacing NAs in taxonomy table (issue with Silva 138.1 in dada2)
# if family missing -> order, if genus missing -> family
tax.df <- as.data.frame(tax)
tax.df$Order <- ifelse(is.na(tax.df$Order), as.character(tax.df$Class), as.character(tax.df$Order))
tax.df$Family <- ifelse(is.na(tax.df$Family), as.character(tax.df$Order), as.character(tax.df$Family))
tax.df$Genus <- ifelse(is.na(tax.df$Genus), as.character(tax.df$Family), as.character(tax.df$Genus))

# reading data into phyloseq
otups = otu_table(st.nochim, taxa_are_rows = FALSE)
taxps = tax_table(as.matrix(tax.df))
metadata <- sample_data(mothmeta)

# modifying sample names in the OTU table to match row names of metadata
sample_names(otups) <- gsub(".fastq", "", sample_names(otups))
sample_names(otups) <- gsub("Lib_1", "Lib1", sample_names(otups))
sample_names(otups) <- gsub("Lib_2", "Lib2", sample_names(otups))
sample_names(otups) <- gsub("Lib_3", "Lib3", sample_names(otups))
sample_names(otups) <- gsub("per_Bc.bc1002R--", "", sample_names(otups))
sample_names(otups) <- gsub("per_Bc.", "", sample_names(otups))
sample_names(otups) <- gsub("--bc1002R", "", sample_names(otups))
sample_names(otups) <- gsub("bc10", "Bc", sample_names(otups))

# checking to make sure sample names and metadata row names match
setdiff(sample_names(otups), row.names(mothmeta))
# should print 'character(0)'

# creating and checking phyloseq object
mothps = phyloseq(otups, taxps, metadata)
mothps

# removing mitochondria, chloroplasts and non-bacteria
ntaxa(mothps) 
mothps1m <- subset_taxa(mothps, (Family!="Mitochondria") | is.na(Family))
ntaxa(mothps1m)

mothps1mc <- subset_taxa(mothps1m, (Order!="Chloroplast") | is.na(Order))
ntaxa(mothps1mc) 

mothps1mcn <- subset_taxa(mothps1mc, (Kingdom!="NA"))
ntaxa(mothps1mcn) 

# removing bad samples
# IMPORTANT NOTE: these samples are not in Genbank - only do this step
# if using the ASV table from Github/ results!:

# repeats with low read numbers
mothps2 = subset_samples(mothps1mcn, sample_names(mothps1mcn) != "Lib3_Bc40")
mothps2 = subset_samples(mothps2, sample_names(mothps2) != "Lib2_Bc41")

# other repeats/dilutions
mothps2 = subset_samples(mothps2, sample3 != "Poo_after_ctrl_WW2_140dil")
mothps2 = subset_samples(mothps2, sample3 != "Poo_after_ctrl_WW4_145")
mothps2 = subset_samples(mothps2, sample3 != "Poo_after_transpl_WW2_166dil")
mothps2 = subset_samples(mothps2, sample3 != "Poo_after_transpl_WW3_164b")
mothps2 = subset_samples(mothps2, sample3 != "Poo_after_transpl_WW4_162a")
mothps2 = subset_samples(mothps2, sample3 != "Poo_after_transpl_yy5_160dil")
mothps2 = subset_samples(mothps2, sample3 != "Whole_Larva_before_WW3_L8")

# removing ASVs with no reads
mothps2 <- filter_taxa(mothps2, function(x) sum(x) >0, TRUE)
# checking what is left
mothps2

# removing mock community and negative controls
mothps3 <- subset_samples(mothps2, transpl != "negative" & transpl != "positive")

#removing taxa with no reads after subsetting
mothps3 <- filter_taxa(mothps3, function(x) sum(x) >0, TRUE)

# check the number of samples and taxa
mothps3

# phyloseq object without washwater samples
mothps3s <- subset_samples(mothps3, tissue != "Water_before")
mothps3s <- filter_taxa(mothps3s, function(x) sum(x) >0, TRUE)
mothps3s # 161 taxa, 61 samples

# copy of the phyloseq object for extracting sequences in 3_phylogenetic_distances_diversity.R
mothps3s_seq <- mothps3s

