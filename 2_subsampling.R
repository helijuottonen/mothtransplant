# Heli Juottonen 2021 heli.juottonen@alumni.helsinki.fi
# R scripts for PacBio sequence analysis of the woodtiger moth frass transplant experiment
#  https://github.com/helijuottonen/mothtransplant

# For the article: xxx
# Juottonen H 1, Moghadam NN 1, Murphy L 1, Mappes J 1,2, Galarza JA 1,2 
# 1 Department of Biological and Environmental Sciences, University of Jyväskylä, Finland
# 2 Organismal and Evolutionary Biology Research Program, University of Helsinki, Finland

# 2: Subsampling the ASV table 
# purpose: to take into account differences in sequence read numbers among samples

# Run before this: 1_basic_processing.R

# Disclaimer: I'm a self-taught R user. I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices. Use at your own risk!


library(vegan)
# v.2.5.7

# replacing long ASV names
taxa_names(mothps3s) <- paste0("ASV", seq(ntaxa(mothps3s)))

# splitting phyloseq object for subsampling
mothps3s.asv = as(otu_table(mothps3s), "matrix")
mothps3s.asv.df = as.data.frame(mothps3s.asv)

mothmeta.ord <- mothmeta[order(row.names(mothmeta)),]
mothps3s.asv.ord <- mothps3s.asv.df[order(row.names(mothps3s.asv.df)),]
mothmeta.ord = mothmeta.ord[row.names(mothmeta.ord) %in% row.names(mothps3s.asv.ord),]

rownames(mothps3s.asv.ord) <- mothmeta.ord$sample3
rownames(mothmeta.ord) <- mothmeta.ord$sample3

# subsampling

min(rowSums(mothps3s.asv.ord))
median(rowSums(mothps3s.asv.ord))

# function for rarefying to median read number
# in samples with less reads than median, all reads are kept
median_raref <- function(x) {
  x_sums <- rowSums(x)
  x_sums2 <- replace(x_sums, x_sums>(median(x_sums)), median(x_sums))
  set.seed(1752)
  x.r <- rrarefy(x, x_sums2)
  x.rdf <- as.data.frame(x.r)
  return(x.rdf)
}

mothps3s.asv.r <- median_raref(mothps3s.asv.ord)
mothps3s.asv.r.df <- as.data.frame(mothps3s.asv.r)

# converting to relative abundances (because median used as basis of subsampling)
# this object used in many further analyses
mothps3s.asv.st <- decostand(mothps3s.asv.r.df, "total", margin=1)


# rebuilding phyloseq object with subsampled OTU table

mothps3s.tax = as(tax_table(mothps3s), "matrix")
mothps3s.tax.df = as.data.frame(mothps3s.tax)

otups2 = otu_table(mothps3s.asv.st, taxa_are_rows = FALSE)
taxps2 = tax_table(as.matrix(mothps3s.tax.df))
metadata2 <- sample_data(mothmeta.ord)

mothps3s2 <- phyloseq(otups2, taxps2, metadata2)
# this phyloseq object used in many further analyses
mothps3s2
# 161 taxa, 61 samples
