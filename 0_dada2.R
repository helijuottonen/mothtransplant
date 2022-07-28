# Heli Juottonen 2021 heli.juottonen@alumni.helsinki.fi
# R scripts for PacBio sequence analysis of the woodtiger moth frass transplant experiment
#  https://github.com/helijuottonen/mothtransplant

# For the article: xxx
# Juottonen H 1, Moghadam NN 1, Murphy L 1, Mappes J 1,2, Galarza JA 1,2 
# 1 Department of Biological and Environmental Sciences, University of Jyväskylä, Finland
# 2 Organismal and Evolutionary Biology Research Program, University of Helsinki, Finland

# 0: Analysis of PacBio 16S rRNA gene sequence data with dada2:
# from fastq files to sequence and taxonomy table

# Disclaimer: I'm a self-taught R user. I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices. Use at your own risk!

# from: https://benjjneb.github.io/LRASManuscript/LRASms_fecal.html

# load dada2 package
library(dada2); packageVersion("dada2")
# v. 1.16.0

# setting up file paths etc.

path <- "data/fastq"
list.files(path)

path.out <- "results/Figures"
path.rds <- "results/RDS"

fns1 <- list.files(path, pattern=".fastq", full.names=TRUE)

# primer sequences:
F27 <- "AGAGTTTGATCMTGGCTCAG"
R1492 <- "CCTTGTTACGACTTCACCCCAG"

rc <- dada2:::rc


# removing primers

nops <- file.path(path, "noprimers", basename(fns1))
prim <- removePrimers(fns1, nops, primer.fwd=F27, primer.rev=dada2:::rc(R1492), orient=TRUE)

saveRDS(nops, file.path(path.rds, "nops.rds"))
nops <- readRDS(file.path(path.rds, "nops.rds"))

# checking the read length histogram after primer removal
lens.fn <- lapply(nops, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)

#filtering based on quality settings etc.
filts <- file.path(path, "noprimers", "filtered", basename(fns1))
track <- filterAndTrim(nops, filts2, minQ=3, minLen=1300, maxLen=1600, maxN=0, rm.phix=FALSE, maxEE=2)

#check how many sequences left
track

# de-replicating reads
drp <- derepFastq(filts, verbose=TRUE)

# learning and inspecting errors
err <- learnErrors(drp, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)
saveRDS(err, file.path(path.rds, "moth_transpl.rds"))
plotErrors(err)
# looks ok based on the dada2 tutorials

# denoising
dd2 <- dada(drp, err=err, BAND_SIZE=32, multithread=TRUE)

saveRDS(dd2, file.path(path.rds, "moth_transpl_dd2.rds"))

cbind(ccs=prim[,1], primers=prim[,2], filtered=track[,2], denoised=sapply(dd2, function(x) sum(x$denoised)))

# making sequence table
st <- makeSequenceTable(dd2); dim(st)

# saving sequence table
saveRDS(st, file.path(path.rds, "moth_transpl_st.rds"))
# reading in the previously saved sequence table
st <- readRDS(file.path(path.rds, "moth_transpl_st.rds"))

# removing chimeras

# checking how many chimeras are detected
# minFoldParentOverAbundance setting from dada2 PacBio example
# default 2, doesn't make a big difference
bimst <- isBimeraDenovo(st, minFoldParentOverAbundance=3.5)
table(bimst)

#removing chimeras
st.nochim <- removeBimeraDenovo(st, method="consensus", minFoldParentOverAbundance=3.5, verbose=TRUE)
dim(st.nochim)

# checking what was removed
sum(st.nochim)/sum(st)

# saving sequence table without chimeras
saveRDS(st.nochim, file.path(path.rds, "moth_transpl_st_nochim.rds"))

# assign taxonomy against Silva 138.1:
tax <- assignTaxonomy(st.nochim, "data/silva_nr99_v138.1_train_set.fa.gz")
# saving taxonomy assignments
saveRDS(tax, file.path(path.rds, "moth_transpl_tax_138.1.rds"))



