# Heli Juottonen 2021 heli.juottonen@alumni.helsinki.fi
# R scripts for PacBio sequence analysis of the woodtiger moth frass transplant experiment
#  https://github.com/helijuottonen/mothtransplant

# For the article: xxx
# Juottonen H 1, Moghadam NN 1, Murphy L 1, Mappes J 1,2, Galarza JA 1,2 
# 1 Department of Biological and Environmental Sciences, University of Jyväskylä, Finland
# 2 Organismal and Evolutionary Biology Research Program, University of Helsinki, Finland

# 7: Core microbiome 7 shared ASVs

# run first:
#   1_basic_processing.R
#   2_subsampling.R

# Disclaimer: I'm a self-taught R user. I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices. Use at your own risk!

# based on https://microbiome.github.io/tutorials/core_venn.html
# package 'microbiome' (needs also eulerr)

# run first:
#   0_dada2.R
#   1_basic_processing.R
#   2_subsampling.R

library(microbiome)
library(eulerr)

# selecting the sample types
tissues <- c("Whole_Larva", "Poo_before", "Poo_after", "AF", "Gut")

list_core <- c() # an empty object to store information

for (n in tissues){ 
  mothps3s2.sub <- subset_samples(mothps3s2, tissue == n) 
  
  core_m <- core_members(mothps3s2.sub, 
                         detection = 0.00001, # 0.001 in at least 90% samples 
                         prevalence = 0.5)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each tissue
  list_core[[n]] <- core_m # add to a list core taxa for each group.
}

# checking which ASVs are the core members per tissue type

# in whole larvae
mothps3.WL <- subset_samples(mothps3s2, tissue=="Whole_Larva")
core.WL <- core_members(mothps3.WL, detection = 0.00001, prevalence = 50/100)
core.WL

# frass before
mothps3.pb <- subset_samples(mothps3s2, tissue=="Poo_before")
core.pb <- core_members(mothps3.pb, detection = 0.00001, prevalence = 50/100)
core.pb

# frass after
mothps3.pa <- subset_samples(mothps3s2, tissue=="Poo_after")
core.pa <- core_members(mothps3.pa, detection = 0.00001, prevalence = 50/100)
core.pa

# abdominal fluid
mothps3.af <- subset_samples(mothps3s2, tissue=="AF")
core.af <- core_members(mothps3.af, detection = 0.00001, prevalence = 50/100)
core.af

# gut
mothps3.g <- subset_samples(mothps3s2, tissue=="Gut")
core.g <- core_members(mothps3.af, detection = 0.00001, prevalence = 50/100)
core.g

# core members overall with 50% prevalence
mothps3.tis <- subset_samples(mothps3s2, tissue %in% tissues)
core.tis <- core_members(mothps3.tis, detection = 0.00001, prevalence = 50/100)
core.tis


