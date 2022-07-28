# Heli Juottonen 2021 heli.juottonen@alumni.helsinki.fi
# R scripts for PacBio sequence analysis of the woodtiger moth frass transplant experiment
#  https://github.com/helijuottonen/mothtransplant

# For the article: xxx
# Juottonen H 1, Moghadam NN 1, Murphy L 1, Mappes J 1,2, Galarza JA 1,2 
# 1 Department of Biological and Environmental Sciences, University of Jyväskylä, Finland
# 2 Organismal and Evolutionary Biology Research Program, University of Helsinki, Finland

# 5: Taxonomy plots: 
# ASV level bubble plot and stacked bar plot for frass,
# supplementary figures?

# run first:
#   1_basic_processing.R
#   2_subsampling.R

# Disclaimer: I'm a self-taught R user. I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices. Use at your own risk!


# library(phyloseq)
library(ggplot2) # v.3.3.5
library(dplyr) # v. 1.0.8
library(forcats) # v.0.5.1
library(tidyr) # v. 1.2.0

#### ASV level bubbleplot

# starting with subsampled data from 2_subsampling.R

# removing diet samples
mothps3s2b <- subset_samples(mothps3s2, tissue != "Control_Diet" & tissue != "WW_diet" & tissue != "yy_Diet")

# removing ASVs that don't occur in more than 2 samples
mothps3sb2 = filter_taxa(mothps3s2b, function(x){sum(x > 0) > 2}, prune = TRUE)
mothps3sb2 #31 taxa

# ASV level

moth_asv <- mothps3sb2 %>%
  psmelt() %>%                              # Melt to long format
  arrange(OTU)  

moth_asv$tissue_tr_gen <- factor(moth_asv$tissue_tr_gen, c("Whole_Larva_before_WW", "Whole_Larva_before_yy", "Poo_before_before_WW", "Poo_before_before_yy", "Poo_after_ctrl_WW", "Poo_after_transpl_WW", "Poo_after_ctrl_yy", "Poo_after_transpl_yy", "Gut_ctrl_WW", "Gut_transpl_WW", "Gut_ctrl_yy", "Gut_transpl_yy","AF_ctrl_WW", "AF_transpl_WW", "AF_ctrl_yy", "AF_transpl_yy"))

# defining which Enterococcus ASVs are E. mundtii (=mund), which E. casseliflavus/gallinarum (=cass)
# checked from tree
mund <- c("ASV47", "ASV43", "ASV44", "ASV53", "ASV57", "ASV45", "ASV12", "ASV38", "ASV23", "ASV3", "ASV40", "ASV65", "ASV74", "ASV75", "ASV84", "ASV100", "ASV113")
cass <- c("ASV48", "ASV46", "ASV42", "ASV36", "ASV68", "ASV41", "ASV59", "ASV5", "ASV1", "ASV62", "ASV50", "ASV51", "ASV66", "ASV79", "ASV88", "ASV89", "ASV93", "ASV94", "ASV95", "ASV96", "ASV97", "ASV98", "ASV101", "ASV103", "ASV106", "ASV109", "ASV110", "ASV111", "ASV112", "ASV114", "ASV115", "ASV117", "ASV118", "ASV119", "ASV121", "ASV122", "ASV124", "ASV130")

moth_asv2 <- mutate(moth_asv, genus2 = ifelse(OTU %in% mund, "Enteroc. mundtii", ifelse(OTU %in% cass, "Enteroc. cass/gall", Genus)))

# defining tissue + treatment labels
blabs <- c("whole larva", "frass before", "frass ctrl", "frass transpl", "gut ctrl", "gut transpl", "ab. fluid ctrl", "ab. fluid transpl")

# defining colours for sample types (=tissue)
tissuecol <- c("#F6C141", "#4EB265", "#882E72", "#882E72", "#5289C7")

# calculating mean ASV abundances for each sample type (=tissue) + genotype + treatment combination
# + removing ASVs with mean = 0, removing repeated entries
moth_asv_sum <- moth_asv2 %>% 
  group_by(OTU, tissue_tr_gen) %>%
  summarize(ASVmean=mean(Abundance), across(c(genotype, tissue, Phylum, genus2))) %>%
  ungroup() %>%
  filter(ASVmean>0) %>%
  distinct()
  
# ordering ASVs by descending mean abundance
moth_asv_sum <- moth_asv_sum %>% mutate(OTU = fct_reorder(OTU, ASVmean))
  
# plotting
mothbubble_asv <- ggplot(moth_asv_sum, aes(x=tissue_tr_gen, y=OTU)) + 
  geom_point(aes(size = ASVmean, colour= tissue)) +
  scale_size_area() +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(Phylum + genus2 ~ genotype, scales="free", space="free") +
  scale_x_discrete(labels=blabs) +
  guides(color = guide_legend(override.aes = list(size = 3), reverse=TRUE)) +
  scale_color_manual(values = tissuecol, labels = c("abdominal fluid", "gut", "frass after", "frass before", "whole larva")) +
  labs(y="ASV", x="tissue/treatment") +
  theme(strip.text = element_text(colour = 'black'), panel.border = element_blank())

mothbubble_asv

pdf("Figures/mothbubble_asv.pdf") 
print(mothbubble_asv)
dev.off()


### frass before after, ASV level

# frass before and after (note that frass is called 'poo' in places in this script!)
mothps3_pba <- subset_samples(mothps3s2, tissue == "Poo_before" | tissue == "Poo_after")
mothps3_pba
# 161 taxa, 13 samples

ps3_pba_asv <- mothps3_pba %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.005) %>%                         # Filter out low abundance taxa
  arrange(Family) 

ps3_pba_asv2 <- mutate(ps3_pba_asv, family2 = ifelse(OTU %in% mund, "Enteroc. mundtii", ifelse(OTU %in% cass, "Enteroc. cass/gall", Family)))

fam_colors2 <- c("0319-6G20" = "#882E72", "Acetobacteraceae" = "#F6C141", "Anaerovoracaceae" = "#F6C141", "Bacillaceae" = "#F1932D", "Beijerinckiaceae" = "#B178A6", "Caulobacteraceae" = "darkgrey", "Chitinophagaceae" = "#D6C1DE", "Comamonadaceae" = "#1965B0", "Defluviicoccaceae" = "#777777", "Enterobacteriaceae" = "#DDD8EF", "Enterococcaceae" = "#D1E5F0", "Erysipelotrichaceae" =  "#4EB265", "env.OPS 17" = "#F1932D", "Lineage IV" = "#E8601C", "Nitrosomonadaceae" = "#777777", "Obscuribacteraceae" = "#F6C141", "Polyangiaceae" = "lightgrey", "Pseudomonadaceae" = "#F1932D", "Sphingomonadaceae" = "#90C987", "Staphylococcaceae" = "#CAE0AB", "Xanthobacteraceae" = "#F7EE55", "Actinomycetaceae" = "#F6C141", "Streptococcaceae" = "#CAE0AB", "Moraxellaceae" = "lightgrey", "Enteroc. mundtii" = "#92C5DE", "Enteroc. cass/gall" = "#6195CF")  

# plot for frass before + after
ps3_pba_asvplot <- ggplot(ps3_pba_asv2, aes(x = sample4, y = Abundance, fill = family2, label=OTU)) + 
  geom_bar(stat = "identity", color="white") +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative abundance (ASV > 0.5%) \n") +
  theme_cowplot(11) + 
  facet_wrap(genotype ~ fct_relevel(transpl, "before", "transpl", "ctrl"), nrow=2, scales="free_x") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) 

ps3_pba_asvplot

pdf("Figures/ps3_pba_asvplot.pdf") 
print(ps3_pba_asvplot)
dev.off()

### supplementary taxonomy plots: whole larvae, diet, abdominal fluid, gut

# ASVs grouped by family but individual ASVs visible

# subsetting phyloseq object by tissue type
# whole larvae
mothps3_wL <- subset_samples(mothps3s2, tissue == "Whole_Larva")
mothps3_wL

# frass before and after (note that frass is called 'poo' in this script!)
mothps3_pba <- subset_samples(mothps3s2, tissue == "Poo_before" | tissue == "Poo_after")
mothps3_pba

# abdominal fluids
mothps3_af <- subset_samples(mothps3s2, tissue == "AF")
mothps3_af

# gut
mothps3_gut <- subset_samples(mothps3s2, tissue == "Gut")
mothps3_gut

# diet
mothps3_diet <- subset_samples(mothps3s2, tissue == "Control_Diet" | tissue == "WW_diet" | tissue == "yy_Diet")
mothps3_diet

# setting specific color to each family
fam_colors <- c("0319-6G20" = "#882E72", "Acetobacteraceae" = "#F6C141", "Anaerovoracaceae" = "#F6C141", "Bacillaceae" = "#F1932D", "Beijerinckiaceae" = "#B178A6", "Caulobacteraceae" = "darkgrey", "Chitinophagaceae" = "#D6C1DE", "Comamonadaceae" = "#1965B0", "Defluviicoccaceae" = "#777777", "Enterobacteriaceae" = "#5289C7", "Enterococcaceae" = "#7BAFDE", "Erysipelotrichaceae" =  "#4EB265", "env.OPS 17" = "#F1932D", "Lineage IV" = "#E8601C", "Nitrosomonadaceae" = "#777777", "Obscuribacteraceae" = "#F6C141", "Polyangiaceae" = "lightgrey", "Pseudomonadaceae" = "#F1932D", "Sphingomonadaceae" = "#90C987", "Staphylococcaceae" = "#CAE0AB", "Xanthobacteraceae" = "#F7EE55", "Actinomycetaceae" = "#F6C141", "Streptococcaceae" = "#CAE0AB", "Moraxellaceae" = "lightgrey")

# whole larvae: transforming data
ps3_wL_asv <- mothps3_wL %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.000) %>%                  # Filter out low abundance taxa
  arrange(Family) 

# plot for whole larvae
ps3_wL_asvplot <- ggplot(ps3_wL_asv, aes(x = sample4, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", color="white") +
  scale_fill_manual(values = fam_colors, limits=force) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative abundance (ASV > 0.5%) \n") +
  xlab("Genotype/family + sample ID") +
  ggtitle("Whole larva") +
  theme_cowplot(11) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  #theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(cols = vars(genotype), scales="free_x", space="free")

ps3_wL_asvplot

# abdominal fluids (AF): transforming data
ps3_af_asv <- mothps3_af %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.005) %>%                         # Filter out low abundance taxa
  arrange(Family) 

ps3_af_asv2 <- mutate(ps3_af_asv, family2 = ifelse(OTU %in% mund, "Enteroc. mundtii", ifelse(OTU %in% cass, "Enteroc. cass/gall", Family)))

# plot for AF
ps3_af_asvplot <- ggplot(ps3_af_asv2, aes(x = sample4, y = Abundance, fill = family2)) + 
  geom_bar(stat = "identity", color="white") +
  scale_fill_manual(values = fam_colors2, limits=force, name="Family") +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative abundance (ASV > 0.5%) \n") +
  xlab("Genotype/family + sample ID") +
  ggtitle("Abdominal fluid") +
  theme_cowplot(11) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_grid(.~ genotype + transpl, scales="free_x", space="free")

ps3_af_asvplot

# gut: transforming data
ps3_gut_asv <- mothps3_gut %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.005) %>%                         # Filter out low abundance taxa
  arrange(Family) 

ps3_gut_asv2 <- mutate(ps3_gut_asv, family2 = ifelse(OTU %in% mund, "Enteroc. mundtii", ifelse(OTU %in% cass, "Enteroc. cass/gall", Family)))

# plot for gut
ps3_gut_asvplot <- ggplot(ps3_gut_asv2, aes(x = sample4, y = Abundance, fill = family2, label=OTU)) + 
  geom_bar(stat = "identity", color="white") +
  scale_fill_manual(values = fam_colors2, limits=force, name="Family") +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative abundance (ASV > 0.5%) \n") +
  xlab("Genotype/family + sample ID") +
  ggtitle("Adult gut") +
  theme_cowplot(11) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_grid(.~ genotype + transpl, scales="free_x", space="free")

ps3_gut_asvplot

# diet: transforming data
ps3_diet_asv <- mothps3_diet %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.000) %>%                         # Filter out low abundance taxa
  arrange(Family) 

# plot for diet
ps3_diet_taxplot <- ggplot(ps3_diet_asv, aes(x = sample4, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", color="white") +
  scale_fill_manual(values = fam_colors, limits=force) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative abundance (ASV > 0.5%) \n") +
  xlab("Genotype/family + sample ID") +
  ggtitle("Diet") +
  theme_cowplot(11) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_grid(cols = vars(tissue), scales="free_x")

ps3_diet_taxplot
