
### Project: IBV intr-host
### Purpose: Plot IBV tree. Use a tree based on consensus sequences from sequence sample number 1 (first replicate) so there are no sample duplicates.

# =========== Load and import ==================

library(tidyverse)
library(wesanderson)
library(magrittr)
library(ggtree)

setwd("/Users/avalesano/Documents/MSTP/LauringLab/Host_level_IBV_evolution/scripts/")
#meta_pairs <- read_csv("../data/processed/household_pairs_wSeqID.csv")
meta_long <- read_csv("../data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv")
tree <- read.tree("../data/processed/RAxML_bestTree.IBV_raxml_PairAnalysis.tree")

# ======================== Some metadata wrangling ==========================

#meta_pairs <- mutate(meta_pairs, pcr_result = ifelse(str_detect(pcr_result, "Y"), "B/YAM", "B/VIC"))
meta_long <- mutate(meta_long, tip_label = paste0(ENROLLID, "_", HOUSE_ID, "_", season, "_", pcr_result))
meta_long <- select(meta_long, tip_label, everything())

#meta_long <- mutate(meta_long, pairs_house_id = ifelse(SeqSampleNumber1 %in% meta_pairs$SeqSampleNumber1, HOUSE_ID, NA)) # want only the ones in the household dataset

# Change tip names in tree
for(r in 1:nrow(meta_long))
{
  row <- meta_long[r,]
  sample_num <- row$SeqSampleNumber1
  new_tip <- row$tip_label
  tree$tip.label[which(tree$tip.label == sample_num)] <- new_tip
}

# ==================== ggtree ==========================

palette <- wesanderson::wes_palette("FantasticFox1", 7, type = "continuous")

raxml_tree <- ggtree(tree) + geom_treescale()
#raxml_tree <- raxml_tree %<+% meta_long + geom_tiplab(aes(color = factor(pairs_house_id)), size = 3) + scale_color_manual(values = palette)
raxml_tree <- raxml_tree %<+% meta_long + geom_tiplab(aes(color = factor(season)), size = 2) #+ scale_color_manual(values = palette)

raxml_tree <- raxml_tree + theme(legend.position = "bottom", legend.text = element_text(size = 11)) + labs(color = "House ID")
raxml_tree <- raxml_tree + guides(colour = guide_legend(override.aes = list(size = 7, shape = 16, alpha = 1)))
raxml_tree

ggsave(plot = raxml_tree, filename = "../results/plots/Tree_RAxML_BySeason.jpg", device = "jpeg", width = 10)
