
### Project: IBV intr-host
### Purpose: Plot IBV tree. Use a tree based on consensus sequences from sequence sample number 1 (first replicate) so there are no sample duplicates.

# =========== Load and import ==================

library(tidyverse)
library(ggtree)

meta_long <- read_csv("data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv")
tree.vic <- read.tree("data/processed/RAxML_bestTree.IBV_VIC_aln_tree.fa")
tree.yam <- read.tree("data/processed/RAxML_bestTree.IBV_YAM_aln_tree.fa")

# ======================== Modify tip labels ==========================

meta_long <- mutate(meta_long, tip_label = paste0(ENROLLID, "_", HOUSE_ID, "_", season, "_", pcr_result))
meta_long <- select(meta_long, tip_label, everything())

for(r in 1:nrow(meta_long))
{
  row <- meta_long[r,]
  sample_num <- row$SeqSampleNumber1
  new_tip <- row$tip_label
  tree.vic$tip.label[which(tree.vic$tip.label == sample_num)] <- new_tip
  tree.yam$tip.label[which(tree.yam$tip.label == sample_num)] <- new_tip
}

# ==================== Plot with ggtree ==========================

tree.vic.plot <- ggtree(tree.vic) + geom_treescale()
tree.vic.plot <- tree.vic.plot %<+% meta_long + geom_tiplab(aes(color = factor(season)), size = 0.5) # Can adjust size for refining in Illustrator

tree.yam.plot <- ggtree(tree.yam) + geom_treescale()
tree.yam.plot <- tree.yam.plot %<+% meta_long + geom_tiplab(aes(color = factor(season)), size = 0.5)

