

### Project: IBV
### Purpose: Generate figure panels for most of the figures. Much of this can be found in `SecondaryAnalysisExploration.R`.
### Working directory: Host_level_IBV_evolution

# =========== Import packages, load data ================

library(tidyverse)
library(wesanderson)
library(ggbeeswarm)
palette <- wesanderson::wes_palette("FantasticFox1")

meta <- read_csv("data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv")
meta_short <- read_csv("data/metadata/flu_b_2010_2017_v4.csv")
metadata_by_seq <- read_csv("data/metadata/Metadata_By_Seq_withALVID.csv")
quality_var <- read_csv("data/processed/qual.snv.csv")
no_freq_cut_var <- read_csv("data/processed/no_freq_cut.qual.snv.csv")
qual_not_coll <- read_csv("data/processed/qual.not.collapsed.snv.csv")
chrs <- read.csv("data/metadata/ibv.segs.csv", stringsAsFactors = TRUE)
IAV_snv_qual_meta <- read_csv("data/processed/IAV_snv_qual_meta.csv")
tree.vic <- read.tree("data/processed/RAxML_bestTree.IBV_VIC_aln_tree.fa")
tree.yam <- read.tree("data/processed/RAxML_bestTree.IBV_YAM_aln_tree.fa")
possible_pairs.dist <- read.csv("data/processed/possible.pairs.dist.csv")
trans_freq <- read.csv("data/processed/trans_freq.csv")
meta_snv <- read.csv("data/processed/meta_snv.csv")

isnv <- filter(quality_var, freq.var < 0.98)
isnv_nomixed <- filter(isnv, !(ENROLLID %in% c("50425")))
isnv_minority <- filter(quality_var, freq.var < 0.5)
isnv_minority_nomixed <- filter(isnv_minority, !(ENROLLID %in% c("50425"))) # remove mixed infections

# ==================== Figure 1 ===================

# Figure 1A
titer.plot <- ggplot(meta, aes(x = as.factor(DPSO), y = genome_copy_per_ul)) + 
  geom_boxplot(notch = TRUE) + 
  ylab(expression(paste(Genomes, "/", mu, L))) + 
  scale_y_log10() + 
  xlab("Days Post Symptom Onset") + 
  theme_bw()

# Note: Figure 1B can be found in `PlotCoverage.R`.
  
# ==================== Figure 2 ===================

# Figure 2A
isnv_by_day.plot <- ggplot(meta_snv, aes(x = as.factor(DPSO), y = iSNV)) + 
  geom_boxplot(outlier.shape = NA, notch = FALSE) + 
  xlab("Day Post Symptom Onset") + 
  geom_jitter(width = 0.3, height = 0.1, size = 1.5) + 
  theme_bw() + 
  ylab("iSNV Per Sample") # save as PDF, 5 by 4

# Figure 2B
snv_by_copynum <- ggplot(meta_snv, aes(x = log(genome_copy_per_ul, 10), y = iSNV)) + 
  geom_point(shape = 19, size = 1.5) + 
  ylab("iSNV Per Sample") + 
  theme_bw() +
  xlab(expression(paste("Log (base 10) of genomes/", mu, L))) # save as PDF, 5 by 4

# Figure 2C
isnv_not_collapsed <- filter(qual_not_coll, freq.var < 0.98)
isnv_not_collapsed_rep1 <- filter(isnv_not_collapsed, SampleNumber == SeqSampleNumber1)
isnv_not_collapsed_rep2 <- filter(isnv_not_collapsed, SampleNumber == SeqSampleNumber2)
isnv_not_collapsed_rep1 <- mutate(isnv_not_collapsed_rep1, mut_id = paste0(mutation, "_", ALV_ID))
isnv_not_collapsed_rep2 <- mutate(isnv_not_collapsed_rep2, mut_id = paste0(mutation, "_", ALV_ID))
merged <- merge(isnv_not_collapsed_rep1, isnv_not_collapsed_rep2, by = "mut_id")
merged_nomixed <- filter(merged, !(ALV_ID.x %in% c("HS1875", "MH10536")))

replicate_concordance_no_mixed <- ggplot(data = merged_nomixed, aes(x = freq.var.x, y = freq.var.y)) +
  geom_point(aes(color = as.factor(floor(log10(genome_copy_per_ul.x))))) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  scale_y_continuous(limits = c(0,0.5)) +
  scale_x_continuous(limits = c(0,0.5)) +
  xlab("Frequency in replicate 1") + 
  ylab("Frequency in replicate 2") +
  scale_color_manual(name = "Log(copies/ul)", values = palette[c(4,3)]) + 
  theme_bw() + 
  theme(legend.position = "")

# ==================== Figure 3 ===================

palette = wesanderson::wes_palette("FantasticFox1")
positions <- read_csv("data/processed/FluSegmentPositions.csv") # segment lengths from reference files
positions_YAM <- filter(positions, Lineage == "B/YAM")

isnv_minority_nomixed_concatpos <- mutate(isnv_minority_nomixed, concat_pos = positions_YAM$AddLength[match(chr, positions_YAM$Segment)] + pos)

freq.by.pos <- ggplot(isnv_minority_nomixed_concatpos, aes(x = concat_pos, y = freq.var, fill = class_factor)) +
  geom_point(size = 3, shape = 21) +
  scale_fill_manual(values = palette[c(4,3)]) +
  theme_bw() + 
  xlab("Concatenated Genome Position") + 
  ylab("Frequency") +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  scale_x_continuous(labels = positions_YAM$Segment, breaks = positions_YAM$AddLength) # save as PDF, 10 by 4.


# ==================== Figure 4 ===================

# Figure 4A
isnv_by_vaccination <- ggplot(meta_snv, aes(y = iSNV, x = as.factor(vaccination_status))) +
  geom_dotplot(stackdir = "center", binaxis = 'y', binwidth = 1, dotsize = 0.25) +
  scale_x_discrete(labels = c("Not Vaccinated", "Vaccinated")) + 
  xlab("") + 
  theme_bw() + 
  ylab("iSNV Per Sample") +
  theme(legend.position = "") # PDF, 5 by 4

# Figure 4B
Vic_Yam_plot <- ggplot(meta_snv, aes(y = iSNV, x = as.factor(pcr_result))) +
  geom_dotplot(stackdir = "center", binaxis = 'y', binwidth = 1, dotsize = 0.3) +
  xlab("") + 
  ylab("iSNV Per Sample") + 
  theme_bw() + 
  theme(legend.position = "") # PDF, 5 by 4

# Figure 4C
IBV_snv <- select(meta_snv, iSNV)
IBV_snv <- mutate(IBV_snv, type = "Influenza B")
IAV_snv <- select(IAV_snv_qual_meta, iSNV)
IAV_snv <- mutate(IAV_snv, type = "Influenza A")
cfIAV_data <- rbind(IBV_snv, IAV_snv)

cfIAV.dotplot <- ggplot(cfIAV_data, aes(y = iSNV, x = as.factor(type), fill = type, color = type)) +
  geom_dotplot(stackdir = "center", binaxis = 'y', binwidth = 1, dotsize = 0.4) + 
  scale_fill_manual(values = c("black", palette[5])) + 
  scale_color_manual(values = c("black", palette[5])) +
  xlab("") + 
  ylab("iSNV Per Sample") + 
  theme_bw() + 
  theme(legend.position = "")

# ==================== Figure 5 ===================

# Figure 5A and 5B

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

tree.vic.plot <- ggtree(tree.vic) + geom_treescale()
tree.vic.plot <- tree.vic.plot %<+% meta_long + geom_tiplab(aes(color = factor(season)), size = 0.5) # Can adjust size for refining in Illustrator

tree.yam.plot <- ggtree(tree.yam) + geom_treescale()
tree.yam.plot <- tree.yam.plot %<+% meta_long + geom_tiplab(aes(color = factor(season)), size = 0.5)

# Figure 5C

cutoffs <- tibble(L1_norm = quantile(possible_pairs.dist$L1_norm[possible_pairs.dist$Household == FALSE], probs = seq(0, 1, 0.05)))
cutoffs$threshold <- seq(0, 1, 0.05)
cutoffs %>% 
  rowwise() %>% 
  mutate(valid_pairs = nrow(possible_pairs.dist[(possible_pairs.dist$valid == TRUE & possible_pairs.dist$L1_norm < L1_norm),])) -> cutoffs
cutoff <- cutoffs$L1_norm[cutoffs$threshold == 0.05]

palette = wesanderson::wes_palette("FantasticFox1")
possible_pairs.dist %>% filter(valid == TRUE | Household == FALSE) %>% select(season, ALV_ID_1, ALV_ID_2, L1_norm, valid, Household) -> L1norm_plot_data

L1norm_plot <- ggplot(L1norm_plot_data, aes(x = L1_norm, fill = as.factor((valid-1)*-1), y = ..ncount..)) +
  geom_histogram(binwidth = 7.5, boundary = 0, position = 'dodge') +
  scale_fill_manual(name = "", labels = c("Household pair", "Community pair"), values = palette[c(3,4)]) +
  xlab("L1 Norm") + 
  ylab("Normalized Count") +
  theme(legend.position = c(0.5, 0.5)) +
  geom_segment(aes(x = cutoff, xend = cutoff, y = 0, yend = 1), linetype = 2, color = palette[5], size = 0.4) + 
  theme_classic() +
  theme(text = element_text(size = 28), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20)) # save as PDF, 15 by 10


# ==================== Figure 6 ===================

snv.betweenhost.plot <- ggplot(trans_freq, aes(x = freq1, y = freq2)) + 
  geom_point(size = 1.5) + 
  xlab("Frequency in donor") + 
  ylab("Frequency in recipient") + 
  theme_classic() # save as PDF, 5 by 4
  
