
### Project: IBV
### Purpose: Analyze variants on original vs. new references.
### Working directory: Host_level_IBV_evolution

# ======================== Import packages, load data ==========================

library(tidyverse)

meta_snv <- read.csv("data/processed/meta_snv.csv")
yam1415 <- read.csv("data/processed/yam1415.minor.quality.snv.csv")
vic1516 <- read.csv("data/processed/vic1516.minor.quality.snv.csv")
yam1617 <- read.csv("data/processed/yam1617.minor.quality.snv.csv")
new_calls <- rbind(yam1415, vic1516, yam1617)

read.csv("data/processed/qual.snv.csv") %>% 
  filter(freq.var < 0.5) %>%
  filter((pcr_result == "B/YAM" & season == "2014-2015") | (pcr_result == "B/VIC" & season == "2015-2016") | (pcr_result == "B/YAM" & season == "2016-2017")) -> old_calls

read_csv("data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv") %>%
  filter((pcr_result == "B/YAM" & season == "2014-2015") | (pcr_result == "B/VIC" & season == "2015-2016") | (pcr_result == "B/YAM" & season == "2016-2017")) %>%
  filter(coverage_qualified == TRUE) -> meta

old_calls %>% group_by(ALV_ID) %>% dplyr::summarize(iSNV = length(unique(mutation))) -> isnv_minority_count
meta_snv <- merge(meta, isnv_minority_count, by = "ALV_ID", all.x = TRUE)
meta_snv$iSNV[is.na(meta_snv$iSNV)] <- 0 # NA means there were no variants; change to 0

new_calls %>% group_by(ALV_ID) %>% dplyr::summarize(iSNV_new = length(unique(mutation))) -> isnv_minority_count_new
meta_snv <- merge(meta_snv, isnv_minority_count_new, by = "ALV_ID", all.x = TRUE)
meta_snv$iSNV_new[is.na(meta_snv$iSNV_new)] <- 0 # NA means there were no variants; change to 0

old_snv <- select(meta_snv, iSNV)
old_snv <- mutate(old_snv, type = "Original reference")
new_snv <- select(meta_snv, iSNV_new)
new_snv <- mutate(new_snv, type = "Season-matched reference")
new_snv <- rename(new_snv, iSNV = iSNV_new)
compare_df <- rbind(old_snv, new_snv)

snv.by.ref.plot <- ggplot(compare_df, aes(y = iSNV, x = as.factor(type))) +
  geom_dotplot(stackdir = "center", binaxis = 'y', binwidth = 1, dotsize = 0.05) + 
  xlab("") + 
  ylab("iSNV Per Sample") + 
  theme_bw() + 
  theme(legend.position = "") # save as PDF, 6 by 5

