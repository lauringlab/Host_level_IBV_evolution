
### Project: IBV intra-host
### Purpose: Rearrange metadata

# ============================= Import packages ================================

library(tidyverse)
library(magrittr)

# ============================= Metadata wrangling ================================

metadata_by_seq <- read_csv("data/metadata/IBV_Sequenced_WithMetadata_InferredTypes.csv")
samples <- data.frame(table(metadata_by_seq$uniqueID))
samples %>% separate(Var1, c("first", "second", "third")) -> samples_split
samples <- cbind(samples, samples_split)
samples <- subset(samples, select=which(!duplicated(names(samples))))
samples <- mutate(samples, ENROLLID = ifelse(str_detect(first, "HS"), second, first))
samples <- filter(samples, third != "<NA>") # remove plasmids
samples <- mutate(samples, is_home_sample = ifelse(str_detect(first, "HS"), TRUE, FALSE))
samples <- rename(samples, uniqueID = Var1)
samples <- rename(samples, seq_replicates = Freq)
samples <- mutate(samples, SPECID = ifelse(!is_home_sample, second, NA))
samples <- mutate(samples, home_spec_accn = ifelse(is_home_sample, first, NA))

# Rearranging the file of metadata by sequencing ID
#metadata_by_seq <- read_csv("../metadata/IBV_Sequenced_WithMetadata_InferredTypes.csv")
metadata_by_seq <- mutate(metadata_by_seq, home_collected = ifelse(str_detect(uniqueID, "HS"), TRUE, FALSE))
read_csv("data/processed/average_coverages.csv") %>% select(Id, coverage) %>% rename(SampleNumber = Id) -> cov_sample
metadata_by_seq <- merge(metadata_by_seq, cov_sample, by = "SampleNumber")
metadata_by_seq <- merge(metadata_by_seq, select(samples, uniqueID, seq_replicates), by = "uniqueID")
dups <- filter(metadata_by_seq, seq_replicates == 2)
single <- filter(metadata_by_seq, seq_replicates == 1)
single <- mutate(single, replicate = "singleton")
dups <- mutate(dups, replicate = ifelse(str_detect(SampleName, "rep_1"), "A", "B"))
metadata_by_seq <- rbind(dups, single)

# Also have each sample (clinic or home) on a separate line. "Long" version.
metadata_v4_long <- read_csv("data/metadata/flu_b_2010_2017_v4LONG.csv")

# Now put them together, with the long form as a base.
metadata_by_seq <- mutate(metadata_by_seq, ALV_ID = ifelse(!home_collected, spec_accn, home_spec_accn)) # not unique here
metadata_v4_long <- mutate(metadata_v4_long, ALV_ID = ifelse(home_spec == 0, SPECID, home_spec_accn)) # unique here
metadata_v4_long <- rename(metadata_v4_long, is_home_spec = home_spec)
metadata_v4_long <- mutate(metadata_v4_long, is_home_spec = ifelse(is_home_spec == 0, FALSE, TRUE))
metadata_by_seq[!duplicated(metadata_by_seq[,c("ALV_ID")]),] %>% select(ALV_ID, seq_replicates) -> unique_metadata_by_seq
metadata_v4_long_seqd <- merge(metadata_v4_long, unique_metadata_by_seq, by = "ALV_ID")

metadata_v4_long <- mutate(metadata_v4_long, sequenced = ifelse(ALV_ID %in% metadata_v4_long_seqd$ALV_ID, TRUE, FALSE))
write_csv(metadata_v4_long, "data/metadata/flu_b_2010_2017_v4LONG_seqdBinary.csv")

# Add sequencing numbers and coverage values.
metadata_v4_long_seqd <- mutate(metadata_v4_long_seqd, SeqSampleNumber1 = rep(NA,nrow(metadata_v4_long_seqd)))
metadata_v4_long_seqd <- mutate(metadata_v4_long_seqd, SeqSampleNumber2 = rep(NA,nrow(metadata_v4_long_seqd)))
metadata_v4_long_seqd <- mutate(metadata_v4_long_seqd, coverage1 = rep(NA,nrow(metadata_v4_long_seqd)))
metadata_v4_long_seqd <- mutate(metadata_v4_long_seqd, coverage2 = rep(NA,nrow(metadata_v4_long_seqd)))
for(r in 1:nrow(metadata_v4_long_seqd))
{
  row <- metadata_v4_long_seqd[r, ]
  alv_id <- row$ALV_ID
  reps <- row$seq_replicates
  
  by_seq <- filter(metadata_by_seq, ALV_ID == alv_id)
  
  if(reps == 1)
  {
    metadata_v4_long_seqd$SeqSampleNumber1[which(as.character(metadata_v4_long_seqd$ALV_ID) == as.character(alv_id))] = by_seq$SampleNumber
    metadata_v4_long_seqd$SeqSampleNumber2[which(as.character(metadata_v4_long_seqd$ALV_ID) == as.character(alv_id))] = NA
    metadata_v4_long_seqd$coverage1[which(as.character(metadata_v4_long_seqd$ALV_ID) == as.character(alv_id))] = by_seq$coverage
    metadata_v4_long_seqd$coverage2[which(as.character(metadata_v4_long_seqd$ALV_ID) == as.character(alv_id))] = NA
  }
  else if(reps == 2)
  {
    rowA <- filter(by_seq, replicate == "A")
    rowB <- filter(by_seq, replicate == "B")
    sA <- rowA$SampleNumber
    sB <- rowB$SampleNumber
    cA <- rowA$coverage
    cB <- rowB$coverage
    metadata_v4_long_seqd$SeqSampleNumber1[which(as.character(metadata_v4_long_seqd$ALV_ID) == as.character(alv_id))] = sA
    metadata_v4_long_seqd$SeqSampleNumber2[which(as.character(metadata_v4_long_seqd$ALV_ID) == as.character(alv_id))] = sB
    metadata_v4_long_seqd$coverage1[which(as.character(metadata_v4_long_seqd$ALV_ID) == as.character(alv_id))] = cA
    metadata_v4_long_seqd$coverage2[which(as.character(metadata_v4_long_seqd$ALV_ID) == as.character(alv_id))] = cB
  }
  
}

metadata_v4_long_seqd <- mutate(metadata_v4_long_seqd, coverage_qualified = ifelse(seq_replicates == 1, ifelse(coverage1 > 1000, TRUE, FALSE), ifelse(coverage1 > 1000 & coverage2 > 1000, TRUE, FALSE)))
metadata_v4_long_seqd$pcr_result[metadata_v4_long_seqd$pcr_result == "B YAM"] <- "B/YAM"
metadata_v4_long_seqd$pcr_result[metadata_v4_long_seqd$pcr_result == "B"] <- "B/YAM" # the two untyped were Yam, via consensus phylogenetic trees
metadata_v4_long_seqd$pcr_result[metadata_v4_long_seqd$pcr_result == "B/Yam"] <- "B/YAM"
metadata_v4_long_seqd$pcr_result[metadata_v4_long_seqd$pcr_result == "B VIC"] <- "B/VIC"
metadata_v4_long_seqd$pcr_result[metadata_v4_long_seqd$pcr_result == "B/Vic"] <- "B/VIC"
metadata_v4_long_seqd <- mutate(metadata_v4_long_seqd, DPSO = ifelse(!is_home_spec, collect-onset, home_spec_date-onset)) # DPSO = days post symptom onset

write_csv(metadata_v4_long_seqd, "data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo.csv")
write_csv(metadata_by_seq, "data/metadata/Metadata_By_Seq_withALVID.csv")

# Got copy number info. Placed in manually, saved as flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv

newmeta <- read_csv("data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv")
newmeta <- rename(newmeta, genome_copy_per_ul = genome_copy_num)
newmeta <- mutate(newmeta, ln_copy_num = log(genome_copy_per_ul))
write_csv(newmeta, "data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv")

