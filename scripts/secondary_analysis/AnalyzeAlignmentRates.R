

### Project: IBV
### Purpose: Analyze alignment rate on original vs. new references.
### Working directory: Host_level_IBV_evolution

# ======================== Import packages, load data ==========================

library(tidyverse)

metadata <- read.csv("data/metadata/IBV_Sequenced_WithMetadata_InferredTypes.csv")
meta_snv <- read.csv("data/processed/meta_snv.csv")


# ================= Get original alignment rate from bowtie2 log files ===================

filenames_2445_VIC <- Sys.glob(paste("data/raw/alignment_logs/original/2445/Vic/*.log", sep = ""))
filenames_2445_YAM <- Sys.glob(paste("data/raw/alignment_logs/original/2445/Yam/*.log", sep = ""))
filenames_2446_VIC <- Sys.glob(paste("data/raw/alignment_logs/original/2446/Vic/*.log", sep = ""))
filenames_2446_YAM <- Sys.glob(paste("data/raw/alignment_logs/original/2446/Yam/*.log", sep = ""))

filenames_VIC <- data.frame("filename" = c(filenames_2445_VIC, filenames_2446_VIC))
filenames_YAM <- data.frame("filename" = c(filenames_2445_YAM, filenames_2446_YAM))
filenames_VIC %>% separate(filename, c("directory1", "directory2", "directory3", "directory4", "run", "lineage", "fourth"), sep = "/") %>% separate(fourth, c("seq_id", "filetype"), sep = "\\.") %>% mutate(filename = paste0(directory1, "/", directory2, "/", directory3, "/", directory4, "/", run, "/", lineage, "/", seq_id, ".", "log")) -> filenames_VIC
filenames_YAM %>% separate(filename, c("directory1", "directory2", "directory3", "directory4", "run", "lineage", "fourth"), sep = "/") %>% separate(fourth, c("seq_id", "filetype"), sep = "\\.") %>% mutate(filename = paste0(directory1, "/", directory2, "/", directory3, "/", directory4, "/", run, "/", lineage, "/", seq_id, ".", "log")) -> filenames_YAM

SeqIDs_VIC <- as.character(filter(metadata, final_result == "B_VIC")$SampleNumber)
SeqIDs_YAM <- as.character(filter(metadata, final_result == "B_YAM")$SampleNumber)
filenames_VIC <- filter(filenames_VIC, seq_id %in% SeqIDs_VIC)
filenames_YAM <- filter(filenames_YAM, seq_id %in% SeqIDs_YAM)
filenames_data <- rbind(filenames_VIC, filenames_YAM)

align_data <- data.frame(filename = filenames_data$filename)
align_data <- mutate(align_data, overall_align_rate = rep(NA,nrow(align_data)))
# go through align_data, read in file, assign alignment rate, merge with filenames_data. then merge with meta_snv.
for(r in 1:nrow(align_data))
{
  row <- align_data[r,]
  filename <- as.character(row$filename)
  
  info <- read.csv(filename)
  strsplit(as.character(info[14,]), split = " ")[[1]][1] %>% strsplit(split = "%") -> rate_info
  rate_info[[1]][1] %>% as.numeric() -> rate
  align_data$overall_align_rate[match(filename, align_data$filename)] <- rate
}

rate_data <- merge(filenames_data, align_data, by = "filename")
rate_data <- rename(rate_data, SeqSampleNumber1 = seq_id)
meta_snv_with_rate <- merge(meta_snv, rate_data, by = "SeqSampleNumber1")

# ==================================== Get alignment rates from new season-matched references ===============================

filenames_new_logs <- Sys.glob(paste("data/raw/alignment_logs/new/*.log", sep = ""))
align_data <- data.frame(filename = filenames_new_logs)
align_data <- mutate(align_data, overall_align_rate_new = rep(NA,nrow(align_data)))
for(r in 1:nrow(align_data))
{
  row <- align_data[r,]
  filename <- as.character(row$filename)
  
  info <- read.csv(filename)
  strsplit(as.character(info[14,]), split = " ")[[1]][1] %>% strsplit(split = "%") -> rate_info
  rate_info[[1]][1] %>% as.numeric() -> rate
  align_data$overall_align_rate_new[match(filename, align_data$filename)] <- rate
}

align_data %>% separate(filename, c("dir1", "dir2", "dir3", "dir4", "file"), sep = "/") %>% separate(file, c("SeqSampleNumber1", "filetype"), sep = "\\.") %>% select(-dir1, -dir2, -dir3, -dir4, -filetype) %>% mutate(SeqSampleNumber1 = as.integer(SeqSampleNumber1)) -> align_data
left_join(meta_snv_with_rate, align_data, by = "SeqSampleNumber1") %>% filter(!is.na(overall_align_rate_new)) -> meta_snv_with_rate_newref

align.rate.plot <- ggplot(meta_snv_with_rate_newref, aes(x = overall_align_rate, y = overall_align_rate_new)) + 
  geom_point() +
  xlab("Alignment rate to original influenza reference") +
  ylab("Alignment rate to new season-matched influenza reference") +
  geom_abline(intercept = 0, slope = 1,linetype = 2, size = 0.5, alpha = 0.5) +
  theme_bw() # PDF, 6 by 5

write.csv(meta_snv_with_rate_newref, "data/processed/meta_snv_withAlignRates.csv")

