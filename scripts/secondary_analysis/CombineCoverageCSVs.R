
### Project: IBV intra-host
### Purpose: Combine the correct variant data from each run into one CSV.

# ======================== Import packages ================================

library(tidyverse)
setwd("/Users/avalesano/Documents/MSTP/LauringLab/Host_level_IBV_evolution/scripts/")

metadata <- read.csv("../data/metadata/IBV_Sequenced_WithMetadata_InferredTypes.csv")

# ======================== Filter the CSVs by samples from the appropriate lineage ================================

vic.2445 <- read_csv("../data/raw/2445_VIC/all.coverage.csv")
yam.2445 <- read_csv("../data/raw/2445_YAM/all.coverage.csv")
vic.2446 <- read_csv("../data/raw/2446_VIC/all.coverage.csv")
yam.2446 <- read_csv("../data/raw/2446_YAM/all.coverage.csv")

metadata_VIC <- filter(metadata, final_result == "B_VIC")
metadata_YAM <- filter(metadata, final_result == "B_YAM")

vic.2445 <- filter(vic.2445, as.character(Id) %in% metadata_VIC$SampleNumber)
yam.2445 <- filter(yam.2445, as.character(Id) %in% metadata_YAM$SampleNumber)
vic.2446 <- filter(vic.2446, as.character(Id) %in% metadata_VIC$SampleNumber)
yam.2446 <- filter(yam.2446, as.character(Id) %in% metadata_YAM$SampleNumber)

# ======================== Combine for one single coverage CSV ================================

vic.all <- rbind(vic.2445, vic.2446)
yam.all <- rbind(yam.2445, yam.2446)

setdiff(unique(vic.all$Id), metadata_VIC$SampleNumber) # should be 0
setdiff(unique(vic.all$Id), metadata_YAM$SampleNumber)
setdiff(unique(yam.all$Id), metadata_YAM$SampleNumber) # should be 0
setdiff(unique(yam.all$Id), metadata_VIC$SampleNumber)

all.coverage <- rbind(vic.all, yam.all)
write_csv(all.coverage, "../data/processed/all.coverage.csv")
