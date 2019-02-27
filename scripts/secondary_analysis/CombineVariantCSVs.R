
### Project: IBV intra-host
### Purpose: combine the correct variant data from each run into one CSV.

# ============================ Packages, read in metadata =======================================

library(tidyverse)
setwd("/Users/avalesano/Documents/MSTP/LauringLab/Host_level_IBV_evolution/scripts/")

metadata <- read.csv("../data/metadata/IBV_Sequenced_WithMetadata_InferredTypes.csv")

# ========================= Read in and filter the CSVs by samples from the appropriate lineage ============================

vic.2445 <- read_csv("../data/raw/2445_VIC/all.variants.csv")
yam.2445 <- read_csv("../data/raw/2445_YAM/all.variants.csv")
vic.2446 <- read_csv("../data/raw/2446_VIC/all.variants.csv")
yam.2446 <- read_csv("../data/raw/2446_YAM/all.variants.csv")

metadata_VIC <- filter(metadata, final_result == "B_VIC")
metadata_YAM <- filter(metadata, final_result == "B_YAM")

vic.2445 <- filter(vic.2445, as.character(Id) %in% metadata_VIC$SampleNumber)
yam.2445 <- filter(yam.2445, as.character(Id) %in% metadata_YAM$SampleNumber)
vic.2446 <- filter(vic.2446, as.character(Id) %in% metadata_VIC$SampleNumber)
yam.2446 <- filter(yam.2446, as.character(Id) %in% metadata_YAM$SampleNumber)

### Combine CSVs by lineage for a full VIC/YAM variant CSV.

vic.all <- rbind(vic.2445, vic.2446)
yam.all <- rbind(yam.2445, yam.2446)
all.variants <- rbind(vic.all, yam.all)

# Get chr and mutation names all the same; some are NA_, others NA_clone1
all.variants$chr <- gsub("NA_clone1", "NR", all.variants$chr)
all.variants$chr <- gsub("NA_", "NR", all.variants$chr)
all.variants$mutation <- gsub("NA_clone1", "NR", all.variants$mutation)
all.variants$mutation <- gsub("NA_", "NR", all.variants$mutation)

write_csv(all.variants, "../data/processed/all.variants.csv")
