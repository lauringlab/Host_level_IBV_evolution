

### Purpose: Make a table linking the sequencing core IDs and the master metadata.
### Here, need to add the metadata information to each line of the sequencing sheet, so that there is one line per sequencing sample number.

# ========================== Load in data =========================

library(tidyverse)

MasterMetadata <- read_csv("data/metadata/flu_b_2010_2017_v3.csv")
SequencingIDs <- read_csv("data/metadata/SeqSubmission.csv")

# Re-label misnamed replicates
SequencingIDs$SampleName[which("301665_M53896_2010_rep_1" == SequencingIDs$SampleName & SequencingIDs$PoolNumber == 114367)] <- "301665_M53896_2010_rep_2" # these five were mis-labeled as rep1, when really they were rep2 on the second run.
SequencingIDs$SampleName[which("320413_MH1033_2011_rep_1" == SequencingIDs$SampleName & SequencingIDs$PoolNumber == 114367)] <- "320413_MH1033_2011_rep_2"
SequencingIDs$SampleName[which("330007_MH2666_2012_rep_1" == SequencingIDs$SampleName & SequencingIDs$PoolNumber == 114367)] <- "330007_MH2666_2012_rep_2"
SequencingIDs$SampleName[which("330108_MH4301_2012_rep_1" == SequencingIDs$SampleName & SequencingIDs$PoolNumber == 114367)] <- "330108_MH4301_2012_rep_2"
SequencingIDs$SampleName[which("331397_MH4247_2012_rep_1" == SequencingIDs$SampleName & SequencingIDs$PoolNumber == 114367)] <- "331397_MH4247_2012_rep_2"

# ========================

# Make empty columns to fill in the loop
SequencingIDs <- mutate(SequencingIDs, SEASON_STUDYID = rep(NA,nrow(SequencingIDs)))
SequencingIDs <- mutate(SequencingIDs, SEASON_HOUSEID = rep(NA,nrow(SequencingIDs)))
SequencingIDs <- mutate(SequencingIDs, MASTER_ID = rep(NA,nrow(SequencingIDs)))
SequencingIDs <- mutate(SequencingIDs, MASTER_HOUSE_ID = rep(NA,nrow(SequencingIDs)))
SequencingIDs <- mutate(SequencingIDs, season = rep(NA,nrow(SequencingIDs)))
SequencingIDs <- mutate(SequencingIDs, onset_date = rep(NA,nrow(SequencingIDs)))
SequencingIDs <- mutate(SequencingIDs, collect_date = rep(NA,nrow(SequencingIDs)))
SequencingIDs <- mutate(SequencingIDs, spec_accn = rep(NA,nrow(SequencingIDs)))
SequencingIDs <- mutate(SequencingIDs, home_spec = rep(NA,nrow(SequencingIDs)))
SequencingIDs <- mutate(SequencingIDs, home_spec_date = rep(NA,nrow(SequencingIDs)))
SequencingIDs <- mutate(SequencingIDs, home_spec_accn = rep(NA,nrow(SequencingIDs)))
SequencingIDs <- mutate(SequencingIDs, final_result = rep(NA,nrow(SequencingIDs)))
SequencingIDs <- mutate(SequencingIDs, flub_count = rep(NA,nrow(SequencingIDs)))
SequencingIDs <- mutate(SequencingIDs, cluster = rep(NA,nrow(SequencingIDs)))
SequencingIDs <- mutate(SequencingIDs, uniqueID = rep(NA,nrow(SequencingIDs)))


for(sample_name in SequencingIDs$SampleName) # these are unique
{
  name_components <- strsplit(sample_name, split = "_")[[1]]
  id_first <- name_components[1]
  id_second <- name_components[2]
  rep_num <- as.character(name_components[5])
  uniqueID <- paste(name_components[1], "_", name_components[2], "_", name_components[3], sep = "")
  
  spec_accn1 <- NA
  season_study_id <- NA
  home_spec_accn1 <- NA
  
  if(str_detect(id_first, "H"))
  {
    home_spec_accn1 <- as.character(id_first)
    season_study_id <- as.character(id_second)
    
  }
  else if(grepl("YAM", id_first))
  {
    next # it's a plasmid control
  }
  else if(grepl("VIC", id_first))
  {
    next # it's a plasmid control
  }
  else
  {
    spec_accn1 <- as.character(id_second)
    season_study_id <- as.character(id_first)
  }
  
  # now we have spec id; either spec_accn or home_spec_accn. Also have season study id.
  # want to find the spec_accn and assign based on that.
  #print(paste(season_study_id, spec_accn1, home_spec_accn1))
  
  metadata_row <- NA
  if(!is.na(spec_accn1))
  {
    metadata_row <- filter(MasterMetadata, SEASON_STUDYID == season_study_id & spec_accn == spec_accn1)
    SequencingIDs[SequencingIDs$SampleName == sample_name, "SEASON_STUDYID" ] <- metadata_row$SEASON_STUDYID
    SequencingIDs[SequencingIDs$SampleName == sample_name, "SEASON_HOUSEID" ] <- metadata_row$SEASON_HOUSEID
    SequencingIDs[SequencingIDs$SampleName == sample_name, "MASTER_ID" ] <- metadata_row$MASTER_ID
    SequencingIDs[SequencingIDs$SampleName == sample_name, "MASTER_HOUSE_ID" ] <- metadata_row$MASTER_HOUSE_ID
    SequencingIDs[SequencingIDs$SampleName == sample_name, "season" ] <- metadata_row$season
    SequencingIDs[SequencingIDs$SampleName == sample_name, "onset_date" ] <- metadata_row$onset_date
    SequencingIDs[SequencingIDs$SampleName == sample_name, "collect_date" ] <- metadata_row$collect_date
    SequencingIDs[SequencingIDs$SampleName == sample_name, "spec_accn" ] <- metadata_row$spec_accn
    SequencingIDs[SequencingIDs$SampleName == sample_name, "home_spec" ] <- metadata_row$home_spec
    SequencingIDs[SequencingIDs$SampleName == sample_name, "home_spec_date" ] <- metadata_row$home_spec_date
    SequencingIDs[SequencingIDs$SampleName == sample_name, "home_spec_accn" ] <- metadata_row$home_spec_accn
    SequencingIDs[SequencingIDs$SampleName == sample_name, "final_result" ] <- metadata_row$final_result
    SequencingIDs[SequencingIDs$SampleName == sample_name, "flub_count" ] <- metadata_row$flub_count
    SequencingIDs[SequencingIDs$SampleName == sample_name, "cluster" ] <- metadata_row$cluster
    SequencingIDs[SequencingIDs$SampleName == sample_name, "uniqueID" ] <- uniqueID
  }
  else if(!is.na(home_spec_accn1))
  {
    metadata_row <- filter(MasterMetadata, SEASON_STUDYID == season_study_id & home_spec_accn == home_spec_accn1)
    SequencingIDs[SequencingIDs$SampleName == sample_name, "SEASON_STUDYID" ] <- metadata_row$SEASON_STUDYID
    SequencingIDs[SequencingIDs$SampleName == sample_name, "SEASON_HOUSEID" ] <- metadata_row$SEASON_HOUSEID
    SequencingIDs[SequencingIDs$SampleName == sample_name, "MASTER_ID" ] <- metadata_row$MASTER_ID
    SequencingIDs[SequencingIDs$SampleName == sample_name, "MASTER_HOUSE_ID" ] <- metadata_row$MASTER_HOUSE_ID
    SequencingIDs[SequencingIDs$SampleName == sample_name, "season" ] <- metadata_row$season
    SequencingIDs[SequencingIDs$SampleName == sample_name, "onset_date" ] <- metadata_row$onset_date
    SequencingIDs[SequencingIDs$SampleName == sample_name, "collect_date" ] <- metadata_row$collect_date
    SequencingIDs[SequencingIDs$SampleName == sample_name, "spec_accn" ] <- metadata_row$spec_accn
    SequencingIDs[SequencingIDs$SampleName == sample_name, "home_spec" ] <- metadata_row$home_spec
    SequencingIDs[SequencingIDs$SampleName == sample_name, "home_spec_date" ] <- metadata_row$home_spec_date
    SequencingIDs[SequencingIDs$SampleName == sample_name, "home_spec_accn" ] <- metadata_row$home_spec_accn
    SequencingIDs[SequencingIDs$SampleName == sample_name, "final_result" ] <- metadata_row$final_result
    SequencingIDs[SequencingIDs$SampleName == sample_name, "flub_count" ] <- metadata_row$flub_count
    SequencingIDs[SequencingIDs$SampleName == sample_name, "cluster" ] <- metadata_row$cluster
    SequencingIDs[SequencingIDs$SampleName == sample_name, "uniqueID" ] <- uniqueID
  }
}

# Standardize the lineage naming
SequencingIDs$final_result[SequencingIDs$final_result == "B YAM"] <- "B/YAM"
SequencingIDs$final_result[SequencingIDs$final_result == "B/Yam"] <- "B/YAM"
SequencingIDs$final_result[SequencingIDs$final_result == "B VIC"] <- "B/VIC"
SequencingIDs$final_result[SequencingIDs$final_result == "B/Vic"] <- "B/VIC"

# Two untyped samples, just labeled "B", were found to cluster with YAM sequences on a phylogeny. Add that information here.
SequencingIDs$final_result[SequencingIDs$final_result == "B"] <- "B/YAM"

write_csv(SequencingIDs, "data/metadata/IBV_Sequenced_WithMetadata_InferredTypes_v2.csv")
