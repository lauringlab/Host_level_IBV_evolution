
### Project: IBV intra-host
### Purpose: Validate putative transmission pairs by L1-norm genetic distance.

# ================== Read in data, import packages =========================

library(tidyverse)

meta_long <- read_csv("data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv") # Each row is a biological sample
qual <- read_csv("data/processed/qual.snv.csv")
meta_wide <- read_csv("data/metadata/flu_b_2010_2017_v4.csv") # Each row is an individual

# ================== Metadata wrangling ===================================

# Dates in meta_wide are YYYY-MM-DD. Change format of onset in meta_long to match.
meta_long <- mutate(meta_long, onset = as.Date(onset, "%m/%d/%y"))

meta_wide <- mutate(meta_wide, home_onset = paste0(home_spec_accn, "_", onset))
meta_wide <- mutate(meta_wide, clinic_onset = paste0(SPECID, "_", onset))

meta_long %>% filter(coverage_qualified == TRUE) %>% mutate(ID_onset = paste0(ALV_ID, "_", onset)) -> meta_long_qualitySeq # all sequenced and have >1000x average coverage

meta_wide <- mutate(meta_wide, haveQualityHomeSample = ifelse(home_onset %in% meta_long_qualitySeq$ID_onset, TRUE, FALSE))
meta_wide <- mutate(meta_wide, haveQualityClinicSample = ifelse(clinic_onset %in% meta_long_qualitySeq$ID_onset, TRUE, FALSE)) 

meta_wide_seqd <- filter(meta_wide, haveQualityHomeSample == TRUE | haveQualityClinicSample == TRUE) # Have a home and/or clinic sample with quality sequence
meta_wide_seqd <- mutate(meta_wide_seqd, ID_onset = paste0(ENROLLID, "_", onset)) # these are individual infections for which we have some sequence data
#write_csv(meta_wide_seqd, "data/metadata/metadata_wide_quality_sequence.csv") # Write to intermediate file for later use.

# ================== Functions for L1-norm =========================

equal_compare <- function(position)
{ 
  # Take in all the variants found in a sample pair at a given posistion and correct for differences in inference. 
  # Meaning the reference base was called here but only in one sample because in the other sample there wasn't a minor allele. 
  # Each position with a minor variant should have at least 2 variants, in both samples after this is done. A minor variant and the reference. 
  # If a position has no variants in either sample then it is excluded from this analysis. In this case only the reference base was present in both samples.  
  # If only major variants are present then the major allele is fixed and there are no other minor variants to infer.
  
  pos_sum_1 <- sum(position$freq1) # What is the sum of all variants at this site in the first sample
  pos_sum_2 <- sum(position$freq2) # What is the sum of all variants at this site in the second sample
  stopifnot(pos_sum_1 + pos_sum_2 > 0) # No variant at this position - something terrible happend somewhere
  if((pos_sum_1 > 0 & min(position$freq1) < 0.95 & pos_sum_1 == min(position$freq1))) warning(paste0("There is only 1 allele in sample", position$ALV_ID1, ", its frequency is : ", position$freq1), ", it may be removed.\n")
  if((pos_sum_2 > 0 & min(position$freq2) < 0.95 & pos_sum_2 == min(position$freq2))) warning(paste0("There is only 1 allele in sample", position$ALV_ID2, ", its frequency is : ", position$freq2), ", it may be removed.\n") # only a minor variant called
  
  if(pos_sum_1 == 0)
  { # there isn't a variant in the first sample. The reference is fixed there. the sum of nothing is 0. Try it.  sum() =0.  sum(c()) =0
    
    x <- which(position$ref == position$var) # Where are the infered calls
    stopifnot(length(x) < 2) # should be 1 or 0.
    if(length(x) == 1)
    { # The reference has been infered in the second sample but there are no varaints in the first so it was not infered
      position$freq1[x] <- 1 # The reference base is fixed here
    } else if(length(x) == 0)
    { # the reference was not infered in the second sample and no variants here. Add the reference now. There must have been another base fixed in the second sample.
      extra <- dplyr::mutate(position, var = ref, mutation = paste0(chr, "_", ref, pos, var), freq1 = 1, freq2 = 0) # adding a line with the reference base fixed in sample 1 but not sample 2
      position <- rbind(position, extra)
    }
    # Do it again if there are no variants in the second sample.
  } else if(pos_sum_2 == 0)
  { # there isn't a variant in the second sample. The reference is fixed.
    x <- which(position$ref == position$var)
    stopifnot(length(x) < 2) # should be 1 or 0
    if(length(x) == 1)
    { # The reference has been inferred in the first sample but there are no variants in the second so it was not inferred
      position$freq2[x] <- 1
    } else if(length(x) == 0)
    { # the reference was not inferred in the first sample and no variants here. Add the reference now
      extra <- dplyr::mutate(position, var = ref, mutation = paste0(chr, "_", ref, pos, var), freq1 = 0, freq2 = 1) # adding a line with the reference base fixed in sample 2 but not sample 1
      position <- rbind(position, extra)
    }
  }
  
  if(nrow(position) > 1)
  { # if after all this there is only 1 row in the position then the variant base is fixed in both samples and not interesting to us. We don't look at all the sites that have the reference base only in both.
    return(position)
  }else
  {
    return(position[F,])
  }
}

get_freqs <- function(pairs, snv)
{
  snv <- subset(snv, ALV_ID %in% pairs, select = c(ALV_ID, mutation, chr, pos, ref, var, freq.var, season, pcr_result))
  
  # just need the snv in the samples we're looking at.
  if(nrow(snv) > 0)
  { # There are mutations.
    # We only compare samples from the same season and strain so in each case the reference base is the same.
    stopifnot(length(unique(snv$season)) == 1, length(unique(snv$pcr_result)) == 1)
    
    mut_table <- tidyr::spread(snv, ALV_ID, freq.var, fill = 0) # a data frame with mutation down the first row and then frequency in either sample in the next 2.
    
    mut_table$ALV_ID1 <- pairs[1] # add column with first sample ID
    mut_table$ALV_ID2 <- pairs[2] # add column with second sample ID
    names(mut_table)[which(names(mut_table) == as.character(pairs[1]))] <-'freq1'
    names(mut_table)[which(names(mut_table) == as.character(pairs[2]))] <-'freq2'
    
    # This function can only be run on samples that qualified for variant identification.
    # If no variants were found in the sample then the SPECID will be missing from mut_table column and so
    # freq1 or freq2 will be missing since nothing was found we set that to 0 here and add the column.
    # equal compare will replace these cases with the reference at 1.
    if(!('freq1' %in% names(mut_table)))
    {
      mut_table$freq1 <- 0
    }
    if(!('freq2' %in% names(mut_table)))
    {
      mut_table$freq2 <- 0
    }
    mut_table <- dplyr::select(mut_table, mutation, chr, pos, ref, var, season, pcr_result, freq1,freq2, ALV_ID1, ALV_ID2)
    mut_table <- mut_table[order(mut_table$chr, mut_table$pos),]
    
    # Fill in differences based on inferred. In this case we are left with only sites that are polymorphic in one or between the 2 samples
    mut_table %>% dplyr::group_by(chr,pos) %>% dplyr::do(equal_compare(.)) -> all.freq
    
    if(nrow(all.freq) > 0)
    { # Some differences exist
      return(all.freq)
    }
    else
    { # some snv were initially present (before frequency filter) but after taking into account differences in inference the samples are the same.
      x <- tibble(mutation = NA, freq1 = NA, freq2 = NA, ALV_ID1 = NA, ALV_ID2 = NA, chr = NA, ref = NA, "pos" = NA, var = NA)
      return(x[F,])
    }
  }
  else
  { # No variants found in either sample
    x <- tibble(mutation = NA, freq1 = NA, freq2 = NA, ALV_ID1 = NA, ALV_ID2 = NA, chr = NA, ref = NA, "pos" = NA, var = NA)
    return(x[F,])
  }
}

dist_tp <- function(pairs, snv)
{
  data.df <- get_freqs(pairs = pairs, snv = snv)
  
  if(nrow(data.df) == 0) # This is the case when there are no variants found in either sample
  {
    d = 0
  } else
  {
    y <- as.matrix(cbind(data.df$freq1, data.df$freq2))
    d = dist(t(y), method = "manhattan")
  }

  return(as.numeric(d))
}

# ================== Get Transmission Pairs: JT's functions, modified for my metadata set =========================

finding_valid <- function(tp)
{
  finding_valid_helper <- function(x)
  {
    x <- dplyr::mutate(x, diff_onset = onset2 - onset1)
    valid_onsets <- x$diff_onset[x$diff_onset > 0] # valid pairs at this point must have a difference on onset date >0.
    
    # The recipient is sick on the house index case and onset1==onset2 and there is no other possible donor date - these are all sick on the first day and need to be handeled in both directions
    if(identical(x$onset1, x$onset2) & all(x$onset2 == min_onset))
    {
      x$valid <- TRUE
      x$double <- TRUE
    }
    else if(length(valid_onsets) > 0)
    {
      # There are some valid donors
      # If there is only 1 valid donor here with onset date closest to that of recipient (i.e donors not sick on the same day)
      if(length(which(x$diff_onset == min(valid_onsets))) == 1)
      {
        x$valid[x$diff_onset == min(valid_onsets)] <- TRUE
      }
    }
    
    x %>% select(-diff_onset) -> x
    return(x)
  } # end finding_valid_helper
  
  # Start main function: finding_valid
  min_onset <- min(tp$onset1)
  
  if(nrow(tp) == 1)
  { 
    tp$valid <- TRUE
    
    # There is only one pair and they are sick on the same day. Valid, but we will need to look at transmission both ways. Set double to TRUE.
    if(tp$onset1 == tp$onset2)
    {
      tp$double <- TRUE
    }
    return(tp)
  }
  else # There are multiple possible pairs here. For this, we will look at each recipeint and validate the pair that has the minimum diff in onset date but whose difference is >0.
  {
    tp %>% dplyr::group_by(ENROLLID2) %>% dplyr::do(finding_valid_helper(.)) -> tp_multi
    return(tp_multi)
  }
}


get_tp <- function(meta_one, interval)
{ 
  
  getting_tp_helper <- function(df, interval)
  {
    house <- unique(df$HOUSE_ID)
    stopifnot(length(house) == 1, length(unique(df$pcr_result)) == 1) # Verify only 1 house here and only 1 strain
    if(length(unique(df$ENROLLID)) != length(unique(df$SPECID))) stop("Please supply a data frame with only one entry/ENROLLID")
    pairs <- tibble(ENROLLID1 = NA, ENROLLID2 = NA, onset1 = NA, onset2 = NA, transmission = NA, valid = NA)[F,]
    
    if(length(unique(df$ENROLLID)) > 1) # Verify more than 1 person in the household
    {
      df <- df[order(df$onset, df$ENROLLID, decreasing = F),] # Order the data frame based on onset date. Earliest is first.
      
      for(i in 1:(nrow(df) - 1)) # Get all pairs
      {
        ENROLLID1 = df$ENROLLID[i]
        onset1 = df$onset[i]
        for(j in (i+1):nrow(df))
        {
          ENROLLID2 = df$ENROLLID[j]
          onset2 = df$onset[j]
          transmission = onset2 - 1
          if((onset2 - onset1) <= interval)
          {
            out <- tibble(ENROLLID1 = ENROLLID1, ENROLLID2 = ENROLLID2, onset1 = onset1, onset2 = onset2, transmission = transmission, valid = F) # all start as not valid
            pairs <- rbind(pairs, out) # We don't want that first line with the NA this just adds the out data frame with the transmission pair if it meets the onset date
          }
        }
      }
      
      pairs <- filter(pairs, !(is.na(onset1)) & !(is.na(onset2)))
      if(nrow(pairs) > 0)
      {
        valid_pairs <- finding_valid(pairs)
        return(valid_pairs)
      } else
      {
        return(pairs[F,])
      }
    } else 
    {
      return(pairs[F,])
    }
  }
  
  # Main function: get_tp
  meta_one %>% dplyr::group_by(HOUSE_ID, pcr_result, season) %>% dplyr::do(getting_tp_helper(., interval)) -> result
  return(result)
}

# ====================== Function to select a single sample per individual ==================

# Method used in IAV paper to select one sample per individual to run the L1-norm comparison. It yields the same results as the method I use below (i.e. take home sample if available).

only_one <- function(meta_df)
{
  
  only_one_helper<-function(x)
  {
    # if there is only one sample then we are done
    if(nrow(x)==1)
    {
      return(x)
    } else if(nrow(x) == 2)
    {
      # There are two samples here with this person and pcr result
      # pick the ones that we sequenced
      pick <- which(x$sequenced == TRUE)
      if(length(pick) == 1)
      {
        #if there is only one then we're done
        return(x[pick, ])
      }
      else if(length(pick) == 0)
      {
        # we didn't sequence any of these
        pick <- 1
        return(x[pick, ]) # take the first one
      }
      else if(length(pick) == 2)
      {
        # we sequenced both
        # See how many qualified for snv calling
        pick <- which(x$coverage_qualified == TRUE)
        # Either 0 or 2 did
        if(length(pick) == 0 | length(pick) == 2)
        {
          # take the one with a higher titer.
          pick <- which(x$genome_copy_per_ul == max(x$genome_copy_per_ul))
        }
        return(x[pick, ])
      }
    } else
    {
      stop("More than 2 samples for this person? That shouldn't happen!")
    }
  }
  
  out <- meta_df %>% dplyr::group_by(ENROLLID,pcr_result,season) %>%
    dplyr::do(only_one_helper(.))
  return(out)
}

meta_long <- read_csv("data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv")
meta_long <- mutate(meta_long, sequenced = TRUE)
meta_long_one <- only_one(meta_long)

meta_long_one %>%
  group_by(season, pcr_result) %>%
  do({expand.grid(.$ALV_ID, .$ALV_ID)}) %>%
  ungroup() -> all_possible_pairs

# ================== Get specimens to compare =========================

# For each individual, we have all infections that have at least one quality sample. That's in meta_wide_seqd. All of these have SNV-quality data.
# Take home sample, unless not available, then just take clinic sample.

meta_wide_seqd <- read_csv("data/metadata/metadata_wide_quality_sequence.csv")

meta_wide_seqd <- mutate(meta_wide_seqd, sampleForDistanceAnalysis = ifelse(haveQualityHomeSample == TRUE, home_spec_accn, SPECID))
meta_wide_unique <- filter(meta_wide_seqd, !(ENROLLID == "50003" & cluster == 0)) # 50003 is in here twice. Infected twice in one season Remove the one that's not in a cluster (no household transmission).

#meta_wide_unique <- mutate(meta_wide_unique, collect_forAnalysis = ifelse(haveQualityHomeSample == TRUE, home_spec_date, collect)) # For some reason, ifelse was changing the format of the date. It resulted in the same thing, but I changed the code here to be more clear.
meta_wide_unique_home <- filter(meta_wide_unique, haveQualityHomeSample == TRUE)
meta_wide_unique_clinic <- filter(meta_wide_unique, haveQualityHomeSample == FALSE)
meta_wide_unique_home <- mutate(meta_wide_unique_home, collect_forAnalysis = home_spec_date)
meta_wide_unique_clinic <- mutate(meta_wide_unique_clinic, collect_forAnalysis = collect)
meta_wide_unique <- rbind(meta_wide_unique_home, meta_wide_unique_clinic)

meta_wide_unique <- mutate(meta_wide_unique, pcr_result = ifelse(str_detect(pcr_result, "Y") | pcr_result == "B", "B/YAM", "B/VIC")) # Format lineage

meta_wide_unique %>%
  group_by(season, pcr_result) %>%
  do({expand.grid(.$sampleForDistanceAnalysis, .$sampleForDistanceAnalysis)}) %>%
  ungroup() -> all_possible_pairs

names(all_possible_pairs) <- c("season", "pcr_result", "ALV_ID_1", "ALV_ID_2") # rename columns

all_possible_pairs <- mutate(all_possible_pairs,
                           ENROLLID1 = meta_wide_unique$ENROLLID[match(x = ALV_ID_1, meta_wide_unique$sampleForDistanceAnalysis)],
                           ENROLLID2 = meta_wide_unique$ENROLLID[match(x = ALV_ID_2, meta_wide_unique$sampleForDistanceAnalysis)],
                           time_onset = as.numeric(meta_wide_unique$onset[match(x = ALV_ID_2, meta_wide_unique$sampleForDistanceAnalysis)]) - 
                             as.numeric(meta_wide_unique$onset[match(ALV_ID_1, meta_wide_unique$sampleForDistanceAnalysis)]),
                           time_collect = as.numeric(meta_wide_unique$collect_forAnalysis[match(x = ALV_ID_2, meta_wide_unique$sampleForDistanceAnalysis)]) - 
                             as.numeric(meta_wide_unique$collect_forAnalysis[match(ALV_ID_1, meta_wide_unique$sampleForDistanceAnalysis)]))

all_possible_pairs <- subset(all_possible_pairs, ALV_ID_1 != ALV_ID_2) # remove self-comparisons
all_possible_pairs <- subset(all_possible_pairs, time_onset >= 0) # only those with ALV_ID_2 sick later than or equal to ALV_ID_1

#write.csv(all_possible_pairs, file = "data/processed/all_possible_pairs.csv")

# ================== Distance calculation =========================

### Remove the mixed infection by ENROLLID, if present.
mixed_infection <- c("50425")
all_possible_pairs_nomixed <- filter(all_possible_pairs, !(ENROLLID1 %in% mixed_infection) & !(ENROLLID2 %in% mixed_infection))

possible_pairs.dist <- plyr::adply(all_possible_pairs_nomixed, 1, summarize, L1_norm = dist_tp(c(ALV_ID_1, ALV_ID_2), snv = qual)) # Get the genetic distance between pairs! Takes the longest.

ggplot(possible_pairs.dist, aes(L1_norm)) + geom_histogram(binwidth = 1) # a quick view.

# ================== Sift by epidemiological linkage =========================

# Add house info for each pair.
possible_pairs.dist <- mutate(possible_pairs.dist,
                            HOUSE_ID1 = meta_wide_unique$HOUSE_ID[match(x = ALV_ID_1, meta_wide_unique$sampleForDistanceAnalysis)],
                            HOUSE_ID2 = meta_wide_unique$HOUSE_ID[match(x = ALV_ID_2, meta_wide_unique$sampleForDistanceAnalysis)],
                            Household = HOUSE_ID1==HOUSE_ID2) # Are these household pairs? True or False

# Get putative transmission pairs
linkage_interval <- 7
all_pairs <- get_tp(meta_wide_unique, interval = linkage_interval) # Within 7 days of each other
all_pairs$pair_id <- 1:nrow(all_pairs)
all_pairs <- ungroup(all_pairs)
valid_pairs <- filter(all_pairs, valid == TRUE) # Here we get 16; one includes the mixed infection
valid_pairs <- mutate(valid_pairs, same_onset_day = ifelse(onset1 == onset2, TRUE, FALSE))
possible_pairs.dist <- left_join(possible_pairs.dist, select(valid_pairs, season, pcr_result, ENROLLID1, ENROLLID2, valid))
possible_pairs.dist$valid[is.na(possible_pairs.dist$valid)] <- FALSE


# ================== Generate distance cutoffs =========================

cutoffs <- tibble(L1_norm = quantile(possible_pairs.dist$L1_norm[possible_pairs.dist$Household == F], probs = seq(0, 1, 0.05)))
cutoffs$threshold <- seq(0, 1, 0.05)
cutoffs %>% 
  rowwise() %>% 
  mutate(valid_pairs = nrow(possible_pairs.dist[(possible_pairs.dist$valid == TRUE & possible_pairs.dist$L1_norm < L1_norm),])) -> cutoffs

cutoff <- cutoffs$L1_norm[cutoffs$threshold == 0.05] # below 5th percentile

possible_pairs.dist <- mutate(possible_pairs.dist, quality_distance = L1_norm <= cutoff)

valid_pairs <- left_join(valid_pairs, select(possible_pairs.dist, ENROLLID1, ENROLLID2, quality_distance, L1_norm))
valid_pairs <- mutate(valid_pairs, ALV_ID_1 = meta_wide_unique$sampleForDistanceAnalysis[match(x = ENROLLID1, meta_wide_unique$ENROLLID)])
valid_pairs <- mutate(valid_pairs, ALV_ID_2 = meta_wide_unique$sampleForDistanceAnalysis[match(x = ENROLLID2, meta_wide_unique$ENROLLID)])

# ============================== Plot histogram of L1-norm ===================================

palette = wesanderson::wes_palette("FantasticFox1")
possible_pairs.dist %>% filter(valid == TRUE | Household == FALSE) %>% select(season, ALV_ID_1, ALV_ID_2, L1_norm, valid, Household) -> L1norm_plot_data

L1norm_plot <- ggplot(L1norm_plot_data, aes(x = L1_norm, fill = as.factor((valid-1)*-1), y = ..ncount..)) +
  geom_histogram(binwidth = 7.5, boundary = 0, position = 'dodge') +
  scale_fill_manual(name = "", labels = c("Household pair", "Community pair"), values = palette[c(3,4)]) +
  xlab("L1 Norm") + ylab("Normalized Count") +
  theme(legend.position = c(0.5, 0.5)) +
  geom_segment(aes(x = cutoff, xend = cutoff, y = 0, yend = 1), linetype = 2, color = palette[5], size = 0.4) + theme_classic() +
  theme(text = element_text(size = 28), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20)) # save as PDF, 15 by 10

#ggsave(plot = L1norm_plot, filename = "results/plots/L1norm_square.pdf", device = "pdf", width = 15, height = 10)

# ============================== Number of SNVs per transmission pair sample ===================================

read_csv("data/processed/snv_qual_meta.o.csv") %>% select(ALV_ID, iSNV) %>% rename(sampleForDistanceAnalysis = ALV_ID) -> snv_by_sample
merged <- merge(snv_by_sample, select(meta_wide_unique, ENROLLID, sampleForDistanceAnalysis), by = "sampleForDistanceAnalysis")
select(merged, iSNV, ENROLLID) %>% rename(ENROLLID1 = ENROLLID) -> merged1
select(merged, iSNV, ENROLLID) %>% rename(ENROLLID2 = ENROLLID) -> merged2
valid_pairs_snv <- left_join(valid_pairs, merged1)
valid_pairs_snv <- rename(valid_pairs_snv, iSNV_in_donor = iSNV)
valid_pairs_snv <- left_join(valid_pairs_snv, merged2)
valid_pairs_snv <- rename(valid_pairs_snv, iSNV_in_recipient = iSNV)

# ============================== Save the output ===============================

write.csv(possible_pairs.dist, file = "data/processed/possible.pairs.dist.csv")
write.csv(valid_pairs_snv, file = "data/processed/transmission_pairs.csv")
write.csv(meta_wide_unique, file = "data/processed/meta_wide_unique_forTransmissionAnalysis.csv")

