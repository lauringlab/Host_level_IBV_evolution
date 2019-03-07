
### Project: IBV intra-host
### Purpose: Explore SNVs across true transmission pairs. 

# ================== Read in data, import packages =========================

library(tidyverse)
library(lubridate)

setwd("/Users/avalesano/Documents/MSTP/LauringLab/Host_level_IBV_evolution/scripts/")

meta_long <- read_csv("../data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv")
meta_wide <- read_csv("../data/metadata/metadata_wide_quality_sequence.csv")
transmission_pairs <- read_csv("../data/processed/transmission_pairs.csv")
var_qual <- read_csv("../data/processed/qual.snv.csv")
var_no_cut <- read_csv("../data/processed/no_freq_cut.qual.snv.csv")

# ================== Get Arranged Data: JT's functions, modified for my metadata set =========================

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
    # Do it again if there are no varaints in the second sample.
  } else if(pos_sum_2 == 0)
  { # there isn't a variant in the second sample. The reference is fixed. also I told you so. (line 196)
    x <- which(position$ref == position$var)
    stopifnot(length(x) < 2) # should be 1 or 0
    if(length(x) == 1)
    { # The reference has been infered in the first sample but there are no varaints in the second so it was not infered
      position$freq2[x] <- 1
    } else if(length(x) == 0)
    { # the reference was not infered in the first sample and no variants here. Add the reference now
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
  # Input: List of SNV calls and a dataframe of one ID pair.
  # Output: Dataframe comparing each SNV's frequency in both samples. Only those where the frequency changes from one sample to another.
  
  snv <- subset(snv, ALV_ID %in% pairs, select = c(ALV_ID, mutation, chr, pos, ref, var, freq.var, season, pcr_result))
  
  if(nrow(snv) > 0)
  { # There are mutations.
    # We only compare samples from the same season and strain so in each case the reference base is the same
    stopifnot(length(unique(snv$season)) == 1, length(unique(snv$pcr_result)) == 1)

    mut_table <- tidyr::spread(snv, ALV_ID, freq.var, fill = 0) # a data frame with mutation down the first row and then frequency in either sample in the next 2.
    
    mut_table$ALV_ID_1 <- pairs[1] # add column with first sample ID
    mut_table$ALV_ID_2 <- pairs[2] # add column with second sample ID
    names(mut_table)[which(names(mut_table) == as.character(pairs[1]))] <- 'freq1' # rename this column as the frequency in the first sample
    names(mut_table)[which(names(mut_table) == as.character(pairs[2]))] <- 'freq2'
    
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
    mut_table <- dplyr::select(mut_table, mutation, chr, pos, ref, var, season, pcr_result, freq1, freq2, ALV_ID_1, ALV_ID_2)
    mut_table <- mut_table[order(mut_table$chr, mut_table$pos),]

    # Fill in differences based on inferred. In this case we are left with only sites that are polymorphic in one or between the 2 samples
    mut_table %>% dplyr::group_by(chr, pos) %>% dplyr::do(equal_compare(.)) -> all.freq

    if(nrow(all.freq) > 0)
    { # Some differences exist
      return(all.freq)
    }
    else
    { # some snv were initially present (before frequency filter) but after taking into account differences in inference the samples are the same.
      x <- tibble(mutation = NA, freq1 = NA, freq2 = NA, ALV_ID_1 = NA, ALV_ID_2 = NA, chr = NA, ref = NA, "pos" = NA, var = NA)
      return(x[F,])
    }
  }
  else
  { # No variants found in either sample
    x <- tibble(mutation = NA, freq1 = NA, freq2 = NA, ALV_ID_1 = NA, ALV_ID_2 = NA, chr = NA, ref = NA, "pos" = NA, var = NA)
    return(x[F,])
  }
}

# Yields mutations that are polymorphic in the donor.
polish_freq <- function(df, relative, min_freq = 0)
{
  rel_quo <- enquo(relative)
  min_freq_quo <- enquo(min_freq)
  m = paste0("Requiring polymophic at least in ", quo_name(rel_quo))
  message(m)
  
  # We apply a freqeuncy cut off. it will be at least 0
  frequency_cut <- rlang::expr(UQ(rel_quo) > UQ(min_freq_quo))
  df <- dplyr::filter(df, UQ(frequency_cut))
  
  # We want polymorphic sites in the given sample. Due to frequency oddities the major allels may not
  # be exactly 1-min_freq i.e. we don't set any frequencies by hand. minor frequencies are set by deepsnv
  # and this includes an estimated error rate. Major frequencies are set by allele counts.
  
  frequency_cut_minor <- rlang::expr(UQ(rel_quo) < 1)
  
  # How many alleles have frequencies in the choosen sample above the min and below 1 at each loci
  df %>% dplyr::group_by(chr, pos, ALV_ID_1, ALV_ID_2) %>% dplyr::summarize(alleles = length(which(UQ(frequency_cut) & UQ(frequency_cut_minor)))) -> counts
  
  #View(counts)
  
  #if(max(counts$alleles > 2)) stop(paste0("More that two alleles found at a position in a sample"))
  df <- merge(df, counts)
  
  return(subset(df, alleles == 2, select = -c(alleles)))
}

# ================== JT's function: get specimens closest to time of transmission ====================

#' Finding SPECID for transmission pairs
#'
#' Identifies the SPECID collected nearest a given date, and on the correct side of
#' the date when specified. The criteria is
#'
#' 1) Sample closest to transmission that are on the "right side" of tranmission
#' (before transmission for donor if possible)
#'
#' 2) Titer is the tie breaker when applicable.
#'
#' 3) If no iSNV were found in a donor sample and iSNV==T we take the sample with iSNV present.
#'
#' @param meta data frame with meta data con
#' @param date The estimated trasmission date
#' @param enrollid The enrollid
#' @param case c("donor","recipient")
#' @param iSNV boolan if TRUE then we take the sample that has iSNV if the
#' "best fit" doesn't. This is the historic use of the function. Recall that we
#' do not have samples before transmission in many cases anyway.
#'
#' @return The SPECID that best fits the critea
#' @examples
#'
#' d<-small_meta$collect[4]+1
#' get_close(small_meta,d,enrollid = "50001",case = "Donor",iSNV = TRUE)
#' @export


get_close <- function(meta, date, enrollid, case, iSNV = TRUE)
{
  id_snv = filter(meta, ENROLLID == enrollid & coverage_qualified == TRUE)
  
  if(nrow(id_snv) == 1) # if there is only one sample from this individual
  { 
    return(id_snv$ALV_ID)
  }
  else if(nrow(id_snv) == 2) # there are two samples for this individual
  { 
    id_snv = dplyr::mutate(id_snv, dtrans = collect_date_master - date)
    
    # For donors we want the largest dtrans possible less than 0 or the smallest
    # possitive dtrans. For the recipient we want the smallest possitive dtrans.
    # It should be impossible for the recipient to have a negative dtrans.
    
    if(all(id_snv$dtrans < 0))
    {
      # all are less than 0 so this puts the sample with least negative dtrans on top. 
      # titer is the tie braker
      id_snv = id_snv[order(id_snv$dtrans, id_snv$genome_copy_per_ul, decreasing = TRUE),]
    } else
    {
      # There is 1 or no dtrans less than 0 so this puts the sample with the
      # lowest dtrans (either the negative 1 or the smallest positive 1) on top.
      # titer is the tie braker
      id_snv = id_snv[order(id_snv$dtrans, -id_snv$genome_copy_per_ul, decreasing = FALSE),]
    }
    if(id_snv$iSNV[1] == 0 & id_snv$iSNV[2] > 0 & iSNV == TRUE)
    {
      # finally if there were no polymorphisms in the first sample and
      # there were a few in the second then we take the second
      warning(paste0(case," ALV_ID ", id_snv$ALV_ID[2], "  being used for ENROLLID ",
                     id_snv$ENROLLID[1], " because no iSNV were found in the prefered
                     ALV_ID" , id_snv$ALV_ID[1]),"\n")
      return(id_snv$ALV_ID[2])
    } else
    {
      if(id_snv$genome_copy_per_ul[1] < id_snv$genome_copy_per_ul[2])
      {
        message(paste0("ALV_ID ", id_snv$ALV_ID[1],
                       "  being used for ENROLLID ", id_snv$ENROLLID[1],
                       " based on the time to the transmission date",
                       "even though this sample has a lower titer than " , id_snv$ALV_ID[2],"\n"))
      }
      return(id_snv$ALV_ID[1])
    }
  }else
  {
    stop("Stopping because neither 1 nor 2 specid found for this sample! Yikes!")
  }
}


# ================== Arrange data ===========================

useful_transmission_pairs <- filter(transmission_pairs, valid == TRUE & quality_distance == TRUE)

# Account for dual directionality of pairs with same onset day
flipped <- filter(useful_transmission_pairs, double == TRUE)
flipped <- plyr::rename(flipped, c("ENROLLID1" = "ENROLLID2", "ENROLLID2" = "ENROLLID1", "onset1" = "onset2", "onset2" = "onset1"))
flipped$pair_id = flipped$pair_id + 0.5
useful_transmission_pairs_withDual <- rbind(useful_transmission_pairs, flipped)

### Want to get the specimens that are closest to the time of transmission (not what was done for L1-norm).
# First prep the metadata.
meta_long <- mutate(meta_long, collect_date_master = mdy(ifelse(is_home_spec, home_spec_date, collect)))
var_qual_minority <- filter(var_qual, freq.var < 0.98)
minority.count <- plyr::ddply(var_qual_minority, ~ALV_ID, summarize, iSNV = length(unique(mutation)))
meta_long <- left_join(meta_long, minority.count)
meta_long$iSNV[is.na(meta_long$iSNV)] <- 0 # No variants in these samples

useful_transmission_pairs_withDual %>% rowwise() %>% mutate(ALV_ID_1 = get_close(meta_long, transmission, ENROLLID1, "donor"), ALV_ID_2 = get_close(meta_long, transmission, ENROLLID2, "recipient")) -> useful_transmission_pairs_withDual
useful_transmission_pairs_withDual <- mutate(useful_transmission_pairs_withDual, collect1 = meta_long$collect_date_master[match(ALV_ID_1, meta_long$ALV_ID)], collect2 = meta_long$collect_date_master[match(ALV_ID_2, meta_long$ALV_ID)])

# Remove the mixed infection ENROLLID, if they are present.
mixed <- c("50425")
useful_transmission_pairs_withDual <- filter(useful_transmission_pairs_withDual, !(ENROLLID1 %in% mixed) & !(ENROLLID2 %in% mixed))

trans_freq <- plyr::adply(useful_transmission_pairs_withDual, 1, function(x) {get_freqs(c(x$ALV_ID_1, x$ALV_ID_2), var_qual)})
write.csv(trans_freq, file = "../data/processed/trans_freq.csv")

# Reduce to sites that are polymorphic in the donor.
trans_freq.comp <- polish_freq(trans_freq, freq1, 0.02)
trans_freq.comp$found <- trans_freq.comp$freq2 > 0.02 # Was it found in the second sample?
write.csv(trans_freq.comp, file = "../data/processed/transmission_pairs_freq.poly_donor.csv")

no_cut_trans_freq <- plyr::adply(useful_transmission_pairs_withDual, 1, function(x){get_freqs(c(x$ALV_ID_1, x$ALV_ID_2), var_no_cut)})
write.csv(no_cut_trans_freq, file = "../data/processed/no_cut_trans_freq.csv")

# Reduce to sites that are polymorphic in the donor.
no_cut_trans_freq.comp <- polish_freq(no_cut_trans_freq, freq1, 0)
no_cut_trans_freq.comp$found <- no_cut_trans_freq.comp$freq2 > 0 # Was it found in the second sample?
write.csv(no_cut_trans_freq.comp, file =  "../data/processed/no_cut_transmission_pairs_freq.poly_donor.csv")


# ================== Plot *donor* SNVs by frequency in recipient vs. donor ===========================

trans_freq.p <- ggplot(trans_freq, aes(x = freq1, y = freq2)) + geom_point() + xlab("Frequency in donor") + ylab("Frequency in recipient") + theme_classic()
trans_freq.comp.p <- ggplot(trans_freq.comp, aes(x = freq1, y = freq2)) + geom_point() + xlab("Frequency in donor") + ylab("Frequency in recipient") + theme_classic()

no_cut_trans_freq.p <- ggplot(no_cut_trans_freq, aes(x = freq1, y = freq2)) + geom_point() + xlab("Frequency in donor") + ylab("Frequency in recipient") + theme_classic()
no_cut_trans_freq.comp.p <- ggplot(no_cut_trans_freq.comp, aes(x = freq1, y = freq2)) + geom_point() + xlab("Frequency in donor") + ylab("Frequency in recipient") + theme_classic()

ggsave(plot = trans_freq.p, filename = "../results/plots/RecipientVsDonorAll.jpg", device = "jpeg")
ggsave(plot = trans_freq.comp.p, filename = "../results/plots/RecipientVsDonor_DonorPolymorphic.jpg", device = "jpeg")
ggsave(plot = no_cut_trans_freq.p, filename = "../results/plots/RecipientVsDonorAll_NoCutoff.jpg", device = "jpeg")
ggsave(plot = no_cut_trans_freq.comp.p, filename = "../results/plots/RecipientVsDonor_DonorPolymorphic_NoCutoff.jpg", device = "jpeg")
