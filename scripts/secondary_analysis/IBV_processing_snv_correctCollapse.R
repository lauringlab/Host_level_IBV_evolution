

### Project: IBV intra-host
### Purpose: Process SNV for downstream analysis.

# ============================ Packages, read in data =======================================

library(magrittr)
library(tidyverse)
library(ggplot2)

setwd("/Users/avalesano/Documents/MSTP/LauringLab/Host_level_IBV_evolution/scripts/")

variants <- read_csv("../data/processed/all.variants.csv")
meta <- read_csv("../data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv")

# =========================== Join variant data and metadata ===============================

metadata_by_seq <- read_csv("../data/metadata/Metadata_By_Seq_withALVID.csv")
variants <- rename(variants, SampleNumber = Id)
variants_withALVID <- merge(variants, select(metadata_by_seq, SampleNumber, ALV_ID), by="SampleNumber")
variants_with_meta <- merge(variants_withALVID, meta, by="ALV_ID")

# =========================== Filter by coverage ===============================

poor_coverage_variants <- filter(variants_with_meta, coverage_qualified == FALSE)
good_coverage_variants <- filter(variants_with_meta, coverage_qualified == TRUE)


# =========================== Function for getting iSNV found in both sequencing replicates ===============================

sift_dups <- function(df)
{
  if(nrow(df) > 2) stop("Too many mutations here")
  
  df <- dplyr::mutate(df, coverage = cov.tst.bw + cov.tst.fw)
  higher_qual <- subset(df, coverage == max(df$coverage))
  if(nrow(higher_qual) > 1)
  { 
    higher_qual <- higher_qual[1,]
  }
  return(higher_qual)
}

# Based on JT's function "quality". Uses local coverage, not average.
collapse_localcov <- function(df)
{
  stopifnot(length(unique(df$ALV_ID))==1)
  
  replicates <- unique(df$seq_replicates)
  
  if(replicates == 2)
  {
    df %>% dplyr::group_by(mutation) %>% dplyr::summarize(found = length(mutation)) -> count_mutations
    
    stopifnot(max(count_mutations$found) < 3)
    
    count_mutations %>% dplyr::filter(found == 2) -> good_mut
      
    df %>% dplyr::filter(mutation %in% good_mut$mutation) -> good_var
      
    good_var %>% dplyr::group_by(mutation) %>% dplyr::do(sift_dups(.)) -> dups_good
    
    return(dplyr::ungroup(dups_good))
    
  } else if(replicates == 1)
  {
    return(df)
  }
}


# =========================== Some diagnostic plots ===============================

freq.good.hist <- ggplot(good_coverage_variants, aes(freq.var)) + geom_histogram(binwidth = 0.01) + ylim(0, 2000) # lots of variants above 0.5 because that is the stringent freq cutoff in options file.
freq.poor.hist <- ggplot(poor_coverage_variants, aes(freq.var)) + geom_histogram(binwidth = 0.01) + ylim(0, 2000) # some minority variants filtered here (832, exactly).
freq.raw.hist <- ggplot(variants_with_meta, aes(freq.var)) + geom_histogram(binwidth = 0.01) + ylim(0, 2000)

# =========================== Filtering: helper functions ===============================

# Get sites that have more than one allele present in a given grouping (chr, season, lineage, etc.)
diverse_sites <- function(df, cutoff,...)
{
  group_var <- rlang::quos(...)
  df %>% dplyr::group_by(!!!group_var) %>% dplyr::summarize(alleles = length(unique(var))) -> df_alleles
  
  df <- dplyr::left_join(df, df_alleles)
  stopifnot(any(!is.na(df$alleles)))
  df %>% dplyr::filter(alleles > cutoff) %>% dplyr::select(-alleles) -> df
  return(df)
}

# After filtering FPs, there could be a monomorphic site that has frequency <1. This sets those sites to 1.
monomorphic <- function(df,...)
{
  group_var <- rlang::quos(...)
  df %>% dplyr::group_by(!!!group_var) %>% dplyr::summarize(alleles = length(unique(var))) -> df_alleles
  
  df <- dplyr::left_join(df, df_alleles)
  stopifnot(any(!is.na(df$alleles)))
  df$freq.var[df$alleles == 1] <- 1
  df %>% dplyr::select(-alleles) -> df
  return(df)
}

# =========================== Filtering ===============================

diverse_variants <- diverse_sites(good_coverage_variants, 1, season, pcr_result, pos, chr)

freq.diverse.hist <- ggplot(diverse_variants, aes(freq.var)) + geom_histogram(binwidth = 0.01) + ylim(0, 4000)

diverse_variants %>% group_by(ALV_ID) %>% do(collapse_localcov(.)) -> qual

qual$class_factor = NA
qual$class_factor[grep("Noncoding", qual$Class)] <- "Noncoding"
qual$class_factor[grep('Syn', qual$Class)] <- "Synonymous"
qual$class_factor[grep('Nonsyn', qual$Class)] <- "Nonsynonymous"
qual$class_factor <- as.factor(qual$class_factor)
qual <- subset(qual, class_factor != "Noncoding")

qual <- diverse_sites(qual, 1, season, pcr_result, pos, chr)

no_freq_cut <- qual
qual %>% filter(freq.var > 0.02) -> qual

qualm <- monomorphic(qual, ALV_ID, season, chr, pos, pcr_result)
no_freq_cut <- monomorphic(no_freq_cut, ALV_ID, season, chr, pos, pcr_result)

qualms <- subset(qualm, !(ref == var & freq.var == 1))
no_freq_cut <- subset(no_freq_cut, !(ref == var & freq.var == 1)) 

write.csv(x = no_freq_cut, file = "../data/processed/no_freq_cut.qual.snv.csv")
write.csv(x = qualms, file = "../data/processed/qual.snv.csv")

# ======================== Same thing, but no resolution by sample duplicates =======================

diverse_variants -> qual
qual$class_factor = NA
qual$class_factor[grep("Noncoding", qual$Class)] <- "Noncoding"
qual$class_factor[grep('Syn', qual$Class)] <- "Synonymous"
qual$class_factor[grep('Nonsyn', qual$Class)] <- "Nonsynonymous"
qual$class_factor <- as.factor(qual$class_factor)
qual <- subset(qual, class_factor != "Noncoding")
qual <- diverse_sites(qual, 1, season, pcr_result, pos, chr)
qual %>% filter(freq.var > 0.02) -> qual
qualm <- monomorphic(qual, ALV_ID, season, chr, pos, pcr_result)
qualms <- subset(qualm, !(ref == var & freq.var == 1))
write.csv(x = qualms, file = "../data/processed/qual.not.collapsed.snv.csv")

