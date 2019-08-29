

### Project: IBV intra-host
### Purpose: Process SNV for downstream analysis.

# ============================ Packages, read in data =======================================

library(tidyverse)

variants <- read_csv("data/processed/all.variants.csv")
meta <- read_csv("data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv")


# =========================== Join variant data and metadata ===============================

metadata_by_seq <- read_csv("data/metadata/Metadata_By_Seq_withALVID.csv")
variants <- rename(variants, SampleNumber = Id)
variants_withALVID <- merge(variants, select(metadata_by_seq, SampleNumber, ALV_ID), by="SampleNumber")
variants_with_meta <- merge(variants_withALVID, meta, by="ALV_ID")

# =========================== Filter by coverage ===============================

poor_coverage_variants <- filter(variants_with_meta, coverage_qualified == FALSE)
good_coverage_variants <- filter(variants_with_meta, coverage_qualified == TRUE)

#freq.good.hist <- ggplot(good_coverage_variants, aes(freq.var)) + geom_histogram(binwidth = 0.01) + ylim(0, 2000)
#freq.poor.hist <- ggplot(poor_coverage_variants, aes(freq.var)) + geom_histogram(binwidth = 0.01) + ylim(0, 2000)
#freq.raw.hist <- ggplot(variants_with_meta, aes(freq.var)) + geom_histogram(binwidth = 0.01) + ylim(0, 2000)

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

# After filtering false positives, there could be a monomorphic site that has frequency < 1. This sets those sites to 1.
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

diverse_variants %>% group_by(ALV_ID) %>% do(collapse_localcov(.)) -> quality_var

quality_var$class_factor = NA
quality_var$class_factor[grep("Noncoding", quality_var$Class)] <- "Noncoding"
quality_var$class_factor[grep('Syn', quality_var$Class)] <- "Synonymous"
quality_var$class_factor[grep('Nonsyn', quality_var$Class)] <- "Nonsynonymous"
quality_var$class_factor <- as.factor(quality_var$class_factor)
quality_var <- subset(quality_var, class_factor != "Noncoding")

quality_var <- diverse_sites(quality_var, 1, season, pcr_result, pos, chr)

no_freq_cut <- quality_var
quality_var %>% filter(freq.var > 0.02) -> quality_var

quality_var_monomorphic <- monomorphic(quality_var, ALV_ID, season, chr, pos, pcr_result)
no_freq_cut_monomorphic <- monomorphic(no_freq_cut, ALV_ID, season, chr, pos, pcr_result)

quality_var_monomorphic <- subset(quality_var_monomorphic, !(ref == var & freq.var == 1))
no_freq_cut_monomorphic <- subset(no_freq_cut_monomorphic, !(ref == var & freq.var == 1)) 

write.csv(x = no_freq_cut_monomorphic, file = "data/processed/no_freq_cut.qual.snv.csv")
write.csv(x = quality_var_monomorphic, file = "data/processed/qual.snv.csv")

# ======================== Same thing, but no resolution by sample duplicates =======================

diverse_variants -> quality_var # Only difference is that there is no collapse step
quality_var$class_factor = NA
quality_var$class_factor[grep("Noncoding", quality_var$Class)] <- "Noncoding"
quality_var$class_factor[grep('Syn', quality_var$Class)] <- "Synonymous"
quality_var$class_factor[grep('Nonsyn', quality_var$Class)] <- "Nonsynonymous"
quality_var$class_factor <- as.factor(quality_var$class_factor)
quality_var <- subset(quality_var, class_factor != "Noncoding")
quality_var <- diverse_sites(quality_var, 1, season, pcr_result, pos, chr)
quality_var %>% filter(freq.var > 0.02) -> quality_var
quality_var <- monomorphic(quality_var, ALV_ID, season, chr, pos, pcr_result)
quality_var <- subset(quality_var, !(ref == var & freq.var == 1))
write.csv(quality_var, file = "data/processed/qual.not.collapsed.snv.csv")

