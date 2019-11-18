
### Project: IBV
### Purpose: Compare iSNV and nucleotide diversity across IAV/IBV with different frequency cutoffs.
### Working directory: Host_level_IBV_evolution

# ======================== Import packages, load data =========================

library(tidyverse)

meta <- read.csv("data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv")
meta_snv <- filter(meta, coverage_qualified == TRUE)
meta_IAV <- read.csv("data/metadata/IAV_all_meta.sequence_success.csv")
meta_IAV_snv <- filter(meta_IAV, snv_qualified == TRUE) %>% mutate(ALV_ID = SPECID)

no_freq_cut_var <- read.csv("data/processed/no_freq_cut.qual.snv.csv")
no_freq_cut_var <- filter(no_freq_cut_var, ENROLLID != "50425")
no_freq_cut_var <- mutate(no_freq_cut_var, site_coverage = cov.tst.fw + cov.tst.bw)
no_freq_cut_var <- mutate(no_freq_cut_var, allele_coverage = n.tst.fw + n.tst.bw)

no_freq_cut_var_IAV <- read.csv("data/processed/no_freq_cut_IAV.qual.snv.csv")
no_freq_cut_var_IAV <- filter(no_freq_cut_var_IAV, !SPECID %in% c("HS1530", "M54062" ,"MH8125", "MH8137" ,"MH8156" ,"MH8390"))
no_freq_cut_var_IAV <- mutate(no_freq_cut_var_IAV, site_coverage = cov.tst.fw + cov.tst.bw)
no_freq_cut_var_IAV <- mutate(no_freq_cut_var_IAV, allele_coverage = n.tst.fw + n.tst.bw)
no_freq_cut_var_IAV <- mutate(no_freq_cut_var_IAV, ALV_ID = SPECID)

plot.median <- function(x) {
  m <- median(x)
  c(y = m, ymin = m, ymax = m)
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

# =================== Functions to calculate nucleotide diversity ================

get_pairwise_distance <- function(df)
{
  if(nrow(df) > 2)
  {
    print(df)
    stopifnot(nrow(df) == 2) # should only have two rows here
  }
  
  minor <- filter(df, freq.var < 0.5)
  major <- filter(df, freq.var > 0.5)
  coverage <- major$site_coverage
  cov_factor <- (coverage * (coverage - 1))
  
  cov_minor <- minor$allele_coverage
  cov_major <- major$allele_coverage
  
  minor_factor <- cov_minor * (cov_minor - 1)
  major_factor <- cov_major * (cov_major - 1)
  
  D <- (cov_factor - (major_factor + minor_factor)) / cov_factor
  
  chr_pos <- paste0(as.character(major$chr), "_", as.character(major$pos))
  
  return(data.frame(site = chr_pos, D = D))
}

get_nt_diversity <- function(df)
{
  stopifnot(length(unique(df$ALV_ID)) == 1) # should be doing one sample at a time
  
  ID <- as.character(df[1,]$ALV_ID)
  
  df %>%
    group_by(chr, pos) %>%
    do(get_pairwise_distance(.)) -> distance_sites
  
  piL <- sum(distance_sites$D) # have not yet corrected for genome length
  
  df_return <- data.frame(ALV_ID = ID, piL = piL)
  
  return(df_return)
}

# ===================== Cutoffs for IBV ==============================

# 1% cutoff
no_freq_cut_var %>% filter(freq.var > 0.01) -> var_IBV
var_IBV_monomorphic <- monomorphic(var_IBV, ALV_ID, season, chr, pos, pcr_result)
var_IBV_monomorphic <- filter(var_IBV_monomorphic, !(ref == var & freq.var == 1))
var_IBV_monomorphic <- filter(var_IBV_monomorphic, freq.var < 1)

var_IBV_monomorphic %>%
  group_by(ALV_ID) %>%
  do(get_nt_diversity(.)) -> nt_diversity

meta_snv <- left_join(meta_snv, nt_diversity, by = "ALV_ID")
meta_snv <- mutate(meta_snv, pi_1 = ifelse(pcr_result == "B/VIC", piL / 14772, piL / 14766)) %>% select(-piL)
meta_snv$pi_1[is.na(meta_snv$pi_1)] <- 0

IBV_minority <- filter(var_IBV_monomorphic, freq.var < 0.5)
IBV_minority %>% group_by(ALV_ID) %>% dplyr::summarize(iSNV_1 = length(unique(mutation))) -> IBV_minority_count
meta_snv <- merge(meta_snv, IBV_minority_count, by = "ALV_ID", all.x = TRUE)
meta_snv$iSNV_1[is.na(meta_snv$iSNV_1)] <- 0 # NA means there were no variants; change to 0

# 2% cutoff
no_freq_cut_var %>% filter(freq.var > 0.02) -> var_IBV
var_IBV_monomorphic <- monomorphic(var_IBV, ALV_ID, season, chr, pos, pcr_result)
var_IBV_monomorphic <- filter(var_IBV_monomorphic, !(ref == var & freq.var == 1))
var_IBV_monomorphic <- filter(var_IBV_monomorphic, freq.var < 1)

var_IBV_monomorphic %>%
  group_by(ALV_ID) %>%
  do(get_nt_diversity(.)) -> nt_diversity

meta_snv <- left_join(meta_snv, nt_diversity, by = "ALV_ID")
meta_snv <- mutate(meta_snv, pi_2 = ifelse(pcr_result == "B/VIC", piL / 14772, piL / 14766)) %>% select(-piL)
meta_snv$pi_2[is.na(meta_snv$pi_2)] <- 0

IBV_minority <- filter(var_IBV_monomorphic, freq.var < 0.5)
IBV_minority %>% group_by(ALV_ID) %>% dplyr::summarize(iSNV_2 = length(unique(mutation))) -> IBV_minority_count
meta_snv <- merge(meta_snv, IBV_minority_count, by = "ALV_ID", all.x = TRUE)
meta_snv$iSNV_2[is.na(meta_snv$iSNV_2)] <- 0 # NA means there were no variants; change to 0

# 0.5% cutoff
no_freq_cut_var %>% filter(freq.var > 0.005) -> var_IBV
var_IBV_monomorphic <- monomorphic(var_IBV, ALV_ID, season, chr, pos, pcr_result)
var_IBV_monomorphic <- filter(var_IBV_monomorphic, !(ref == var & freq.var == 1))
var_IBV_monomorphic <- filter(var_IBV_monomorphic, freq.var < 1)

var_IBV_monomorphic %>%
  group_by(ALV_ID) %>%
  do(get_nt_diversity(.)) -> nt_diversity

meta_snv <- left_join(meta_snv, nt_diversity, by = "ALV_ID")
meta_snv <- mutate(meta_snv, pi_05 = ifelse(pcr_result == "B/VIC", piL / 14772, piL / 14766)) %>% select(-piL)
meta_snv$pi_05[is.na(meta_snv$pi_05)] <- 0

IBV_minority <- filter(var_IBV_monomorphic, freq.var < 0.5)
IBV_minority %>% group_by(ALV_ID) %>% dplyr::summarize(iSNV_05 = length(unique(mutation))) -> IBV_minority_count
meta_snv <- merge(meta_snv, IBV_minority_count, by = "ALV_ID", all.x = TRUE)
meta_snv$iSNV_05[is.na(meta_snv$iSNV_05)] <- 0 # NA means there were no variants; change to 0

# ===================== Cutoffs for IAV ==============================

# Genome lengths from reference files:
# if it's H1N1, it's cali09. 13417.
# If it's 10-11 or 11-12 H3N2, it's perth. 13417.
# if it's H3N2 12-13 or 13-14, it's victoria. 13549.
# for 14-15 H3N2, use HK. 13549.

GetGenomeLength <- function(pcr_result, season)
{
  if(pcr_result == "A/H1N1")
  {
    x = 13417
  } else if(season %in% c("10-11", "11-12"))
  {
    x = 13417
  } else
  {
    x = 13549
  }
}
GetGenomeLength.v <- Vectorize(GetGenomeLength)

# 1% cutoff
no_freq_cut_var_IAV %>% filter(freq.var > 0.01) -> var_IAV
var_IAV_monomorphic <- monomorphic(var_IAV, ALV_ID, season, chr, pos, pcr_result)
var_IAV_monomorphic <- filter(var_IAV_monomorphic, !(ref == var & freq.var == 1))
var_IAV_monomorphic <- filter(var_IAV_monomorphic, freq.var < 1)

var_IAV_monomorphic %>%
  group_by(ALV_ID) %>%
  do(get_nt_diversity(.)) -> nt_diversity

meta_IAV_snv <- left_join(meta_IAV_snv, nt_diversity, by = "ALV_ID")
meta_IAV_snv <- mutate(meta_IAV_snv, pi_1 =  piL / GetGenomeLength.v(pcr_result, season)) %>% select(-piL)
meta_IAV_snv$pi_1[is.na(meta_IAV_snv$pi_1)] <- 0

IAV_minority <- filter(var_IAV_monomorphic, freq.var < 0.5)
IAV_minority %>% group_by(ALV_ID) %>% dplyr::summarize(iSNV_1 = length(unique(mutation))) -> IAV_minority_count
meta_IAV_snv <- merge(meta_IAV_snv, IAV_minority_count, by = "ALV_ID", all.x = TRUE)
meta_IAV_snv$iSNV_1[is.na(meta_IAV_snv$iSNV_1)] <- 0 # NA means there were no variants; change to 0

# 2% cutoff
no_freq_cut_var_IAV %>% filter(freq.var > 0.02) -> var_IAV
var_IAV_monomorphic <- monomorphic(var_IAV, ALV_ID, season, chr, pos, pcr_result)
var_IAV_monomorphic <- filter(var_IAV_monomorphic, !(ref == var & freq.var == 1))
var_IAV_monomorphic <- filter(var_IAV_monomorphic, freq.var < 1)

var_IAV_monomorphic %>%
  group_by(ALV_ID) %>%
  do(get_nt_diversity(.)) -> nt_diversity

meta_IAV_snv <- left_join(meta_IAV_snv, nt_diversity, by = "ALV_ID")
meta_IAV_snv <- mutate(meta_IAV_snv, pi_2 =  piL / GetGenomeLength.v(pcr_result, season)) %>% select(-piL)
meta_IAV_snv$pi_2[is.na(meta_IAV_snv$pi_2)] <- 0

IAV_minority <- filter(var_IAV_monomorphic, freq.var < 0.5)
IAV_minority %>% group_by(ALV_ID) %>% dplyr::summarize(iSNV_2 = length(unique(mutation))) -> IAV_minority_count
meta_IAV_snv <- merge(meta_IAV_snv, IAV_minority_count, by = "ALV_ID", all.x = TRUE)
meta_IAV_snv$iSNV_2[is.na(meta_IAV_snv$iSNV_2)] <- 0 # NA means there were no variants; change to 0

# 0.5% cutoff
no_freq_cut_var_IAV %>% filter(freq.var > 0.005) -> var_IAV
var_IAV_monomorphic <- monomorphic(var_IAV, ALV_ID, season, chr, pos, pcr_result)
var_IAV_monomorphic <- filter(var_IAV_monomorphic, !(ref == var & freq.var == 1))
var_IAV_monomorphic <- filter(var_IAV_monomorphic, freq.var < 1)

var_IAV_monomorphic %>%
  group_by(ALV_ID) %>%
  do(get_nt_diversity(.)) -> nt_diversity

meta_IAV_snv <- left_join(meta_IAV_snv, nt_diversity, by = "ALV_ID")
meta_IAV_snv <- mutate(meta_IAV_snv, pi_05 =  piL / GetGenomeLength.v(pcr_result, season)) %>% select(-piL)
meta_IAV_snv$pi_05[is.na(meta_IAV_snv$pi_05)] <- 0

IAV_minority <- filter(var_IAV_monomorphic, freq.var < 0.5)
IAV_minority %>% group_by(ALV_ID) %>% dplyr::summarize(iSNV_05 = length(unique(mutation))) -> IAV_minority_count
meta_IAV_snv <- merge(meta_IAV_snv, IAV_minority_count, by = "ALV_ID", all.x = TRUE)
meta_IAV_snv$iSNV_05[is.na(meta_IAV_snv$iSNV_05)] <- 0 # NA means there were no variants; change to 0

# ======================= Plot and compare ==========================

IAV.hist <- ggplot(no_freq_cut_var_IAV, aes(freq.var)) + 
  geom_histogram(binwidth = 0.0001) + 
  xlim(0, 0.02) + 
  geom_vline(xintercept = 0.02, linetype = "dotted", color = "black", size = 0.5) +
  geom_vline(xintercept = 0.01, linetype = "dotted", color = "black", size = 0.5) + 
  geom_vline(xintercept = 0.005, linetype = "dotted", color = "black", size = 0.5)

IBV.hist <- ggplot(no_freq_cut_var, aes(freq.var)) + 
  geom_histogram(binwidth = 0.0001) + 
  xlim(0, 0.02) + 
  geom_vline(xintercept = 0.02, linetype = "dotted", color = "black", size = 0.5) +
  geom_vline(xintercept = 0.01, linetype = "dotted", color = "black", size = 0.5) + 
  geom_vline(xintercept = 0.005, linetype = "dotted", color = "black", size = 0.5)

# iSNV richness
select(meta_snv, iSNV_2) %>% mutate(type = "IBV", cutoff = "2%") %>% mutate(iSNV = iSNV_2) %>% select(-iSNV_2) -> IBV_snv_2
select(meta_snv, iSNV_1) %>% mutate(type = "IBV", cutoff = "1%") %>% mutate(iSNV = iSNV_1) %>% select(-iSNV_1) -> IBV_snv_1
select(meta_snv, iSNV_05) %>% mutate(type = "IBV", cutoff = "0.5%") %>% mutate(iSNV = iSNV_05) %>% select(-iSNV_05) -> IBV_snv_05

select(meta_IAV_snv, iSNV_2) %>% mutate(type = "IAV", cutoff = "2%") %>% mutate(iSNV = iSNV_2) %>% select(-iSNV_2) -> IAV_snv_2
select(meta_IAV_snv, iSNV_1) %>% mutate(type = "IAV", cutoff = "1%") %>% mutate(iSNV = iSNV_1) %>% select(-iSNV_1) -> IAV_snv_1
select(meta_IAV_snv, iSNV_05) %>% mutate(type = "IAV", cutoff = "0.5%") %>% mutate(iSNV = iSNV_05) %>% select(-iSNV_05) -> IAV_snv_05

compare_cutoffs_iSNV <- rbind(IBV_snv_2, IBV_snv_1, IBV_snv_05, IAV_snv_2, IAV_snv_1, IAV_snv_05)

palette <- wesanderson::wes_palette("FantasticFox1")
dotplot <- ggplot(compare_cutoffs_iSNV, aes(y = iSNV, x = as.factor(type), fill = type, color = type)) +
  geom_dotplot(stackdir = "center", binaxis = 'y', binwidth = 1, dotsize = 0.4) + 
  stat_summary(fun.data = "plot.median", geom = "errorbar", colour = "red", width = 0.95, size = 0.3) +
  xlab("") + 
  ylab("iSNV Per Sample") + 
  theme_bw() + 
  theme(legend.position = "") +
  facet_wrap(~cutoff) + 
  scale_color_manual(name = "", values = palette[c(3,4)]) +
  scale_fill_manual(name = "", values = palette[c(3,4)])

# pi diversity
select(meta_snv, pi_2) %>% mutate(type = "IBV", cutoff = "2%") %>% mutate(pi = pi_2) %>% select(-pi_2) -> IBV_pi_2
select(meta_snv, pi_1) %>% mutate(type = "IBV", cutoff = "1%") %>% mutate(pi = pi_1) %>% select(-pi_1) -> IBV_pi_1
select(meta_snv, pi_05) %>% mutate(type = "IBV", cutoff = "0.5%") %>% mutate(pi = pi_05) %>% select(-pi_05) -> IBV_pi_05

select(meta_IAV_snv, pi_2) %>% mutate(type = "IAV", cutoff = "2%") %>% mutate(pi = pi_2) %>% select(-pi_2) -> IAV_pi_2
select(meta_IAV_snv, pi_1) %>% mutate(type = "IAV", cutoff = "1%") %>% mutate(pi = pi_1) %>% select(-pi_1) -> IAV_pi_1
select(meta_IAV_snv, pi_05) %>% mutate(type = "IAV", cutoff = "0.5%") %>% mutate(pi = pi_05) %>% select(-pi_05) -> IAV_pi_05

compare_cutoffs_pi <- rbind(IBV_pi_2, IBV_pi_1, IBV_pi_05, IAV_pi_2, IAV_pi_1, IAV_pi_05)

palette <- wesanderson::wes_palette("FantasticFox1")
pi.plot <- ggplot(compare_cutoffs_pi, aes(y = pi, x = as.factor(type), fill = type, color = type)) +
  geom_jitter() +
  stat_summary(fun.data = "plot.median", geom = "errorbar", colour = "red", width = 0.95, size = 0.3) +
  xlab("") + 
  ylab("Nucleotide diversity (pi)") + 
  theme_bw() + 
  theme(legend.position = "") +
  facet_wrap(~cutoff) + 
  scale_color_manual(name = "", values = palette[c(3,4)]) +
  scale_fill_manual(name = "", values = palette[c(3,4)])

# statistical tests: iSNV
wilcox.test(IBV_snv_2$iSNV, IAV_snv_2$iSNV)
wilcox.test(IBV_snv_1$iSNV, IAV_snv_1$iSNV)
wilcox.test(IBV_snv_05$iSNV, IAV_snv_05$iSNV)

# statistical tests: pi
wilcox.test(IBV_pi_2$pi, IAV_pi_2$pi)
wilcox.test(IBV_pi_1$pi, IAV_pi_1$pi)
wilcox.test(IBV_pi_05$pi, IAV_pi_05$pi)

