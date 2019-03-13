
### Project: IBV
### Purpose: Some initial analysis of iSNVs. Explore the data.

# =========== Import packages, load data ================

setwd("/Users/avalesano/Documents/MSTP/LauringLab/Host_level_IBV_evolution/scripts/")

library(tidyverse)
library(magrittr)
library(wesanderson)
library(ggbeeswarm)
palette <- wesanderson::wes_palette("FantasticFox1")

meta <- read_csv("../data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv")
meta_short <- read_csv("../data/metadata/flu_b_2010_2017_v4.csv")
metadata_by_seq <- read_csv("../data/metadata/Metadata_By_Seq_withALVID.csv")
quality_var <- read_csv("../data/processed/qual.snv.csv")
no_freq_cut_var <- read_csv("../data/processed/no_freq_cut.qual.snv.csv")
qual_not_coll <- read_csv("../data/processed/qual.not.collapsed.snv.csv")
chrs <- read.csv("../data/metadata/ibv.segs.csv", stringsAsFactors = T) # May need to check that this is correct.
IAV_snv_qual_meta <- read_csv("../data/processed/IAV_snv_qual_meta.csv")

# =========== How severe was the 2% cutoff? ================

nocut_isnv <- filter(no_freq_cut_var, freq.var < 0.98) # 8693 variants
isnv <- filter(quality_var, freq.var < 0.98) # 136 variants
isnv_nomixed <- filter(isnv, !(ENROLLID %in% c("50425")))
isnv_minority <- filter(quality_var, freq.var < 0.5) # 68 variants below 50%

ggplot(isnv, aes(freq.var)) + geom_histogram(binwidth = 0.01)
ggplot(nocut_isnv, aes(freq.var)) + geom_histogram(binwidth = 0.0001) + xlim(0, 0.02) + geom_vline(xintercept = 0.02, linetype = "dotted", color = "black", size = 1) # Mostly with very low frequencies, even though seen in both samples. Almost all below 0.5%.

# =========== Histogram of iSNV counts per sample ================

titer.plot <- ggplot(meta, aes(x = as.factor(DPSO), y = genome_copy_per_ul)) + geom_boxplot(notch = TRUE) + ylab(expression(paste(Genomes,"/" ,mu,L))) + scale_y_log10() + xlab("Days Post Symptom Onset") + theme_bw()
titer.plot.scatter <- ggplot(meta, aes(x = DPSO, y = genome_copy_per_ul)) + geom_point() + ylab(expression(paste(Genomes,"/" , mu, L))) + scale_y_log10() + xlab("Days Post Symptom Onset") + theme_bw()

ggsave(plot = titer.plot, filename = "../results/plots/TiterByDPSO_box.jpg", device = "jpeg")
ggsave(plot = titer.plot.scatter, filename = "../results/plots/TiterByDPSO_scatter.jpg", device = "jpeg")

# =========== iSNV by DPSO ================

isnv_minority_nomixed <- filter(isnv_minority, !(ENROLLID %in% c("50425"))) # remove mixed infections

isnv_minority %>% group_by(ALV_ID) %>% dplyr::summarize(iSNV = length(unique(mutation)), HA_iSNV = length(which(chr == "HA"))) -> isnv_minority_count
meta_snv <- filter(meta, coverage_qualified == TRUE)
meta_snv <- merge(meta_snv, isnv_minority_count, by="ALV_ID", all.x = TRUE)
meta_snv$iSNV[is.na(meta_snv$iSNV)] <- 0 # NA means there were no variants; change to 0
meta_snv$HA_iSNV[is.na(meta_snv$HA_iSNV)] <- 0

write_csv(meta_snv, "../data/processed/meta_snv.csv")

isnv_by_day.plot <- ggplot(meta_snv, aes(x = as.factor(DPSO), y = iSNV)) + geom_boxplot(outlier.shape = NA, notch = FALSE) + xlab("Day Post Symptom Onset") + geom_jitter(width = 0.3, height = 0.1) + theme_bw()
ggsave(plot = isnv_by_day.plot, filename = "../results/plots/SNVbyDPSO.jpg", device = "jpeg")

# =========== Histogram of iSNV counts per sample ================

isnv.per.sample <- ggplot(meta_snv, aes(x = iSNV)) + geom_histogram(binwidth = 1, fill = palette[5], color = "white") + xlab("Number of iSNV (Bin Width = 1)") + ylab("Number of samples") + theme_bw()
ggsave(plot = isnv.per.sample, filename = "../results/plots/SNVperSample.jpg", device = "jpeg")

# Compare to IAV data
IBV_snv <- select(meta_snv, iSNV)
IBV_snv <- mutate(IBV_snv, type = "Influenza B")
IAV_snv <- select(IAV_snv_qual_meta, iSNV)
IAV_snv <- mutate(IAV_snv, type = "Influenza A")
cfIAV_data <- rbind(IBV_snv, IAV_snv)

cfIAV.plot.identity <- ggplot(cfIAV_data, aes(x = iSNV, fill = type)) + geom_histogram(binwidth = 1, alpha = 0.5, color = "white", position = "identity") + scale_fill_manual(name = "Type", values = c(palette[3], palette[5])) + xlab("Number of iSNV") + ylab("Number of samples") + theme_bw()
cfIAV.plot.dodge <- ggplot(cfIAV_data, aes(x = iSNV, fill = type)) + geom_histogram(binwidth = 1, alpha = 1, color = "white", position = "dodge") + scale_fill_manual(name = "Type", values = c(palette[3], palette[5])) + xlab("Number of iSNV") + ylab("Number of samples") + theme_bw()

ggsave(plot = IAV.isnv.per.sample, filename = "../results/plots/SNVperSample_cfIAV_identity.jpg", device = "jpeg")
ggsave(plot = IAV.isnv.per.sample, filename = "../results/plots/SNVperSample_cfIAV_dodge.jpg", device = "jpeg")

# =========== iSNV counts per sample by genome copy number ================

snv_by_copynum <- ggplot(meta_snv, aes(x = log(genome_copy_per_ul, 10), y = iSNV)) + geom_point(shape = 19) + xlab("Log (base 10) of genome copies/uL") + geom_vline(xintercept = 5, linetype = "dotted", color = palette[5], size = 1.5) + theme_bw()
snv_by_gc <- lm(data = meta_snv, formula = iSNV ~ log(genome_copy_per_ul, 10))
summary(snv_by_gc)
ggsave(plot = snv_by_copynum, filename = "../results/plots/SNVbyCopyNumber.jpg", device = "jpeg")

# =========== iSNV counts by vaccination status ================

plot.median <- function(x) {
  m <- median(x)
  c(y = m, ymin = m, ymax = m)
}

isnv_by_vaccination <- ggplot(meta_snv, aes(y = iSNV, x = as.factor(vaccination_status))) +
  geom_dotplot(stackdir = "center", binaxis = 'y', binwidth = 1, dotsize = 0.2) +
  stat_summary(fun.data = "plot.median", geom = "errorbar", colour = "red", width = 0.95, size = 0.3) +
  scale_x_discrete(labels = c("Not Vaccinated", "Vaccinated")) + xlab("") + theme_bw()

ggsave(plot = isnv_by_vaccination, filename = "../results/plots/SNVbyVaccinationStatus.jpg", device = "jpeg")

# =========== iSNV across genome segments ================

isnv_minority_nomixed$chr <- factor(isnv_minority_nomixed$chr, levels = rev(c("PB2","PB1","PA","HA","NP","NR","M1","NS")))
chrs$chr <- factor(chrs$chr, levels = levels(isnv_minority_nomixed$chr))

# I think this IS correct. Just excludes the mixed infection. 68 minority variants in all, 40 if you exclude the mixed infection. Then 4 HA NonSyn variants (all non antigenic).
genome_location.plot <- ggplot(isnv_minority_nomixed, aes(x = pos, y = chr)) +
  geom_point(aes(color = class_factor), shape = 108, size = 5) +
  geom_segment(data = chrs, aes(x = start, y = chr, xend = stop, yend = chr)) +
  ylab("") +
  xlab("") +
  scale_color_manual(name = "", values = palette[c(4,3)]) +
  theme(axis.ticks = element_blank(),
        axis.line.x = element_blank(), axis.line.y = element_blank()) +
  scale_x_continuous(breaks = c()) +
  theme(legend.position = "right") + theme_minimal() + theme(panel.grid.major = element_line(colour = "white"))

ggsave(plot = genome_location.plot, filename = "../results/plots/SNVbyGenomeLocation.jpg", device = "jpeg")

# Mark minority SNV that are found in multiple individuals.

isnv_minority_nomixed %>% group_by(chr, pos, var, pcr_result, season) %>% dplyr::summarize(counts = length(unique(ENROLLID))) -> isnv_minority_nomixed_count
isnv_minority_nomixed <- mutate(isnv_minority_nomixed, scratch = paste(chr, pos, chr, var, season, pcr_result, sep = "."))
isnv_minority_nomixed_count <- mutate(isnv_minority_nomixed_count, scratch = paste(chr, pos, chr, var, season, pcr_result, sep = "."))
multiple <- subset(isnv_minority_nomixed, scratch %in% isnv_minority_nomixed_count$scratch[isnv_minority_nomixed_count$counts > 1])

genome_location.plot.multiple <- genome_location.plot + geom_point(data = multiple, aes(x = pos, y = as.numeric(chr) + 0.2, color = class_factor), size = 1, shape = 25)

# =========== Frequency by Nonsyn/Syn ================

freq_histogram <- ggplot(isnv_minority_nomixed, aes(x = freq.var, fill = class_factor)) + geom_histogram(color = "white", binwidth = 0.05, position = "dodge", boundary = 0.02) +
  xlab("Frequency") + ylab("Number of iSNV") +
  scale_fill_manual(name = "" , values = palette[c(4,3)])+
  theme(legend.position = c(0.5, 0.5)) + theme_classic()

ggsave(plot = freq_histogram, filename = "../results/plots/SNVbyFrequency.jpg", device = "jpeg")

# =========== Distribution of sampling times infections with home and clinic samples ================

meta_homeAndClinic <- filter(meta_short, home_spec == 1) # This is all flu B samples from HIVE. Need to filter by those that which we sequenced.
IDs_cov_qualified <- filter(meta, coverage_qualified == TRUE)
meta_homeAndClinic <- filter(meta_homeAndClinic, SPECID %in% IDs_cov_qualified$ALV_ID & home_spec_accn %in% IDs_cov_qualified$ALV_ID)

# Now calculate time between sampling.
meta_homeAndClinic %>% mutate(DPS1 = home_spec_date-onset, DPS2 = collect-onset) -> meta_homeAndClinic
meta_homeAndClinic <- meta_homeAndClinic[order(meta_homeAndClinic$DPS1, meta_homeAndClinic$DPS2, decreasing = TRUE),]
meta_homeAndClinic <- mutate(meta_homeAndClinic, DPS2 = ifelse(DPS1 == DPS2, yes = DPS2 + 0.3, no = DPS2)) 
meta_homeAndClinic$sort_order <- 1:nrow(meta_homeAndClinic)

sampling_distribution <- ggplot(meta_homeAndClinic, aes(x = DPS1, xend = DPS2, y = sort_order, yend = sort_order)) + geom_segment(color = palette[3]) + ylab("") +
  xlab("Day post symptom onset") +
  geom_point(aes(y = sort_order, x = DPS1), color = "black") +
  geom_point(aes(y = sort_order, x = DPS2), color = "black") + 
  scale_x_continuous(breaks = -2:6) + theme_classic() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())

ggsave(plot = sampling_distribution, filename = "../results/plots/SampleCollectionTimes.jpg", device = "jpeg")

# =========== Concordance between sequencing replicates ================

# Get all variants from sequencing replicates of the same sample (whether it's a home or clinic sample). Draw variants from qual_not_coll.

isnv_not_collapsed <- filter(qual_not_coll, freq.var < 0.98)
#isnv_not_collapsed <- filter(isnv_not_collapsed, !(ENROLLID %in% c("50425"))) # remove mixed infection

isnv_not_collapsed_rep1 <- filter(isnv_not_collapsed, SampleNumber == SeqSampleNumber1)
isnv_not_collapsed_rep2 <- filter(isnv_not_collapsed, SampleNumber == SeqSampleNumber2)
isnv_not_collapsed_rep1 <- mutate(isnv_not_collapsed_rep1, mut_id = paste0(mutation, "_", ALV_ID))
isnv_not_collapsed_rep2 <- mutate(isnv_not_collapsed_rep2, mut_id = paste0(mutation, "_", ALV_ID))
merged <- merge(isnv_not_collapsed_rep1, isnv_not_collapsed_rep2, by = "mut_id")

merged_nomixed <- filter(merged, !(ALV_ID.x %in% c("HS1875", "MH10536")))

replicate_concordance_all <- ggplot(data = merged, aes(x = freq.var.x, y = freq.var.y)) +
  geom_point(aes(color = as.factor(floor(log10(genome_copy_per_ul.x))))) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  scale_y_continuous(limits = c(0,0.5)) +
  scale_x_continuous(limits = c(0,0.5)) +
  xlab("Frequency in replicate 1") + ylab("Frequency in replicate 2") +
  scale_color_manual(name = "Log(copies/ul)", values = palette[c(4,3)]) + theme_bw()

replicate_concordance_no_mixed <- ggplot(data = merged_nomixed, aes(x = freq.var.x, y = freq.var.y)) +
  geom_point(aes(color = as.factor(floor(log10(genome_copy_per_ul.x))))) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  scale_y_continuous(limits = c(0,0.5)) +
  scale_x_continuous(limits = c(0,0.5)) +
  xlab("Frequency in replicate 1") + ylab("Frequency in replicate 2") +
  scale_color_manual(name = "Log(copies/ul)", values = palette[c(4,3)]) + theme_bw()

ggsave(plot = replicate_concordance_all, filename = "../results/plots/ReplicateConcordance_All.jpg", device = "jpeg")
ggsave(plot = replicate_concordance_no_mixed, filename = "../results/plots/ReplicateConcordance_NoMixed.jpg", device = "jpeg")

# =========== Concordance between frequency in home and clinic isolates ================ 

### Is this just because there are so few variants? Probably. Few variants in each sample, few samples that have home and clinic sample on the same day.

meta_homeAndClinic <- filter(meta_short, home_spec == 1) # only those from which we have both home and clinic samples
meta_homeAndClinic <- mutate(meta_homeAndClinic, within_host_time = collect - home_spec_date)
#meta_homeAndClinic <- filter(meta_homeAndClinic, within_host_time == 0) # filter to specify within-host time
variants_home <- filter(isnv_minority_nomixed, ALV_ID %in% meta_homeAndClinic$home_spec_accn)
variants_clinic <- filter(isnv_minority_nomixed, ALV_ID %in% meta_homeAndClinic$SPECID)

variants_home %>% group_by(ALV_ID) %>% filter(mutation %in% variants_clinic$mutation) -> home
variants_clinic %>% group_by(ALV_ID) %>% select(mutation, freq.var, ALV_ID) -> clinic
merge(home, clinic, by = "mutation") %>% select(mutation, ALV_ID.x, SampleNumber, freq.var.x, ALV_ID.y, freq.var.y) -> home_vs_clinic

home_clinic_concordance <- ggplot(data = home_vs_clinic, aes(x=freq.var.x,y=freq.var.y)) + geom_point()+
  xlab("Frequency in home isolate") + ylab("Frequency in clinic isolate") + 
  geom_abline(slope=1,intercept = 0,lty=2)+
  scale_x_continuous(limits = c(0,0.5))+scale_y_continuous(limits = c(0,0.5)) + theme_bw()

# =========== Changes in frequency across paired home and clinic isolates ================ 

variants_for_paired_sample_analysis <- isnv_nomixed

meta_homeAndClinic <- filter(meta_short, home_spec == 1) # only those from which we have both home and clinic samples
meta_homeAndClinic <- mutate(meta_homeAndClinic, within_host_time = collect - home_spec_date)
variants_home <- filter(variants_for_paired_sample_analysis, ALV_ID %in% meta_homeAndClinic$home_spec_accn)
variants_clinic <- filter(variants_for_paired_sample_analysis, ALV_ID %in% meta_homeAndClinic$SPECID)
variants_home <- mutate(variants_home, mut_id = paste0(mutation, "_", ENROLLID))
variants_clinic <- mutate(variants_clinic, mut_id = paste0(mutation, "_", ENROLLID))

intra <- merge(variants_home, variants_clinic, by = "mut_id", all = TRUE)
intra <- mutate(intra, freq1 = ifelse(is.na(freq.var.x), 0, freq.var.x))
intra <- mutate(intra, freq2 = ifelse(is.na(freq.var.y), 0, freq.var.y))
intra <- mutate(intra, ENROLLID = ifelse(!is.na(ENROLLID.x), ENROLLID.x, ENROLLID.y))

# By enroll ID, give it the withinhost time from the metadata.
intra <- merge(intra, select(meta_homeAndClinic, ENROLLID, within_host_time), by = "ENROLLID")

intra %>% mutate(Endpoint = "Persistent") %>%
  mutate(Endpoint = if_else(freq1 == 0,"Arisen", Endpoint)) %>%
  mutate(Endpoint = if_else(freq2 == 0,"Lost", Endpoint)) -> intra

intra$Endpoint <- factor(intra$Endpoint, levels = c("Persistent","Arisen","Lost"), ordered = TRUE)

# Maybe frame by syn/non-syn?
intra.plot <- ggplot(intra, aes(x = as.factor(within_host_time),
                             y = freq2 - freq1,
                             fill = Endpoint)) +
  geom_quasirandom(pch=21,color='black',size=2) +
  scale_fill_manual(values=palette[c(1,3,5)],name="") +
  xlab("Time within host (days)") + ylab("Change in frequency") + ggtitle("SNV across paired home and clinic samples") + theme_bw()

ggsave(plot = intra.plot, filename = "../results/plots/IntraHostSNV_Longitudinal.jpg", device = "jpeg")

# =========== SNVs in HA antigenic vs non-antigenic sites ================ 

antigenic_sites <- read_csv("../data/processed/antigenic_positions.csv")

VIC_start <- filter(antigenic_sites, Lineage == "B/VIC")$Start
VIC_end <- filter(antigenic_sites, Lineage == "B/VIC")$End
VIC_sites <- c()
for(r in 1:length(VIC_start))
{
  range <- c(VIC_start[r]:VIC_end[r])
  VIC_sites <- c(VIC_sites, range)
}

YAM_start <- filter(antigenic_sites, Lineage == "B/YAM")$Start
YAM_end <- filter(antigenic_sites, Lineage == "B/YAM")$End
YAM_sites <- c()
for(r in 1:length(YAM_start))
{
  range <- c(YAM_start[r]:YAM_end[r])
  YAM_sites <- c(YAM_sites, range)
}

variants_HA_nonsyn <- filter(isnv_minority_nomixed, class_factor == "Nonsynonymous" & chr == "HA") # only four! None in antigenic sites, by inspection.
variants_HA_nonsyn <- mutate(variants_HA_nonsyn, is_antigenic = ifelse(pcr_result == "B/VIC", pos %in% VIC_sites, pos %in% YAM_sites))

#antigenic.plot <- ggplot(variants_HA_nonsyn, aes(y = freq.var, x = is_antigenic, color = class_factor)) +
  geom_quasirandom(varwidth = TRUE) +
  ylab("iSNV Frequency") + xlab(label = "") +
  scale_x_discrete(labels = c("Antigenic site","Nonantigenic site")) +
  scale_y_continuous(limits = c(0,1)) + scale_color_manual(name = "Mutation Class", values = palette[c(5,4)]) + ggtitle("Mutations in HA") + theme_bw()

