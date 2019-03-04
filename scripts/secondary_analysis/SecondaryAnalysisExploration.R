
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
qual <- read_csv("../data/processed/qual.snv.csv")
no_freq_cut <- read_csv("../data/processed/no_freq_cut.qual.snv.csv")
qual_not_coll <- read_csv("../data/processed/qual.not.collapsed.snv.csv")
chrs <- read.csv("../data/metadata/ibv.segs.csv", stringsAsFactors = T) # May need to check that this is correct.

# =========== How severe was the 2% cutoff? ================

nocut_isnv <- filter(no_freq_cut, freq.var < 0.98) # 8693 variants
std_isnv <- filter(qual, freq.var < 0.98) # 136 variants

ggplot(std_isnv, aes(freq.var)) + geom_histogram(binwidth = 0.01)
ggplot(nocut_isnv, aes(freq.var)) + geom_histogram(binwidth = 0.0001) + xlim(0, 0.02) + geom_vline(xintercept = 0.02, linetype = "dotted", color = "black", size = 1) # Mostly with very low frequencies, even though seen in both samples. Almost all below 0.5%.

# =========== Histogram of iSNV counts per sample ================

titer.plot <- ggplot(meta, aes(x=as.factor(DPSO), y=genome_copy_per_ul)) + geom_boxplot(notch = FALSE) + ylab(expression(paste(Genomes,"/" ,mu,L))) + scale_y_log10() + xlab("Days Post Symptom Onset") + theme_bw()
ggsave(plot = titer.plot, filename = "../results/plots/TiterByDPSO.jpg", device = "jpeg")

# =========== iSNV by DPSO ================

min.qual.o <- subset(qual, freq.var < 0.98)

#write_csv(min.qual.o, "../data/processed/min.qual.o.csv")

min.qual <- subset(min.qual.o, !(ENROLLID %in% c("50425"))) # remove mixed infections
min.count.sample.o <- min.qual.o %>% group_by(ALV_ID) %>% dplyr::summarize(iSNV = length(unique(mutation)), HA_iSNV = length(which(chr=="HA"))) 
snv_qual_meta.o <- subset(meta, coverage_qualified == TRUE)
snv_qual_meta.o <- merge(snv_qual_meta.o, min.count.sample.o, by="ALV_ID", all.x=T)
snv_qual_meta.o$iSNV[is.na(snv_qual_meta.o$iSNV)] <- 0 # these are the ones with no diversity
snv_qual_meta.o$HA_iSNV[is.na(snv_qual_meta.o$HA_iSNV)] <- 0

write_csv(snv_qual_meta.o, "../data/processed/snv_qual_meta.o.csv")

isnv_by_day.p <- ggplot(snv_qual_meta.o, aes(x = as.factor(DPSO), y = iSNV)) + geom_boxplot() + xlab("Day Post Symptom Onset") + geom_jitter() + theme_bw() #+ ylim(0, 10)
ggsave(plot = isnv_by_day.p, filename = "../results/plots/SNVbyDPSO.jpg", device = "jpeg")

# =========== Histogram of iSNV counts per sample ================

isnv.per.sample <- ggplot(snv_qual_meta.o, aes(x = iSNV)) + geom_histogram(binwidth = 1, fill = "maroon") + xlab("Number of iSNV (Bin Width = 1)") + ylab("Number of samples") + theme_bw()
ggsave(plot = isnv.per.sample, filename = "../results/plots/SNVperSample.jpg", device = "jpeg")

# =========== iSNV counts per sample by genome copy number ================

snv_qual_meta.o.above_E5 <- filter(snv_qual_meta.o, genome_copy_per_ul > 100000)
isnv.per.sample.aboveE5 <- ggplot(snv_qual_meta.o.above_E5, aes(x = iSNV)) + geom_histogram(binwidth = 1, color = "white", fill = "maroon") + xlab("Number of iSNV") + ylab("Number of samples") + ggtitle("iSNV in samples above E5 copies/uL") + theme_bw()

snv_by_copynum <- ggplot(snv_qual_meta.o, aes(x = log(genome_copy_per_ul, 10), y = iSNV)) + geom_point(shape = 19) + xlab("Log (base 10) of genome copies/uL") + geom_vline(xintercept = 5, linetype = "dotted", color = "maroon", size = 1.5) + theme_bw()
snv_by_gc <- lm(data = snv_qual_meta.o, formula = iSNV ~ log(genome_copy_per_ul, 10))
summary(snv_by_gc)
ggsave(plot = snv_by_copynum, filename = "../results/plots/SNVbyCopyNumber.jpg", device = "jpeg")

# =========== iSNV counts by vaccination status ================

plot.median <- function(x) {
  m <- median(x)
  c(y = m, ymin = m, ymax = m)
}

isnv_by_vaccination <- ggplot(snv_qual_meta.o, aes(y = iSNV, x = as.factor(vaccination_status))) +
  geom_dotplot(stackdir = "center", binaxis = 'y', binwidth = 1, dotsize = 0.5) +
  stat_summary(fun.data = "plot.median", geom = "errorbar", colour = "red", width = 0.95, size = 0.3) +
  scale_x_discrete(labels = c("Not Vaccinated", "Vaccinated")) + xlab("") + theme_bw()
ggsave(plot = isnv_by_vaccination, filename = "../results/plots/SNVbyVaccinationStatus.jpg", device = "jpeg")

# =========== iSNV across genome segments ================

# Using the list excluding the outlier individual
min.qual$chr <- factor(min.qual$chr, levels = rev(c("PB2","PB1","PA","HA","NP","NR","M1","NS")))

chrs$chr <- factor(chrs$chr, levels = levels(min.qual$chr)) # set factors on the is meta data

genome_loc.p <- ggplot(min.qual, aes(x = pos, y = chr)) +
  geom_point(aes(color = class_factor), shape = 108, size = 5 )+
  geom_segment(data = chrs, aes(x = start, y = chr, xend = stop, yend = chr)) +
  ylab("") +
  xlab("") +
  scale_color_manual(name = "", values = palette[c(4,3)]) +
  theme(axis.ticks = element_blank(),
        axis.line.x = element_blank(), axis.line.y = element_blank()) +
  scale_x_continuous(breaks = c()) +
  theme(legend.position = "right") + theme_minimal()

ggsave(plot = genome_loc.p, filename = "../results/plots/SNVbyGenomeLocation.jpg", device = "jpeg")

# add multiple iSNV
# In how many people (ENROLLID) is the muation found. It would have to be minor in both 
# people to count in this plot. The chr pos var is used so that called and infered variants 
# of the same var are counted together.

min.count <- min.qual %>% group_by(chr,pos,var,pcr_result,season) %>%
  dplyr::summarize(counts=length(unique(ENROLLID))) 
min.qual <- mutate(min.qual,scratch = paste(chr,pos,chr,var,season,pcr_result,sep="."))
min.count <- mutate(min.count,scratch = paste(chr,pos,chr,var,season,pcr_result,sep="."))
multiple <- subset(min.qual,scratch %in% min.count$scratch[min.count$counts>1])

min.qual <- subset(min.qual,select=-c(scratch))
min.count <- subset(min.count,select=-c(scratch))

genome_loc.p.dots <- genome_loc.p + geom_point(data = multiple, aes(x=pos, y=as.numeric(chr) + 0.3, color=class_factor), size=1, shape=6)

# =========== Frequency by Nonsyn/Syn ================

min.qual.low <- subset(min.qual, freq.var < 0.5) # Exclude outlier individual, and only look at those below 50%.

freq_hist.p <- ggplot(min.qual.low, aes(x=freq.var,fill=class_factor)) + geom_histogram(color="white",binwidth=.05,position=position_dodge(),boundary = 0.02) +
  xlab("Frequency") + ylab("iSNV") +
  scale_fill_manual(name="" ,values=palette[c(4,3)] )+
  theme(legend.position = c(0.5, 0.5)) + theme_classic()

ggsave(plot = freq_hist.p, filename = "../results/plots/SNVbyFrequency.jpg", device = "jpeg")

# =========== Distribution of sampling times infections with home and clinic samples ================

meta_homeAndClinic <- filter(meta_short, home_spec == 1)

meta_homeAndClinic %>% mutate(DPS1 = home_spec_date-onset, DPS2 = collect-onset) -> meta_homeAndClinic
meta_homeAndClinic <- meta_homeAndClinic[order(meta_homeAndClinic$DPS1,meta_homeAndClinic$DPS2,decreasing = T),]
meta_homeAndClinic <- mutate(meta_homeAndClinic, DPS2 = ifelse(DPS1 == DPS2, yes = DPS2 + 0.3, no = DPS2)) 
meta_homeAndClinic$sort_order <- 1:nrow(meta_homeAndClinic)

sampling_distribution <- ggplot(meta_homeAndClinic, aes(x = DPS1, xend = DPS2, y = sort_order, yend = sort_order)) + geom_segment(color = palette[1]) + ylab("") +
  xlab("Day post symptom onset") + 
  theme(axis.line.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_point(aes(y = sort_order, x = DPS1), color=palette[1])+
  geom_point(aes(y = sort_order, x = DPS2), color=palette[1]) + 
  scale_x_continuous(breaks = -2:6) + theme_classic()

ggsave(plot = sampling_distribution, filename = "../results/plots/SampleCollectionTimes.jpg", device = "jpeg")

# =========== Concordance between sequencing replicates ================

# Get all variants from sequencing replicates of the same sample (whether it's a home or clinic sample). Draw variants from qual_not_coll.

min.qual_not_coll.o <- subset(qual_not_coll, freq.var < 0.98)
min.qual_not_coll <- min.qual_not_coll.o
#min.qual_not_coll <- subset(min.qual_not_coll.o, !(ENROLLID %in% c("50425"))) # remove mixed infections

min.qual_not_coll.1 <- filter(min.qual_not_coll, SampleNumber == SeqSampleNumber1)
min.qual_not_coll.2 <- filter(min.qual_not_coll, SampleNumber == SeqSampleNumber2)
min.qual_not_coll.1 <- mutate(min.qual_not_coll.1, mut_id = paste0(mutation, "_", ALV_ID))
min.qual_not_coll.2 <- mutate(min.qual_not_coll.2, mut_id = paste0(mutation, "_", ALV_ID))
merged <- merge(min.qual_not_coll.1, min.qual_not_coll.2, by = "mut_id")

replicate_concordance <- ggplot(data = merged, aes(x = freq.var.x, y = freq.var.y)) +
  geom_point(aes(color = as.factor(floor(log10(genome_copy_per_ul.x))))) +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  scale_y_continuous(limits = c(0,0.5)) +
  scale_x_continuous(limits = c(0,0.5)) +
  xlab("Frequency in replicate 1") + ylab("Frequency in replicate 2") +
  scale_color_manual(name = "Log(copies/ul)", values = palette[c(4,3)]) + theme_bw()

ggsave(plot = replicate_concordance, filename = "../results/plots/ReplicateConcordance.jpg", device = "jpeg")

# =========== Concordance between frequency in home and clinic isolates ================ 

### Is this just because there are so few variants? ###

# Draw variants from min.qual.low (collapsed version, removed mixed infection, 2-50%).

meta_homeAndClinic <- filter(meta_short, home_spec == 1) # only those from which we have both home and clinic samples
meta_homeAndClinic <- mutate(meta_homeAndClinic, within_host_time = collect - home_spec_date)
#meta_homeAndClinic <- filter(meta_homeAndClinic, within_host_time == 0) # filter to specify within-host time
variants_home <- filter(min.qual.low, ALV_ID %in% meta_homeAndClinic$home_spec_accn)
variants_clinic <- filter(min.qual.low, ALV_ID %in% meta_homeAndClinic$SPECID)

variants_home %>% group_by(ALV_ID) %>% filter(mutation %in% variants_clinic$mutation) -> home
variants_clinic %>% group_by(ALV_ID) %>% select(mutation, freq.var, ALV_ID) -> clinic
merge(home, clinic, by = "mutation") %>% select(mutation, ALV_ID.x, SampleNumber, freq.var.x, ALV_ID.y, freq.var.y) -> home_vs_clinic

home_clinic_concordance <- ggplot(data = home_vs_clinic, aes(x=freq.var.x,y=freq.var.y)) + geom_point()+
  xlab("Frequency in home isolate") + ylab("Frequency in clinic isolate") + 
  geom_abline(slope=1,intercept = 0,lty=2)+
  scale_x_continuous(limits = c(0,0.5))+scale_y_continuous(limits = c(0,0.5)) + theme_bw()

# =========== Changes in frequency across paired home and clinic isolates ================ 

# for each variant, get frequency from home and clinic sample, and within host time.
# min.qual.o for all verified SNV from 2-98%
# min.qual for all verified SNV from 2-98%. Excludes outlier infection.
# min.qual.low for all verified SNV from 2-50%. Excludes outlier infection.

variants_for_paired_sample_analysis <- min.qual

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

# by enroll ID, give it the withinhost time from the metadata
intra <- merge(intra, select(meta_homeAndClinic, ENROLLID, within_host_time), by = "ENROLLID")

intra %>% mutate(Endpoint = "Persistent") %>%
  mutate(Endpoint = if_else(freq1 == 0,"Arisen", Endpoint)) %>%
  mutate(Endpoint = if_else(freq2 == 0,"Lost", Endpoint)) -> intra

intra$Endpoint <- factor(intra$Endpoint, levels = c("Persistent","Arisen","Lost"), ordered = T)

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

variants_HA_nonsyn <- filter(min.qual.o, class_factor == "Nonsynonymous" & chr == "HA") # only six! None in antigenic sites, by inspection.
variants_HA_nonsyn <- mutate(variants_HA_nonsyn, is_antigenic = ifelse(pcr_result == "B/VIC", pos %in% VIC_sites, pos %in% YAM_sites))

antigenic.plot <- ggplot(variants_HA_nonsyn, aes(y = freq.var, x = is_antigenic, color = class_factor)) +
  geom_quasirandom(varwidth = TRUE) +
  ylab("iSNV frequency") + xlab(label = "") +
  scale_x_discrete(labels = c("Antigenic site","Nonantigenic site")) +
  scale_y_continuous(limits=c(0,1)) + scale_color_manual(name="Mutation Class",values=palette[c(1,4)]) + ggtitle("Mutations in HA")

