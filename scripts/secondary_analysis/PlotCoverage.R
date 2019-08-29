
### Project: IBV intra-host
### Purpose: Plot coverage over the genome

# ======================== Import packages ================================

require(tidyverse)
require(HIVEr)

# ======================== Read in data, filter by lineage ================================

vic.cov.files = c("data/raw/2445_VIC/all.coverage.csv", "data/raw/2446_VIC/all.coverage.csv")
yam.cov.files = c("data/raw/2445_YAM/all.coverage.csv", "data/raw/2446_YAM/all.coverage.csv")

# read_rbind is in HIVEr
vic.cov <- read_rbind(vic.cov.files,
                3,
                cols = cols(
                  chr = col_character(),
                  chr.pos = col_integer(),
                  coverage = col_integer(),
                  concat.pos = col_integer(),
                  Id = col_character()
                ))

yam.cov <- read_rbind(yam.cov.files,
                    3,
                    cols = cols(
                      chr = col_character(),
                      chr.pos = col_integer(),
                      coverage = col_integer(),
                      concat.pos = col_integer(),
                      Id = col_character()
                    ))

# Filter CSVs by lineage type.
metadata <- read_csv("data/metadata/IBV_Sequenced_WithMetadata_InferredTypes.csv")
VIC_IDs <- filter(metadata, final_result == "B_VIC")
YAM_IDs <- filter(metadata, final_result == "B_YAM")

vic.cov <- filter(vic.cov, Id %in% VIC_IDs$SampleNumber)
yam.cov <- filter(yam.cov, Id %in% YAM_IDs$SampleNumber)

# change chr names for NA_ and NA_clone1 to NR here.
vic.cov <- mutate(vic.cov, chr = gsub("NA_", "NR", chr))
yam.cov <- mutate(yam.cov, chr = gsub("NA_clone1", "NR", chr))


# --------------------------------- Functions --------------------------------
#   Read in the csv files used in the data analysis below
# -----------------------------------------------------------------------------
# coverage

slide <- function(cov.df, setup.df)
{
  coverage = rep(NA, nrow(setup.df))
  for(i in 1:nrow(setup.df))
  {
    s = setup.df$starts[i]
    e = setup.df$ends[i]
    subset(cov.df, concat.pos >= s & concat.pos < e, select = c(coverage)) -> position
    mean(position$coverage) -> coverage[i]
  }
  out <- data.frame(mean = coverage, concat.pos = setup.df$concat.pos, chr = setup.df$chr)
  out$Id = unique(cov.df$Id)
  out$run = unique(cov.df$run)
  return(out)
}

cov_plot <- function(cov.df, title)
{
  cov.df %>% group_by(chr) %>% summarize(first = min(concat.pos), last = max(concat.pos)) %>% plyr::adply(1,function(x) data.frame(starts = seq(x$first,x$last,by=100))) %>% mutate(ends = ifelse(starts + 200 < last, starts + 200, last)) -> setup
  setup %>% select(starts, ends) -> setup_means
  setup$concat.pos <- apply(setup_means, 1, function(x) mean(x))
  plyr::ddply(cov.df, ~Id + run, slide, setup) -> cov.slid.df
  
  x.labels <- plyr::ddply(cov.slid.df,~chr,plyr::summarize,concat.pos=concat.pos[which(abs(concat.pos-mean(concat.pos))==(min(abs(concat.pos-mean(concat.pos)))))])
  x.labels <- plyr::ddply(x.labels, ~chr, function(x) return(x[1,]))
  
  cov.plot <- ggplot(cov.slid.df, mapping = aes(x = as.factor(concat.pos), y = mean)) + geom_boxplot(fill="white")
  cov.plot <- cov.plot + ggtitle(title) + ylab("Read depth") + scale_x_discrete(labels = x.labels$chr, breaks = x.labels$concat.pos) + xlab("Concatenated Genome Position")
  cov.plot <- cov.plot + theme(axis.title.y = element_text(vjust=1.2))
  cov.plot <- cov.plot + theme(legend.position = "none") + theme_classic()
  cov.plot <- cov.plot + theme(text = element_text(size = 35), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))
  return(cov.plot)
}

# ======================== Plot ===================

# Combine and reorder so that PB2 is first (PB1 was first in the reference file)
all.cov <- rbind(vic.cov, yam.cov)
pb2.cov <- filter(all.cov, chr == "PB2")
pb2.pos.max <- max(pb2.cov$chr.pos)
pb1.cov <- filter(all.cov, chr == "PB1")
pb1.concat.max <- max(pb1.cov$concat.pos)
pb2.cov <- mutate(pb2.cov, concat.pos = concat.pos - pb1.concat.max)
pb1.cov <- mutate(pb1.cov, concat.pos = concat.pos + pb2.pos.max)
all.cov.rest <- filter(all.cov, chr != "PB2" & chr != "PB1")
all.cov.reorder <- rbind(pb2.cov, pb1.cov, all.cov.rest)

cov_plot(all.cov.reorder, title = "") -> coverage.plot.all

# ========================= Figure 1B ==============

ggsave(filename = "results/plots/coverage_plot.pdf", plot = all.cov.reorder, device = "pdf", width = 6, height = 4)

