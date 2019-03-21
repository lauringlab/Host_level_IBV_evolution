
### Project: IBV intra-host
### Purpose: Get average coverage per sample.

# ======================== Import packages ==========================

library(magrittr)
library(tidyverse)

# ======================== Get average coverage ==========================

# Here we calculate the average coverage for each sequenced sample. As an aside this coverage does not include reads that had gaps at the position in question.

cov <- read_csv("data/processed/all.coverage.csv")

cov %>% group_by(Id) %>% summarize(coverage = mean(coverage)) %>% select(Id, coverage) -> cov_by_sample

write.csv(cov_by_sample, "data/processed/average_coverages.csv")


