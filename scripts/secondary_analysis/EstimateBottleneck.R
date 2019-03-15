

### Project: IBV intra-host
### Purpose: Estimate the bottleneck using the presence-absence model and the beta-binomial model.

# ================== Read in data, import packages =========================

library(tidyverse)
library(wesanderson)
library(bbmle)
library(HIVEr)
detach(package:plyr)

palette <- wes_palette("Darjeeling1")
setwd("/Users/avalesano/Documents/MSTP/LauringLab/Host_level_IBV_evolution/scripts/")

meta_long <- read_csv("../data/metadata/flu_b_2010_2017_v4LONG_withSeqInfo_gc.csv")
possible_pairs_dist <- read_csv("../data/processed/possible.pairs.dist.csv")
transmission_pairs <- read_csv("../data/processed/transmission_pairs.csv")

transmission_freqs <- read_csv("../data/processed/trans_freq.csv")
transmission_freqs_poly_in_donor <- read_csv("../data/processed/transmission_pairs_freq.poly_donor.csv")
transmission_freqs_no_cut <- read_csv("../data/processed/no_cut_trans_freq.csv")
transmission_freqs_no_cut_poly_in_donor <- read_csv("../data/processed/no_cut_transmission_pairs_freq.poly_donor.csv")


# ================== Quick and Dirty: Use all of JT's support functions in HIVEr/transmission_models.R to do PA model for qual SNV. =========================

### Select which dataset to use: 2% frequency cutoff or no cutoff
trans_freq.comp <- transmission_freqs_poly_in_donor
trans_freq.comp <- transmission_freqs_no_cut_poly_in_donor

### Run the maximum likelihood optimization. // ALV: why isn't it grouping by pair_id?
pa_total_fit <- trans_fit(trans_freq.comp, Nb_max = 100, model = "PA", threshold = NULL, acc = NULL, pair_id)
trans_freq.comp %>% group_by(pair_id) %>% summarize(donor_mutants = length(which(freq1 > 0 & freq1 < 0.5))) %>% mutate(weight_factor_kk = max(donor_mutants)/donor_mutants, weight_factor = 1) -> counts
zdpois_fit <- dist_prob_wrapper(ddist = "dzpois", params = "lambda")
dzpois_model_fit <- bbmle::mle2(minuslogl = zdpois_fit, start = list(lambda = 1), data = list(data = pa_total_fit, weight = counts))

### Summarize the model fit
fit <- profile(dzpois_model_fit)
plot(fit)
con_int <- confint(dzpois_model_fit)
AIC(dzpois_model_fit)
summary(dzpois_model_fit)

### Plot the fit with the data by window. // May need to change rowwise...
w <- 0.05
step <- 0.025
windows <- tibble(s = seq(0.02, 1-w, by = step), end = s + w)

windows %>% rowwise() %>%
  mutate(iSNV = nrow(filter(trans_freq.comp, freq1 >= s, freq1 < end)),
         transmitted = nrow(filter(trans_freq.comp, freq1 >= s, freq1 < end, found == TRUE)),
         freq = mean(trans_freq.comp$freq1[which(trans_freq.comp$freq1 >= s & trans_freq.comp$freq1 < end)]),
         prob = transmitted/iSNV,
         error_bottom = qbinom(c(0.025), iSNV, prob)/iSNV,
         error_top = qbinom(c(0.975), iSNV, prob)/iSNV,
         many = iSNV > 5) -> out

Pt_PA <- function(x, l, max_Nb) # Math given in McCrone et al. 2018. Probability of transmission in the presence-absence model given a bottleneck size Nb.
{
  s <- 0
  for(i in 1:max_Nb)
  {
    c <- ((1 - x)^i) * (l^i)/((exp(l) - 1)*factorial(i))
    s <- s + c
  }
  return(1 - s)
}

model <- tibble(s = seq(0, 1, 0.01))
model <- mutate(model, prob = Pt_PA(s, dzpois_model_fit@coef, 100),
                lower = Pt_PA(s, con_int[1], 100),
                upper = Pt_PA(s, con_int[2], 100),
                diff = upper - lower)

window_data.plot <- ggplot() + geom_point(data = out, aes(x = freq, y = prob, alpha = many)) +
  geom_errorbar(data = out, aes(x = freq, ymin = error_bottom, ymax = error_top, alpha = many)) +
  geom_line(data = model, aes(x = s, y = prob), color = palette[1]) +
  geom_ribbon(data = model, aes(x = s, ymin = lower, ymax = upper), alpha = 0.5, fill = palette[5]) +
  scale_alpha_manual(values = c(0, 1)) + theme(legend.position = 'none') +
  xlab("Frequency in Donor") + ylab("Probability Transmitted") +
  geom_point(data = trans_freq.comp,
             aes(x = freq1, y = as.numeric(found) + (as.numeric(found) - 0.5)/10), alpha = 0.5)+
  scale_y_continuous(breaks = seq(0, 1, 0.25))
window_data.plot

# Summary: for the 2% cutoff data, we get an estimate of -30 with no confidence interval. 
# For the no-cutoff data, we get an estimate of 87 with a confidence interval of 81 - 92.
# This is a very fragile dataset. Let's see what happens with the beta-binomial model.

# ===================== Do the same thing for the beta-binomial model ===========================

### Select which dataset to use: 2% frequency cutoff or no cutoff
trans_freq.comp <- transmission_freqs_poly_in_donor
#trans_freq.comp <- transmission_freqs_no_cut_poly_in_donor

### Need to add genome copy info
trans_freq.comp <- mutate(trans_freq.comp, gc_ul1 = meta_long$genome_copy_per_ul[match(ALV_ID_1, meta_long$ALV_ID)], gc_ul2 = meta_long$genome_copy_per_ul[match(ALV_ID_2, meta_long$ALV_ID)])

### Do the optimization
# // ALV: error here. object Nb not found.
beta_total_fit <- trans_fit(filter(trans_freq.comp, freq1 < 0.5), Nb_max = 100, model = "BetaBin", threshold = 0.02, acc = accuracy_stringent, pair_id)

zdpois_fit <- dist_prob_wrapper(ddist = "dzpois", params = "lambda")

trans_freq.comp %>% group_by(pair_id) %>%
  summarize(donor_mutants = length(which(freq1 > 0 & freq1 < 0.5))) %>%
  mutate(weight_factor_kk = max(donor_mutants)/donor_mutants,
         weight_factor = 1) -> counts

dzpois_model_fit_bb <- bbmle::mle2(minuslogl = zdpois_fit, start = list(lambda = 1),
                                 data = list(data = beta_total_fit,
                                             weight = counts))

# Characterize the model fit
conf_int_BB <- bbmle::confint(dzpois_model_fit_bb)
summary(dzpois_model_fit_bb)
AIC(dzpois_model_fit_bb)

### Plot // ALV: need to modify

model_betaBin<-tibble(s = seq(0,1,0.01))
model_betaBin<- model_betaBin %>% rowwise() %>%
  mutate(prob = Pt_BetaBin(s,dzpois_model_fit_bb@coef,100),
         lower = Pt_BetaBin(s,conf_int_BB[1],100),
         upper = Pt_BetaBin(s,conf_int_BB[2],100))

window_data_bb.p<-ggplot()+geom_point(data = out,aes(x=freq,y=prob,alpha=many))+
  geom_errorbar(data=out,aes(x=freq,ymin=error_bottom,ymax=error_top,alpha=many))+
  geom_line(data=model_betaBin,aes(x=s,y=prob),color=cbPalette[5])+
  geom_ribbon(data=model_betaBin,aes(x=s,ymin=lower,ymax=upper),alpha=0.5,fill=cbPalette[5])+
  scale_alpha_manual(values=c(0,1))+theme(legend.position = 'none')+
  xlab("Frequency in Donor")+ylab("Probability of transmission")+
  geom_point(data=trans_freq.comp,
             aes(x=freq1,y=as.numeric(found)+(as.numeric(found)-0.5)/10),alpha=0.5)+
  scale_y_continuous(breaks = seq(0,1,0.25))

