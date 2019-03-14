

### Project: IBV intra-host
### Purpose: Estimate the bottleneck using the presence-absence model and the beta-binomial model.

# ================== Read in data, import packages =========================

library(tidyverse)
library(wesanderson)
library(bbmle)
library(HIVEr)

palette <- wes_palette("Darjeeling1")
setwd("/Users/avalesano/Documents/MSTP/LauringLab/Host_level_IBV_evolution/scripts/")

possible_pairs_dist <- read_csv("../data/processed/possible.pairs.dist.csv")
transmission_pairs <- read_csv("../data/processed/transmission_pairs.csv")

transmission_freqs <- read_csv("../data/processed/trans_freq.csv")
transmission_freqs_poly_in_donor <- read_csv("../data/processed/transmission_pairs_freq.poly_donor.csv")
transmission_freqs_no_cut <- read_csv("../data/processed/no_cut_trans_freq.csv")
transmission_freqs_no_cut_poly_in_donor <- read_csv("../data/processed/no_cut_transmission_pairs_freq.poly_donor.csv")

# ================== JT's bottleneck helper functions in HIVEr/transmission_models.R ==================

is_wholenumber <-function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}

dzpois<-function(x,lambda){ # the pmf of the zero truncated poisson distribution
  if(x<=0 | is_wholenumber(x)==F){
    return(0)
  }else{
    (exp(-1*lambda)*lambda^x)/(factorial(x)*(1-exp(-1*lambda)))
    
  }
}
dzpois<-Vectorize(dzpois,vectorize.args = c("x"))

p_all<-function(p,n){ # probability all success - only finding the variant at frequency p in n draws
  p^n
}
p_all<-Vectorize(p_all,vectorize.args="n")

math_fit=function(data,Nb_max,model,threshold,acc){
  # this  calculation is for each position in the genome.
  stopifnot(length(unique(data$chr))==1, length(unique(data$pos))==1)
  # and each model.
  if(model=="PA"){#| model=="PA-straight"){
    #  In this fit we take the minority frequency to be  correct
    # and set the major frequency to 1-minority. This accounts for
    # the fact that frequencies are related by the expression :
    # minor allele + major allele + errror =1.
    #Here we make the major allele frequency = major allele + error.
    #The error is always small. if it exceeds 1% then we through an error here.
    
    if(1-sum(data$freq1)>0.01){
      warning("The sum of the frequencies is less than 99% make sure the assumption major freq = 1-minor freq is valid")
    }
    data$freq1[data$freq1==max(data$freq1)]<-1-min(data$freq1)
    
    
    found<-data[data$found==T,] # only alleles that were transmitted
    Nb<-1:Nb_max # Here are the bottlenecks
    
    if(nrow(found)==0 | nrow(found)>2){
      stop(paste0("No variant transmitted for this site or",
                  "there are more than 2 variants here"))
    }else if(nrow(found)==1){ # one variant found here. All successes
      prob<-p_all(p=found$freq1,n=Nb)
      # this is a vector of probabilities for each
      # n prob[i]= the probability of only getting that
      # variant in Nb[i] (i.e. draws)
    }else if(nrow(found)==2){
      # if at least on of each allele was transmitted
      
      first_var<-p_all(p=found$freq1[1],n=Nb) # all this one
      second_var<-p_all(p=found$freq1[2],n=Nb) # all the other one
      one_each<-1-(first_var+second_var)
      # at least one of each -
      # This is a vector as above since R adds and subtracts the
      # elements of the vectors as expected
      prob<-one_each
    }
    
  }else if(model=="BetaBin"){# | model =="BetaBin-straight"){
    if(nrow(data[data$freq1>0.5,])>0){
      data<-data[df$freq1<0.5,]
      warning("The beta binomials model only uses minor alleles. Subsetting the data now.")
    }
    # 2 is the recipient
    v_r = data$freq2
    v_d = data$freq1
    gc_ul = data$gc_ul2
    threshold = 0.02
    
    prob = L.Nb.beta(v_r,v_d,Nb,gc_ul,threshold,acc)
  }
  return(tibble(Nb=Nb,prob=prob))
  
}

trans_fit<-function(data,Nb_max,model,threshold,acc,...){
  
  group <- rlang::quos(...,Nb)
  original_group <- rlang::quos(...)
  probs<-data %>% dplyr::group_by(chr,pos,pair_id) %>%
    dplyr::do(math_fit(.data,Nb_max,model,threshold,acc))
  # For each genomic position in question
  
  #control for number of the mutations in each donor.
  counts <- data %>% dplyr::group_by(!!!original_group) %>%
    dplyr::summarise(donor_mutants = length(which(freq1>0 & freq1<0.5)))
  max_mut<-max(counts$donor_mutants)
  LL.df<-probs %>% dplyr::left_join(counts,by = "pair_id") %>%
    dplyr::group_by(!!!group) %>%
    dplyr::summarize(LL=sum(log(prob)),
                     weighted_LL = (max_mut/unique(donor_mutants)) * LL)
  # Get the  log likelihood of the each lambda for this pair - we can sum accros pairs later.
  return(LL.df)
}

dist_prob <- function(data,weight,ddist,...){
  params <- rlang::enquos(...)
  ddist <-rlang::enexpr(ddist)
  
  Nb = unique(data$Nb) # probability of data
  data<- dplyr::mutate(data,prob_D = exp(LL))
  # I don't know the details of quoting ect. But this works.
  # prob of Nb
  prob_Nb<-rlang::eval_tidy(quo(Nb %>% purrr::map_dbl(.f = !!ddist,!!!params)))
  data<-data %>%
    dplyr::left_join(dplyr::tibble(Nb=Nb,prob_Nb = prob_Nb),by="Nb") %>%
    dplyr::mutate(prob_D_and_Nb = prob_D*prob_Nb)
  
  l_by_pair<- data %>% dplyr::group_by(pair_id) %>%
    dplyr::summarize(prob_D_given_dist = sum(prob_D_and_Nb),
                     LL_D_given_dist = log(prob_D_given_dist)) %>%
    dplyr::rowwise()%>%
    dplyr::mutate(weighted_total_LL  = weight$weight_factor[weight$pair_id==pair_id] * LL_D_given_dist)
  return(-sum(l_by_pair$weighted_total_LL))
  
}

dist_prob_wrapper<-function(ddist,params){
  f_string<- paste0("function(data,weight,",params,"){dist_prob(data,weight,",ddist,",", params,")}")
  eval(parse(text = f_string))
}

model_summary<-function(data){
  Nb<-data$Nb[which(data$LL==max(data$LL))] # Get the  max lambda
  good_range<-subset(data,LL> (max(LL)-1.92)) # get the bottlenecks that fall in this region the 95% confidence intereval
  lower<-good_range$Nb[1]
  upper <- good_range$Nb[nrow(good_range)]
  return(tibble(Nb=Nb,lower_95=lower,
                upper_95=upper))
}

pa_sim<-function(data,lambda){
  pair<-unique(data$pair_id)
  if(length(pair)>1){
    warning(paste0("Running on ",length(pair)," pairs. All will have the same bottleneck."))
  }
  Nb<-rzpois(1,lambda)
  out<-data %>% dplyr::group_by(chr,pos,pair_id) %>%
    dplyr::do(pa_sim_helper(.,Nb))
  return(out)
}

pa_sim_helper<-function(data,Nb){
  # data refers to the 2 mutations at this position, bottlenecks-
  #  a data frame with the bottle_neck,
  #model - the name of the column in the bottlenecks that
  # contains the bottleneck size we want to use.
  if(nrow(data)!=2){
    stop(print("There are not 2 mutations at this point."))
  }
  pair<-unique(data$pair_id)
  
  if(length(pair)!=1){
    stop(print("There should only be one pair id at this point."))
  }
  freq_success<-1-min(data$freq1)
  # This is the major allele frequency - calculated the same as in the model
  
  success<-rbinom(1,Nb,freq_success)
  # the number of success in n trials.
  #  How many of the major alleles in the draw
  #
  data$found<-F
  if(success==Nb){ # only found major variants
    data$found[data$freq1==max(data$freq1)]=T
  } else if(success==0){ # only the minor
    data$found[data$freq1==min(data$freq1)]=T
  } else if(success>0 & success<Nb){   # both were found
    data$found=T
  } else{
    stop("Error!")
  }
  return(data)
}

rzpois<-function(n,lambda){
  #from http://giocc.com/zero_truncated_poisson_sampling_algorithm.html
  out = vector(mode="double",length=n)
  for(i in 1:n){
    k = 1
    t = exp(-lambda) / (1 - exp(-lambda)) * lambda
    s = t
    u = runif(1)
    while(s < u){
      k = k+1
      t = t*(lambda / k)
      s = s+t
    }
    out[i] <- k
  }
  return(out)
}

simulations<-function(data,runs,lambda,FUN,threshold=NULL,acc=NULL,...){
  # iSNV data, how many iterations,
  # what is lambda value and what function is used to simulate the data.
  
  group_by<-rlang::quos(...)
  rows_needed<-nrow(data)
  model.df<-tibble(freq1=rep(NA,rows_needed*runs),
                   trial=rep(1:runs,each=rows_needed),
                   prob=NA)
  pairs<-length(unique(data$pair_id))
  for (i in 1:runs){
    
    if(is.null(threshold)){
      trial.df<-data %>% dplyr::group_by(!!!group_by)%>%
        dplyr::do(FUN(.,lambda))
    }else{
      trial.df<-data %>% dplyr::group_by(!!!group_by)%>%
        dplyr::do(FUN(.,lambda,threshold,acc))
    }
    logit<-glm(formula =found~freq1,family=binomial(logit),data=trial.df) # Fit a logit model to the data
    trial.df$prob<-logit$fitted.values
    trial.df$trial<-i
    ind<-which(model.df$trial==i)
    model.df$prob[ind]<-trial.df$prob# add to the final output
    model.df$freq1[ind]<-trial.df$freq1
  }
  return(model.df)
}


# ================== Quick and Dirty: Use all of JT's support functions in HIVEr/transmission_models.R to do PA model for qual SNV. =========================

trans_freq.comp <- transmission_freqs_poly_in_donor
trans_freq.comp <- transmission_freqs_no_cut_poly_in_donor

pa_total_fit <- trans_fit(trans_freq.comp, Nb_max = 100, model = "PA", threshold = NULL, acc = NULL, pair_id)

counts<-trans_freq.comp %>% group_by(pair_id) %>%
  summarize(donor_mutants = length(which(freq1>0 & freq1<0.5))) %>%
  mutate(weight_factor_kk = max(donor_mutants)/donor_mutants,
         weight_factor = 1)

zdpois_fit<-dist_prob_wrapper(ddist = "dzpois",params = "lambda")

dzpois_model_fit<-bbmle::mle2(minuslogl = zdpois_fit,start = list(lambda = 1),
                              data = list(data = pa_total_fit,
                                          weight = counts))

mean_zpois<-function(l) l/(1-exp(-1*l))
con_int<-confint(dzpois_model_fit)

AIC(dzpois_model_fit)
summary(dzpois_model_fit)

w<-0.05
step = 0.025
windows<-tibble(s=seq(0.02,1-w,by=step),end=s+w)

out <- windows %>% rowwise() %>%
  mutate(iSNV = nrow(filter(trans_freq.comp,freq1>=s,freq1<end)),
         transmitted = nrow(filter(trans_freq.comp,freq1>=s,freq1<end,found==T)),
         freq = mean(trans_freq.comp$freq1[which(trans_freq.comp$freq1>=s & trans_freq.comp$freq1<end)]),
         prob = transmitted/iSNV,
         error_bottom = qbinom(c(0.025),iSNV,prob)/iSNV,
         error_top = qbinom(c(0.975),iSNV,prob)/iSNV,
         many = iSNV>5)

Pt_PA<-function(x,l,max_Nb){
  s<-0
  for(i in 1:max_Nb){
    c<- ((1-x)^i) *  (l^i)/((exp(l)-1)*factorial(i) )
    s<-s+c
  }
  return(1-s)
}

model<-tibble(s = seq(0,1,0.01))
model<- mutate(model,prob = Pt_PA(s,dzpois_model_fit@coef,100),
               lower = Pt_PA(s,con_int[1],100),
               upper = Pt_PA(s,con_int[2],100))

window_data.p<-ggplot()+geom_point(data = out,aes(x=freq,y=prob,alpha=many))+
  geom_errorbar(data=out,aes(x=freq,ymin=error_bottom,ymax=error_top,alpha=many))+
  geom_line(data=model,aes(x=s,y=prob),color=palette[5])+
  geom_ribbon(data=model,aes(x=s,ymin=lower,ymax=upper),alpha=0.5,fill=palette[5])+
  scale_alpha_manual(values=c(0,1))+theme(legend.position = 'none')+
  xlab("Frequency in Donor")+ylab("Probability of transmission")+
  geom_point(data=trans_freq.comp,
             aes(x=freq1,y=as.numeric(found)+(as.numeric(found)-0.5)/10),alpha=0.5)+
  scale_y_continuous(breaks = seq(0,1,0.25))

### Ok, this has problems. Weird-looking curve. I will need to go through each and every step, or just try to do this analysis from scratch.
