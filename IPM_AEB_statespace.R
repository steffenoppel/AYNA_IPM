#### LOAD LIBRARIES ####
library(nimble)
library(jagsUI)
library(coda)
library(doParallel) 
library(tidyverse)
library(tidybayes)
library(strex)
library(beepr)
library(postpack)

#### LOAD DATA ####
load("IPM_AEB_dat_stateSpace.RData") # TODO fix
const <- c(const, maxAge = const$n.occasions)

#### LOAD INITS #####
load("IPM_AEB_inits_stateSpace.Rdata") # TODO fix

# TODO
# run this in both jags and nimble for a sanity check
# limited variability in goodyear - do we want to change threshold?
# how do we want to deal with max age

#### MODEL CODE ####
code <- nimbleCode({
  
  #-------------------------------------------------  
  # 1. PRIORS FOR ALL DATA SETS
  #-------------------------------------------------
  
  # -------------------------------------------------        
  # 1.1. Priors and constraints FOR FECUNDITY
  # -------------------------------------------------
  
  for (t in 1:n.years.fec){  
    ann.fec[t] ~ dbeta(10,10) # AEB change, more variation in data
  } # t
  
  # -------------------------------------------------        
  # 1.2. Priors and constraints FOR POPULATION COUNTS
  # -------------------------------------------------
  
  # option 2 - index by s and goodyear/badyear
  # for (s in 1:n.sites.count){
  #     sigma.obs[s] ~ dexp(0.1)
  # } # s
  
  # -------------------------------------------------        
  # 1.3. Priors and constraints FOR SURVIVAL
  # -------------------------------------------------
  
  mean.p.in ~ dbeta(1,1)
  mean.p.propensity ~ dbeta(1,1)
  mean.p.det ~ dbeta(1,1)

  ### RECAPTURE PROBABILITY
  #mean.p.ad[1] ~ dbeta(8.3, 25)	# low detection
  mean.p.ad ~ dbeta(1, 1)	# low detection

  mu.p.juv ~ dnorm(-4, sd = 0.25) 
  mu.p.ad <- log(mean.p.ad / (1-mean.p.ad)) 
  
  # Prior for shape of increase in juvenile recapture probability with age
  #agebeta ~ dnorm(2, sd = 0.25) # AEB change
  agebeta ~ dnorm(2, sd = 1) # AEB change
  
  sigma.p ~ dexp(10) # AEB change, reduce variability 
  
  ## RANDOM TIME EFFECT ON RESIGHTING PROBABILITY...
  for (t in 1:(n.occasions-1)){
    # ...FOR JUVENILES
    for (j in 1:t){ 
      p.juv[t,j] <- 0
    }
    for (j in (t+1):(n.occasions)){
      logit(p.juv[t,j])  <- mu.p.juv + agebeta*sqrt(j - t) + eps.p[j-1]
      
      mean.p.recruit[1, t, j] <- p.juv[t, j]
      mean.p.recruit[2, t, j] <- p.ad[t]
    }
    # ...FOR ADULTS
    logit(p.ad[t])  <- mu.p.ad + eps.p[t]  
    eps.p[t] ~ dnorm(0, sd = sigma.p)
  }
  
  ### SURVIVAL PROBABILITY
  #mean.phi.juv ~ dbeta(10.75, 3.5)   
  mean.phi.juv ~ dbeta(1, 1)            
  #mean.phi.ad ~ dbeta(50,2)  
  mean.phi.ad ~ dbeta(1,1)              
  
  mu.juv <- log(mean.phi.juv / (1-mean.phi.juv)) # Logit transformation
  mu.ad <- log(mean.phi.ad / (1-mean.phi.ad)) # Logit transformation
  
  sigma.phi ~ dexp(10) # AEB changed, don't want this huge
  
  ## RANDOM TIME EFFECT ON SURVIVAL
  for (j in 1:(n.occasions-1)){
    logit(phi.juv[j]) <- mu.juv + eps.phi[j]*juv.poss[j] # TODO - does this even matter?
    # there is already a lot of variability in phi.juv, don't think eps.phi matters as much for this age class
    # more important for adults
    logit(phi.ad[j]) <- mu.ad + eps.phi[j] 
    eps.phi[j] ~ dnorm(0, sd = sigma.phi) 
    
    mean.phi[1, j] <- phi.juv[j]
    mean.phi[2, j] <- phi.ad[j]
  }
  
  #-------------------------------------------------  
  # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
  #-------------------------------------------------
  
  # -------------------------------------------------        
  # 2.1. System process: female based matrix model
  # -------------------------------------------------
  
  ### RECRUIT PROBABILITY ###
  for (age in 1:maxAge) {
    logit(p.juv.recruit.f[age]) <- mu.p.juv + agebeta*sqrt(age)
  }
  
  #IM[1,1,1] ~ T(dnorm(263/2,sd = 20), 0, Inf)
  IM[1,1,1] ~ dpois(263/2)
  IM[1,1,2] <- 0
  #IM[1,1,3] <- round(IM[1,1,1]) - IM[1,1,2]
  IM[1,1,3] <- IM[1,1,1] - IM[1,1,2]
  
  #IM[1,2,1] ~ T(dnorm(275/2,sd = 20), 0, Inf)
  IM[1,2,1] ~ dpois(275/2)
  IM[1,2,2] ~ dbin(p.juv.recruit.f[2], round(IM[1,2,1]))
  #IM[1,2,3] <- round(IM[1,2,1]) - IM[1,2,2]
  IM[1,2,3] <- IM[1,2,1] - IM[1,2,2]
  
  #IM[1,3,1] ~ T(dnorm(264/2,sd = 20), 0, Inf)
  IM[1,3,1] ~ dpois(264/2)
  IM[1,3,2] ~ dbin(p.juv.recruit.f[3], round(IM[1,3,1]))
  #IM[1,3,3] <- round(IM[1,3,1]) - IM[1,3,2]
  IM[1,3,3] <- IM[1,3,1] - IM[1,3,2]
  
  #IM[1,4,1] ~ T(dnorm(177/2,sd = 20), 0, Inf) 
  IM[1,4,1] ~ dpois(177/2)
  IM[1,4,2] ~ dbin(p.juv.recruit.f[4], round(IM[1,4,1]))
  #IM[1,4,3] <- round(IM[1,4,1]) - IM[1,4,2]
  IM[1,4,3] <- IM[1,4,1] - IM[1,4,2]
  
  #IM[1,5,1] ~ T(dnorm(290/2,sd = 20), 0, Inf)
  IM[1,5,1] ~ dpois(290/2)
  IM[1,5,2] ~ dbin(p.juv.recruit.f[5], round(IM[1,5,1]))
  #IM[1,5,3] <- round(IM[1,5,1]) - IM[1,5,2]
  IM[1,5,3] <- IM[1,5,1] - IM[1,5,2]
  
  #IM[1,6,1] ~ T(dnorm(90/2,sd = 20), 0, Inf) 
  IM[1,6,1] ~ dpois(90/2)
  IM[1,6,2] ~ dbin(p.juv.recruit.f[6], round(IM[1,6,1]))
  #IM[1,6,3] <- round(IM[1,6,1]) - IM[1,6,2]
  IM[1,6,3] <- IM[1,6,1] - IM[1,6,2]
  
  #IM[1,7,1] ~ T(dnorm(158/2,sd = 20), 0, Inf)
  IM[1,7,1] ~ dpois(158/2)
  IM[1,7,2] ~ dbin(p.juv.recruit.f[7], round(IM[1,7,1]))
  #IM[1,7,3] <- round(IM[1,7,1]) - IM[1,7,2]
  IM[1,7,3] <- IM[1,7,1] - IM[1,7,2]
  
  for(age in 8:maxAge) {
    IM[1,age,1] ~ dbin(mean.phi.ad, IM[1,age-1,3])
    IM[1,age,2] ~ dbin(p.juv.recruit.f[age], IM[1,age,1])
    IM[1,age,3] <- IM[1,age,1] - IM[1,age,2]
  }
  
  N.recruits[1] <- sum(IM[1,1:maxAge,2]) 
  N.ad.surv[1] <- 0
  N.breed.ready[1] <- 0
  #Ntot.breed[1] ~ T(dnorm(640/2,sd = 20), 0, Inf)  
  Ntot.breed[1] ~ dpois(640/2)  
  #N.atsea[1] ~ T(dnorm(224/2,sd = 20), 0, Inf)    
  N.atsea[1] ~ dpois(224/2)    
  nestlings[1] <- 0
  #JUV[1] ~ T(dnorm(232/2, sd = 20), 0, Inf)         
  JUV[1] ~ dpois(232/2)         
  
  #Ntot[1]<-sum(IM[1,1:maxAge,3]) + round(Ntot.breed[1])+round(N.atsea[1]) 
  Ntot[1]<-sum(IM[1,1:maxAge,3]) + Ntot.breed[1] + N.atsea[1] 
  
  ### FOR EVERY SUBSEQUENT YEAR POPULATION PROCESS
  
  for (tt in 2:n.years.fec){
    for (age in 1:maxAge) {
      logit(p.juv.recruit[age,tt]) <- mu.p.juv + eps.p[tt+offset-1] + agebeta * sqrt(age)
    }
    
    ## IMMATURE MATRIX WITH 3 columns:
    # 1: survivors from previous year
    # 2: recruits in current year
    # 3: unrecruited in current year (available for recruitment next year)
    
    #IM[tt,1,1] ~ dbin(phi.juv[tt+offset-1], round(JUV[tt-1]))                                
    IM[tt,1,1] ~ dbin(phi.juv[tt+offset-1], JUV[tt-1])                                
    IM[tt,1,2] <- 0 
    IM[tt,1,3] <- IM[tt,1,1] - IM[tt,1,2]
    for(age in 2:maxAge) {
      IM[tt,age,1] ~ dbin(phi.ad[tt+offset-1], IM[tt-1,age-1,3])
      IM[tt,age,2] ~ dbin(p.juv.recruit[age,tt], IM[tt,age,1])
      IM[tt,age,3] <- IM[tt,age,1] - IM[tt,age,2]
    }
    
    N.recruits[tt] <- sum(IM[tt,1:maxAge,2])  ### number of this years recruiters
    #N.ad.surv[tt] ~ dbin(phi.ad[tt+offset-1], round(Ntot.breed[tt-1])+round(N.atsea[tt-1]))  ### previous year's adults that survive
    N.ad.surv[tt] ~ dbin(phi.ad[tt+offset-1], Ntot.breed[tt-1]+N.atsea[tt-1])  ### previous year's adults that survive
    N.breed.ready[tt] ~ dbin(p.ad[tt+offset-1], N.ad.surv[tt]) ### number of available breeders is proportion of survivors that returns
    Ntot.breed[tt]<- N.breed.ready[tt]+N.recruits[tt]  ### number of counted breeders is sum of old breeders returning and first recruits
    N.atsea[tt] <- N.ad.surv[tt]-N.breed.ready[tt] ### potential breeders that remain at sea
    nestlings[tt] <- round(ann.fec[tt] * 0.5 * Ntot.breed[tt]) ### nestlings produced, half are female                                                 
    JUV[tt] ~ dpois(nestlings[tt]) ### juveniles produced
    
    ### THE TOTAL POPULATION ###
    Ntot[tt]<-sum(IM[tt,1:maxAge,3]) + Ntot.breed[tt]+N.atsea[tt]  ## total population size is all the immatures plus adult breeders and adults at sea
    
  } # tt
  
  # -------------------------------------------------
  # 2.2. Observation process for population counts: state-space model of annual counts
  # -------------------------------------------------
  
  # for (s in 1:n.sites.count){
  #   for (t in 1:n.years.fec){
  #     y.count[t,s] ~ T(dnorm(Ntot.breed[t]*prop.sites[t,s]*2, sd = sigma.obs[s]), 0, Inf)
  #     y.count[t,s] ~ dpois(Ntot.breed[t]*prop.sites[t,s]*2)
  #   }	# t
  # }	# s
  
  # truncated normal looks bad. let's just try the poisson, which is going to be more stable if we sum over sites
  for (t in 1:n.years.fec){
    y.count[t] ~ dpois(Ntot.breed[t]*2)
  }	# t
  
  # -------------------------------------------------
  # 2.3. Likelihood for fecundity: Logistic regression from the number of surveyed broods
  # -------------------------------------------------
  
  for (t in 1:(n.years.fec)){ 
    J[t] ~ dbin(ann.fec[t], R[t])
  } 
  
  # -------------------------------------------------
  # 2.4. Likelihood for adult and juvenile survival from CMR
  # -------------------------------------------------
  #
  for (i in 1:n.inds) {
    for (t in (first[i]+1):n.occasions) {
      # phi[i,t] related to individual age and year per previous model
      phi[i,t] <- mean.phi[ageCat[i,t-1], t-1]
      
      p.recruit[i,t] <- mean.p.recruit[ageCat[i,t-1], t-1, t]
      
      p.in[i,t] <- mean.p.in
      p.propensity[i,t] <- mean.p.propensity
      
      trans.mat[i,t,1,1] <- phi[i,t]*(1-p.recruit[i,t])
      trans.mat[i,t,1,2] <- 0
      trans.mat[i,t,1,3] <- phi[i,t]*p.recruit[i,t]*p.in[i,t]
      trans.mat[i,t,1,4] <- phi[i,t]*p.recruit[i,t]*(1-p.in[i,t])
      trans.mat[i,t,1,5] <- 0
      trans.mat[i,t,1,6] <- 0
      trans.mat[i,t,1,7] <- 1-phi[i,t]

      trans.mat[i,t,2,1] <- 0
      trans.mat[i,t,2,2] <- phi[i,t]*(1-p.propensity[i,t])
      trans.mat[i,t,2,3] <- 0
      trans.mat[i,t,2,4] <- 0
      trans.mat[i,t,2,5] <- phi[i,t]*p.propensity[i,t]*p.in[i,t]
      trans.mat[i,t,2,6] <- phi[i,t]*p.propensity[i,t]*(1-p.in[i,t])
      trans.mat[i,t,2,7] <- 1-phi[i,t]
      
      trans.mat[i,t,3,1] <- 0
      trans.mat[i,t,3,2] <- phi[i,t]*(1-p.propensity[i,t])
      trans.mat[i,t,3,3] <- 0
      trans.mat[i,t,3,4] <- 0
      trans.mat[i,t,3,5] <- phi[i,t]*p.propensity[i,t]*p.in[i,t]
      trans.mat[i,t,3,6] <- phi[i,t]*p.propensity[i,t]*(1-p.in[i,t])
      trans.mat[i,t,3,7] <- 1-phi[i,t]
      
      trans.mat[i,t,4,1] <- 0
      trans.mat[i,t,4,2] <- phi[i,t]*(1-p.propensity[i,t])
      trans.mat[i,t,4,3] <- 0
      trans.mat[i,t,4,4] <- 0
      trans.mat[i,t,4,5] <- phi[i,t]*p.propensity[i,t]*p.in[i,t]
      trans.mat[i,t,4,6] <- phi[i,t]*p.propensity[i,t]*(1-p.in[i,t])
      trans.mat[i,t,4,7] <- 1-phi[i,t]
      
      trans.mat[i,t,5,1] <- 0
      trans.mat[i,t,5,2] <- phi[i,t]*(1-p.propensity[i,t])
      trans.mat[i,t,5,3] <- 0
      trans.mat[i,t,5,4] <- 0
      trans.mat[i,t,5,5] <- phi[i,t]*p.propensity[i,t]*p.in[i,t]
      trans.mat[i,t,5,6] <- phi[i,t]*p.propensity[i,t]*(1-p.in[i,t])
      trans.mat[i,t,5,7] <- 1-phi[i,t]
      
      trans.mat[i,t,6,1] <- 0
      trans.mat[i,t,6,2] <- phi[i,t]*(1-p.propensity[i,t])
      trans.mat[i,t,6,3] <- 0
      trans.mat[i,t,6,4] <- 0
      trans.mat[i,t,6,5] <- phi[i,t]*p.propensity[i,t]*p.in[i,t]
      trans.mat[i,t,6,6] <- phi[i,t]*p.propensity[i,t]*(1-p.in[i,t])
      trans.mat[i,t,6,7] <- 1-phi[i,t]
      
      trans.mat[i,t,7,1] <- 0
      trans.mat[i,t,7,2] <- 0
      trans.mat[i,t,7,3] <- 0
      trans.mat[i,t,7,4] <- 0
      trans.mat[i,t,7,5] <- 0
      trans.mat[i,t,7,6] <- 0
      trans.mat[i,t,7,7] <- 1
      
      Z[i, t] ~ dcat(trans.mat[i, t, Z[i, t-1], 1:7])
      Y[i, t] ~ dcat(obs.mat[t, Z[i, t], 1:2])
    }
  }
  
  for (t in 1:n.occasions) {
    p.det[t] <- mean.p.det
    
    # alive at sea, never bred
    obs.mat[t, 1, 1] <- 1
    obs.mat[t, 1, 2] <- 0
    
    # alive at sea, has bred
    obs.mat[t, 2, 1] <- 1
    obs.mat[t, 2, 2] <- 0
    
    # recruit in colony
    obs.mat[t, 3, 1] <- 1-p.det[t]
    obs.mat[t, 3, 2] <- p.det[t]
    
    # recruit out of colony
    obs.mat[t, 4, 1] <- 1
    obs.mat[t, 4, 2] <- 0
    
    # breed in colony
    obs.mat[t, 5, 1] <- 1-p.det[t]
    obs.mat[t, 5, 2] <- p.det[t]
    
    # breed out of colony
    obs.mat[t, 6, 1] <- 1
    obs.mat[t, 6, 2] <- 0
    
    # dead
    obs.mat[t, 7, 1] <- 1
    obs.mat[t, 7, 2] <- 0
  }
  
})

#### PARAMETERS TO MONITOR ####
params <- c(
  # FECUNDITY
  "ann.fec",
  
  # COUNTS
  #"sigma.obs", 
  
  # SURVIVAL
  "mean.phi.juv", "mean.phi.ad", "sigma.phi", "eps.phi",
  "mu.p.juv", "mean.p.ad", "agebeta", "sigma.p", "eps.p",
  "mean.p.in",
  "mean.p.propensity",
  "mean.p.det",
  
  # ABUNDANCE
  "Ntot", "Ntot.breed", "N.atsea", "N.ad.surv",
  "N.breed.ready", "N.recruits", "nestlings", "JUV", "IM"
)

#### MCMC SETTINGS ####
nb <- 0 #burn-in
ni <- 200#000 #total iterations
nt <- 1 # thin
nc <- 3  #chains
adaptInterval = 200
#maxContractions = 1000
#scale = 1
#sliceWidth = 1

#### COMPILE CONFIGURE AND BUILD ####
t.start <- Sys.time()
Rmodel <- nimbleModel(code = code, constants = const, 
                      data = dat, 
                      inits = inits,
                      check = FALSE, calculate = FALSE)
conf <- configureMCMC(Rmodel, monitors = params, thin = nt, 
                      useConjugacy = FALSE,
                      control = list(adaptInterval = adaptInterval#,
                                     #maxContractions = maxContractions, 
                                     #scale = scale, 
                                     #sliceWidth = sliceWidth
                      )) 
Rmcmc <- buildMCMC(conf)  
Cmodel <- compileNimble(Rmodel, showCompilerOutput = FALSE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
t.end <- Sys.time()
(runTime <- t.end - t.start)
beep(sound = 1)


#### RUN MCMC ####

# with inits - does not mix

t.start <- Sys.time()
#sink("oops.txt")
out1 <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = nc, inits = inits,
                setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
#sink()
t.end <- Sys.time()
(runTime <- t.end - t.start)

rhat <- gelman.diag(out1, multivariate = FALSE)
beep(sound = 1)
tmp <- rhat$psrf %>% 
  as.data.frame() %>% 
  filter(is.infinite(`Point est.`)|is.nan(`Point est.`)) %>% 
  rownames_to_column() %>% filter(!str_detect(rowname, "IM")) %>% 
  filter(!(str_first_number(rowname) == 1))

# TODO
# check correlations and apply block samplers
# also add back informed priors

summ <- t(post_summ(out1, get_params(out1, type = "base_index"), 
                    neff = TRUE, Rhat = TRUE, probs = c(0.025, 0.5, 0.975))) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name")
beep(sound = 1)

IMsamps <- post_subset(out1, "IM", matrix = TRUE, chains = TRUE, iters = TRUE)
Nbreedreadysamps <- post_subset(out1, "N.breed.ready", matrix = TRUE, chains = TRUE, iters = TRUE)
Nrecruitssamps <- post_subset(out1, "N.recruits", matrix = TRUE, chains = TRUE, iters = TRUE)
Ntotbreed <- post_subset(out1, "Ntot.breed", matrix = TRUE, chains = TRUE, iters = TRUE)

allNodes <- Cmodel$getNodeNames()
monitoredNodes <- allNodes[allNodes %in% colnames(out1[[1]])]
stoch <- Cmodel$isStoch(monitoredNodes)
stochMonitoredNodes <- monitoredNodes[stoch]

corr_df <- cor(out1[[1]], out1[[1]], method = "pearson")
corr_df[lower.tri(corr_df, diag = T)] <- NA
corr_df <- corr_df %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(-1) %>% 
  filter(!is.na(value)) %>% 
  filter(!((abs(value) - 1) > -1e-05)) %>% 
  filter(abs(value) > 0.5) %>% 
  arrange(desc(abs(value))) %>% 
  filter((rowname %in% stochMonitoredNodes) & (name %in% stochMonitoredNodes))


