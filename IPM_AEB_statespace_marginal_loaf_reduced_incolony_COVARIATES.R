#### LOAD LIBRARIES ####
library(nimble)
#library(jagsUI)
library(coda)
#library(doParallel) 
library(tidyverse)
library(tidybayes)
library(strex)
library(beepr)
library(postpack)
#library(nimbleEcology)
library(here)

# REVISE FUNCTIONS ####

dDHMMo_mod <- nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),##
                 probObs = double(3),
                 probTrans = double(3),
                 mult = double(0),
                 len = double(),## length of x (needed as a separate param for rDHMM)
                 checkRowSums = double(0, default = 1),
                 log = integer(0, default = 0)) {
    if (length(init) != dim(probObs)[1]) stop("In dDHMMo: Length of init does not match ncol of probObs in dDHMMo.")
    if (length(init) != dim(probTrans)[1]) stop("In dDHMMo: Length of init does not match dim(probTrans)[1] in dDHMMo.")
    if (length(init) != dim(probTrans)[2]) stop("In dDHMMo: Length of init does not match dim(probTrans)[2] in dDHMMo.")
    if (length(x) != len) stop("In dDHMMo: Length of x does not match len in dDHMM.")
    if (len - 1 > dim(probTrans)[3]) stop("In dDHMMo: dim(probTrans)[3] does not match len - 1 in dDHMMo.")
    if (len != dim(probObs)[3]) stop("In dDHMMo: dim(probObs)[3] does not match len in dDHMMo.")
    if (abs(sum(init) - 1) > 1e-6) stop("In dDHMMo: Initial probabilities must sum to 1.")
    
    if (checkRowSums) {
      transCheckPasses <- TRUE
      for (i in 1:dim(probTrans)[1]) {
        for (k in 1:dim(probTrans)[3]) {
          thisCheckSum <- sum(probTrans[i,,k])
          if (abs(thisCheckSum - 1) > 1e-6) {
            ## Compilation doesn't support more than a simple string for stop()
            ## so we provide more detail using a print().
            print("In dDHMMo: Problem with sum(probTrans[i,,k]) with i = ", i, " k = ", k, ". The sum should be 1 but is ", thisCheckSum)
            transCheckPasses <- FALSE
          }
        }
      }
      obsCheckPasses <- TRUE
      for (i in 1:dim(probObs)[1]) {
        for (k in 1:dim(probObs)[3]) {
          thisCheckSum <- sum(probObs[i,,k])
          if (abs(thisCheckSum - 1) > 1e-6) {
            print("In dDHMMo: Problem with sum(probObs[i,,k]) with i = ", i, " k = ", k, ". The sum should be 1 but is ", thisCheckSum)
            obsCheckPasses <- FALSE
          }
        }
      }
      if(!(transCheckPasses | obsCheckPasses))
        stop("In dDHMMo: probTrans and probObs were not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
      if(!transCheckPasses)
        stop("In dDHMMo: probTrans was not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
      if(!obsCheckPasses)
        stop("In dDHMMo: probObs was not specified correctly. Probabilities in each row must sum to 1.")
    }
    
    pi <- init # State probabilities at time t=1
    logL <- 0
    nObsClasses <- dim(probObs)[2]
    lengthX <- length(x)
    for (t in 1:lengthX) {
      if (x[t] > nObsClasses | x[t] < 1) stop("In dDHMMo: Invalid value of x[t].")
      Zpi <- probObs[, x[t], t] * pi # Vector of P(state) * P(observation class x[t] | state)
      sumZpi <- sum(Zpi)    # Total P(observed as class x[t])
      logL <- logL + log(sumZpi) * mult  # Accumulate log probabilities through time
      if (t != lengthX) pi <- ((Zpi %*% probTrans[,,t])/sumZpi)[1, ] # State probabilities at t+1
    }
    returnType(double())
    if (log) return(logL)
    return(exp(logL))
  }
)

rDHMMo_mod <- nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),
                 probObs = double(3),
                 probTrans = double(3),
                 mult = double(0),
                 len = double(),
                 checkRowSums = double(0, default = 1)) {
    nStates <- length(init)
    if (nStates != dim(probObs)[1]) stop("In rDHMMo: Length of init does not match nrow of probObs in dDHMM.")
    if (nStates != dim(probTrans)[1]) stop("In rDHMMo: Length of init does not match dim(probTrans)[1] in dDHMM.")
    if (nStates != dim(probTrans)[2]) stop("In rDHMMo: Length of init does not match dim(probTrans)[2] in dDHMM.")
    if (len - 1 > dim(probTrans)[3]) stop("In rDHMMo: len - 1 does not match dim(probTrans)[3] in dDHMM.")
    if (abs(sum(init) - 1) > 1e-6) stop("In rDHMMo: Initial probabilities must sum to 1.")
    if (checkRowSums) {
      transCheckPasses <- TRUE
      for (i in 1:dim(probTrans)[1]) {
        for (k in 1:dim(probTrans)[3]) {
          thisCheckSum <- sum(probTrans[i,,k])
          if (abs(thisCheckSum - 1) > 1e-6) {
            ## Compilation doesn't support more than a simple string for stop()
            ## so we provide more detail using a print().
            print("In rDHMMo: Problem with sum(probTrans[i,,k]) with i = ", i, " k = ", k, ". The sum should be 1 but is ", thisCheckSum)
            transCheckPasses <- FALSE
          }
        }
      }
      obsCheckPasses <- TRUE
      for (i in 1:dim(probObs)[1]) {
        for (k in 1:dim(probObs)[3]) {
          thisCheckSum <- sum(probObs[i,,k])
          if (abs(thisCheckSum - 1) > 1e-6) {
            print("In rDHMMo: Problem with sum(probObs[i,,k]) with i = ", i, " k = ", k, ". The sum should be 1 but is ", thisCheckSum)
            obsCheckPasses <- FALSE
          }
        }
      }
      if(!(transCheckPasses | obsCheckPasses))
        stop("In rDHMMo: probTrans and probObs were not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
      if(!transCheckPasses)
        stop("In rDHMMo: probTrans was not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
      if(!obsCheckPasses)
        stop("In rDHMMo: probObs was not specified correctly. Probabilities in each row must sum to 1.")
    }
    
    returnType(double(1))
    ans <- numeric(len)
    
    trueInit <- 0
    
    r <- runif(1, 0, 1)
    j <- 1
    while (r > sum(init[1:j])) j <- j + 1
    trueState <- j
    
    for (i in 1:len) {
      # Detect based on the true state
      r <- runif(1, 0, 1)
      j <- 1
      while (r > sum(probObs[trueState, 1:j, i])) j <- j + 1
      ans[i] <- j
      
      # Transition to a new true state
      if (i != len) {
        r <- runif(1, 0, 1)
        j <- 1
        while (r > sum(probTrans[trueState, 1:j, i])) j <- j + 1
        trueState <- j
      }
    }
    return(ans)
  })

#### LOAD DATA ####
load("IPM_AEB_dat_stateSpace_marginal_loaf_reduced_incolony_COVARIATES.RData") 

#### LOAD INITS #####
load("IPM_AEB_inits_stateSpace_marginal_loaf_reduced_incolony_COVARIATES.Rdata") 

#### MODEL CODE ####
code <- nimbleCode({
  
  #-------------------------------------------------  
  # 1. PRIORS FOR ALL DATA SETS
  #-------------------------------------------------
  
  # -------------------------------------------------        
  # 1.1. Priors and constraints FOR FECUNDITY
  # -------------------------------------------------
  
  mean.fec ~ dbeta(10,10)
  mu.fec <- log(mean.fec / (1-mean.fec)) # Logit transformation
  sigma.fec ~ dexp(10)
  for (t in 1:n.years.fec){  
    eps.fec[t] ~ dnorm(0, sd = sigma.fec)
    logit(ann.fec[t]) <- mu.fec + eps.fec[t] # AEB change, more variation in data
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
  
  # mean.p.fidelity[1] ~ dbeta(35/2,10/2) # immature fidelity
  # mean.p.fidelity[2] ~ dbeta(60,4) # adult fidelity
  # mean.p.fidelity ~ dbeta(60,4) # adult fidelity
  
  mean.p.propensity ~ dbeta(10,3) # adult breeding prob, if on land
  
  mean.p.atsea ~ dbeta(4,10) # prob of remaining at sea, breeders
  
  mean.p.det ~ dnorm(3,sd = 0.5) # probability of detecting an adult, given they are there
  
  #sigma.p.det ~ dexp(10)
  
  #for (t in 1:n.occasions) {
  #  eps.det[t] ~ dnorm(0, sd = sigma.p.det)
  #  logit(p.det[t]) <- mean.p.det + eps.det[t]
  #}
  
  ### RECAPTURE PROBABILITY
  #mean.p.ad[1] ~ dbeta(8.3, 25)	# low detection
  #mean.p.ad ~ dbeta(1, 1)	# low detection, # no longer really relevant, because we have propensity now
  
  mu.p.juv ~ dnorm(-3, sd = 0.25) 
  #mu.p.ad <- log(mean.p.ad / (1-mean.p.ad)) 
  
  # Prior for shape of increase in juvenile recapture probability with age
  #agebeta ~ dnorm(2, sd = 0.25) # AEB change
  agebeta ~ dnorm(0.75, sd = 0.1) # AEB change
  
  sigma.p ~ dexp(10) # AEB change, reduce variability 
  
  ## RANDOM TIME EFFECT ON RECRUITMENT PROBABILITY...
  for (t in 1:(n.occasions-1)){
    # ...FOR JUVENILES
    for (j in 1:t){ 
      p.juv[t,j] <- 0
    }
    for (j in (t+1):(n.occasions)){
      logit(p.juv[t,j])  <- mu.p.juv + agebeta*(j - t) #+ eps.p[j-1]
      
      #mean.p.recruit[t, j] <- p.juv[t, j]
      #mean.p.recruit[2, t, j] <- p.ad[t]
    }
    # ...FOR ADULTS
    #logit(p.ad[t])  <- mu.p.ad + eps.p[t]  
    eps.p[t] ~ dnorm(0, sd = sigma.p)
    logit(p.det[t]) <- mean.p.det + eps.p[t]
  }
  
  ### SURVIVAL PROBABILITY
  mean.phi.juv ~ dbeta(30,40) 
  mean.phi.imm ~ dbeta(40,10)
  #mean.phi.juv ~ dbeta(1, 1)            
  mean.phi.ad ~ dbeta(50,4)  
  #mean.phi.ad ~ dbeta(1,1)
  
  inflation.factor ~ dbeta(3, 50)
  
  mu.juv <- log(mean.phi.juv / (1-mean.phi.juv)) # Logit transformation
  mu.imm <- log(mean.phi.im) / (1-mean.phi.im))
  mu.ad <- log(mean.phi.ad / (1-mean.phi.ad)) # Logit transformation
  
  sigma.phi ~ dexp(10) # AEB changed, don't want this huge
  
  ## RANDOM TIME EFFECT ON SURVIVAL
  
  w1 ~ dbeta(1,1)
  w2 ~ dbeta(1,1)
  w3 ~ dbeta(1,1)
  w4 ~ dbeta(1,1)
  ones ~ dconstraint(w1 + w2 + w3 + w4 == 1)
  
  sigma.beta1 ~ dexp(10)
  # sigma.beta2 ~ dexp(1)
  # sigma.beta3 ~ dexp(1)
  # sigma.beta4 ~ dexp(1)
  # sigma.beta5 ~ dexp(1)
  
  # beta.ICCAT.ll.e[1] ~ dnorm(0, sd = sigma.beta1)  
  # beta.ICCAT.ll.mit[1] ~ dnorm(0, sd = sigma.beta2) 
  # beta.Nam.ll.mit[1] ~ dnorm(0, sd = sigma.beta3) 
  # beta.SA.ll.mit[1] ~ dnorm(0, sd = sigma.beta4) 
  # beta.Uru.ll.mit[1] ~ dnorm(0, sd = sigma.beta5)
  # #beta.ICCAT.ll.e[2] ~ dnorm(0, sd = sigma.beta1)  
  # beta.ICCAT.ll.mit[2] ~ dnorm(0, sd = sigma.beta2) 
  # beta.Nam.ll.mit[2] ~ dnorm(0, sd = sigma.beta3) 
  # beta.SA.ll.mit[2] ~ dnorm(0, sd = sigma.beta4) 
  # beta.Uru.ll.mit[2] ~ dnorm(0, sd = sigma.beta5) 
  beta.mit[1] ~ dnorm(0, sd = sigma.beta1)
  beta.mit[2] ~ dnorm(0, sd = sigma.beta1)
  
  for (j in 1:(n.occasions-1)){
    # logit(phi.juv[j]) <- mu.juv + eps.phi[j]*juv.poss[j] + beta.ICCAT.ll.e[1]*ICCAT.ll.e[j] + beta.ICCAT.ll.mit[1]*ICCAT.ll.mit[j] + beta.Nam.ll.mit[1]*Nam.ll.mit[j] + beta.SA.ll.mit[1]*SA.ll.mit[j] + beta.Uru.ll.mit[1]*Uru.ll.mit[j]
    # logit(phi.ad[j]) <- mu.ad + eps.phi[j] + beta.ICCAT.ll.e[2]*ICCAT.ll.e[j] + beta.ICCAT.ll.mit[2]*ICCAT.ll.mit[j] + beta.Nam.ll.mit[2]*Nam.ll.mit[j] + beta.SA.ll.mit[2]*SA.ll.mit[j] + beta.Uru.ll.mit[2]*Uru.ll.mit[j]
    
    mit.index[j] <- w1*ICCAT.ll.mit[j] + w2*Nam.ll.mit[j] + w3*SA.ll.mit[j] + w4*Uru.ll.mit[j]
    
    logit(phi.juv[j]) <- mu.juv + eps.phi[j] + beta.mit[1]*mit.index[j]
    logit(phi.im[j]) <-  mu.im  + eps.phi[j] + beta.mit[1]*mit.index[j]
    logit(phi.ad[j])  <- mu.ad  + eps.phi[j] + beta.mit[2]*mit.index[j]
    
    eps.phi[j] ~ dnorm(0, sd = sigma.phi) 
    
    mean.phi[1, j] <- phi.juv[j]
    mean.phi[2, j] <- phi.im[j] 
    mean.phi[3, j] <- phi.ad[j]
  }
  
  #-------------------------------------------------  
  # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
  #-------------------------------------------------
  
  # -------------------------------------------------        
  # 2.1. System process: female based matrix model
  # -------------------------------------------------
  
  ### RECRUIT PROBABILITY ###
  for (age in 1:maxAge) {
    logit(p.juv.recruit.f[age]) <- mu.p.juv + agebeta*(age)
  }
  
  #IM[1,1,1] ~ T(dnorm(263/2,sd = 20), 0, Inf)
  IM[1,1,1] ~ dpois(263/2*0.9)
  IM[1,1,2] <- 0
  #IM[1,1,3] <- round(IM[1,1,1]) - IM[1,1,2]
  IM[1,1,3] <- IM[1,1,1] - IM[1,1,2]
  
  #IM[1,2,1] ~ T(dnorm(275/2,sd = 20), 0, Inf)
  IM[1,2,1] ~ dpois(275/2*0.9)
  IM[1,2,2] ~ dbin(p.juv.recruit.f[2], round(IM[1,2,1]))
  #IM[1,2,3] <- round(IM[1,2,1]) - IM[1,2,2]
  IM[1,2,3] <- IM[1,2,1] - IM[1,2,2]
  
  #IM[1,3,1] ~ T(dnorm(264/2,sd = 20), 0, Inf)
  IM[1,3,1] ~ dpois(264/2*0.9)
  IM[1,3,2] ~ dbin(p.juv.recruit.f[3], round(IM[1,3,1]))
  #IM[1,3,3] <- round(IM[1,3,1]) - IM[1,3,2]
  IM[1,3,3] <- IM[1,3,1] - IM[1,3,2]
  
  #IM[1,4,1] ~ T(dnorm(177/2,sd = 20), 0, Inf) 
  IM[1,4,1] ~ dpois(177/2*0.9)
  IM[1,4,2] ~ dbin(p.juv.recruit.f[4], round(IM[1,4,1]))
  #IM[1,4,3] <- round(IM[1,4,1]) - IM[1,4,2]
  IM[1,4,3] <- IM[1,4,1] - IM[1,4,2]
  
  #IM[1,5,1] ~ T(dnorm(290/2,sd = 20), 0, Inf)
  IM[1,5,1] ~ dpois(290/2*0.9)
  IM[1,5,2] ~ dbin(p.juv.recruit.f[5], round(IM[1,5,1]))
  #IM[1,5,3] <- round(IM[1,5,1]) - IM[1,5,2]
  IM[1,5,3] <- IM[1,5,1] - IM[1,5,2]
  
  #IM[1,6,1] ~ T(dnorm(90/2,sd = 20), 0, Inf) 
  IM[1,6,1] ~ dpois(90/2*0.9)
  IM[1,6,2] ~ dbin(p.juv.recruit.f[6], round(IM[1,6,1]))
  #IM[1,6,3] <- round(IM[1,6,1]) - IM[1,6,2]
  IM[1,6,3] <- IM[1,6,1] - IM[1,6,2]
  
  #IM[1,7,1] ~ T(dnorm(158/2,sd = 20), 0, Inf)
  IM[1,7,1] ~ dpois(158/2*0.9)
  IM[1,7,2] ~ dbin(p.juv.recruit.f[7], round(IM[1,7,1]))
  #IM[1,7,3] <- round(IM[1,7,1]) - IM[1,7,2]
  IM[1,7,3] <- IM[1,7,1] - IM[1,7,2]
  
  for(age in 8:maxAge) {
    IM[1,age,1] ~ dbin(mean.phi.im*(1+inflation.factor), IM[1,age-1,3])
    IM[1,age,2] ~ dbin(p.juv.recruit.f[age], IM[1,age,1])
    IM[1,age,3] <- IM[1,age,1] - IM[1,age,2]
  }
  
  N.recruits[1] <- sum(IM[1,1:maxAge,2]) 
  N.ad.surv[1] <- 0
  N.breed.ready[1] <- 0
  #Ntot.breed[1] ~ T(dnorm(640/2,sd = 20), 0, Inf)  
  Ntot.breed[1] ~ dpois(640/2*0.9)  
  #N.atsea[1] ~ T(dnorm(224/2,sd = 20), 0, Inf)    
  N.atsea[1] ~ dpois(40/2*0.9) # 640*0.05/2, with a little wiggle
  N.loaf[1] ~ dpois(120/2*0.9) # taken from count of loafers in 2010
  nestlings[1] <- 0
  #JUV[1] ~ T(dnorm(232/2, sd = 20), 0, Inf)         
  JUV[1] ~ dpois(232/2*0.9)
  
  #Ntot[1]<-sum(IM[1,1:maxAge,3]) + round(Ntot.breed[1])+round(N.atsea[1]) 
  Ntot[1]<-sum(IM[1,1:maxAge,3]) + Ntot.breed[1] + N.atsea[1] + N.loaf[1]
  
  ### FOR EVERY SUBSEQUENT YEAR POPULATION PROCESS
  
  for (tt in 2:n.years.fec){
    for (age in 1:maxAge) {
      logit(p.juv.recruit[age,tt]) <- mu.p.juv + agebeta*(age) # + eps.p[tt+offset-1]
    }
    
    ## IMMATURE MATRIX WITH 3 columns:
    # 1: survivors from previous year
    # 2: recruits in current year
    # 3: unrecruited in current year (available for recruitment next year)
    
    #IM[tt,1,1] ~ dbin(phi.juv[tt+offset-1], round(JUV[tt-1]))                                
    IM[tt,1,1] ~ dbin(phi.juv[tt+offset-1]*(1+inflation.factor), JUV[tt-1])                                
    IM[tt,1,2] <- 0 
    IM[tt,1,3] <- IM[tt,1,1] - IM[tt,1,2]
    for(age in 2:maxAge) {
      IM[tt,age,1] ~ dbin(phi.im[tt+offset-1]*(1+inflation.factor), IM[tt-1,age-1,3])
      IM[tt,age,2] ~ dbin(p.juv.recruit[age,tt], IM[tt,age,1])
      IM[tt,age,3] <- IM[tt,age,1] - IM[tt,age,2]
    }
    
    N.recruits[tt] <- sum(IM[tt,1:maxAge,2])  ### number of this years recruiters
    #N.ad.surv[tt] ~ dbin(phi.ad[tt+offset-1], round(Ntot.breed[tt-1])+round(N.atsea[tt-1]))  ### previous year's adults that survive
    N.ad.surv[tt] ~ dbin(phi.ad[tt+offset-1], Ntot.breed[tt-1]+N.atsea[tt-1]+N.loaf[tt-1])  ### previous year's adults that survive
    N.breed.ready[tt] ~ dbin(1-mean.p.atsea, N.ad.surv[tt]) ### number of available breeders is proportion of survivors that returns
    N.atsea[tt] <- N.ad.surv[tt]-N.breed.ready[tt] ### potential breeders that remain at sea
    N.loaf[tt] ~ dbin(1-mean.p.propensity, N.breed.ready[tt])
    Ntot.breed[tt]<- N.breed.ready[tt]-N.loaf[tt]+N.recruits[tt]  ### number of counted breeders is sum of old breeders returning and first recruits
    nestlings[tt] <- round(ann.fec[tt] * 0.5 * Ntot.breed[tt]) ### nestlings produced, half are female                                                 
    JUV[tt] ~ dpois(nestlings[tt]) ### juveniles produced
    
    ### THE TOTAL POPULATION ###
    Ntot[tt]<-sum(IM[tt,1:maxAge,3]) + Ntot.breed[tt] + N.atsea[tt] + N.loaf[tt] ## total population size is all the immatures plus adult breeders and adults at sea
    
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
  # bummer that we can't estimate detection probability, though
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
    # Z[i, t] ~ dcat(trans.mat[i, t, Z[i, t-1], 1:7])
    # Y[i, t] ~ dcat(obs.mat[t, Z[i, t], 1:2])
    
    Y[i, first[i]:n.occasions] ~ dDHMMo_mod(init = init[i, 1:6], # initial state probabilities, provide this as data (1s and 0s, assuming first state known)
                                            probObs = obs.mat[i, 1:6, 1:2, first[i]:n.occasions], # time dependent 3d array (index by i) [nstates, nevents, t]; assume known at first and provide as data
                                            probTrans = trans.mat[i, 1:6, 1:6, (first[i]+1):n.occasions], # time dependent 3d array (index by i) [nstates, nstates, t]
                                            mult = mult[i],
                                            len = n.occasions - first[i] + 1, # length of observations
                                            checkRowSums = 0) 
    
    obs.mat[i,1:6, 1:2, first[i]] <- obs.mat.init[i, 1:6, 1:2]
    
    for (t in (first[i]+1):n.occasions) {
      # phi[i,t] related to individual age and year per previous model
      phi[i,t] <- mean.phi[ageCat[i,t-1], t-1] 
      
      #mean.p.recruit[t, j] <- p.juv[t, j]
      #mean.p.recruit[2, t, j] <- p.ad[t]
      
      #p.recruit[i,t] <- mean.p.recruit[ageCat[i,t-1], t-1, t]
      p.recruit[i,t] <- p.juv[t-1, t]
      
      #p.fidelity[i,t] <- mean.p.fidelity
      p.propensity[i,t] <- mean.p.propensity
      p.atsea[i,t] <- mean.p.atsea
      
      # 1. ALIVE AT SEA - NEVER BRED
      trans.mat[i,1,1,t] <- phi[i,t]*(1-p.recruit[i,t])
      trans.mat[i,1,2,t] <- 0
      trans.mat[i,1,3,t] <- 0
      trans.mat[i,1,4,t] <- phi[i,t]*(p.recruit[i,t])
      trans.mat[i,1,5,t] <- 0
      trans.mat[i,1,6,t] <- 1-phi[i,t]
      
      # 2. ALIVE AT SEA - LAST BRED IN COLONY
      trans.mat[i,2,1,t] <- 0 
      trans.mat[i,2,2,t] <- mean.phi[3,t-1]*p.atsea[i,t]
      trans.mat[i,2,3,t] <- mean.phi[3,t-1]*(1-p.atsea[i,t])*(1-p.propensity[i,t])
      trans.mat[i,2,4,t] <- 0
      trans.mat[i,2,5,t] <- mean.phi[3,t-1]*(1-p.atsea[i,t])*(p.propensity[i,t])
      trans.mat[i,2,6,t] <- 1-mean.phi[3,t-1]
      
      # 3. LOAF IN - HAS BRED
      trans.mat[i,3,1,t] <- 0 
      trans.mat[i,3,2,t] <- mean.phi[3,t-1]*p.atsea[i,t]
      trans.mat[i,3,3,t] <- mean.phi[3,t-1]*(1-p.atsea[i,t])*(1-p.propensity[i,t])
      trans.mat[i,3,4,t] <- 0
      trans.mat[i,3,5,t] <- mean.phi[3,t-1]*(1-p.atsea[i,t])*(p.propensity[i,t])
      trans.mat[i,3,6,t] <- 1-mean.phi[3,t-1]
      
      # 4. RECRUIT IN
      trans.mat[i,4,1,t] <- 0 
      trans.mat[i,4,2,t] <- mean.phi[3,t-1]*p.atsea[i,t]
      trans.mat[i,4,3,t] <- mean.phi[3,t-1]*(1-p.atsea[i,t])*(1-p.propensity[i,t])
      trans.mat[i,4,4,t] <- 0
      trans.mat[i,4,5,t] <- mean.phi[3,t-1]*(1-p.atsea[i,t])*(p.propensity[i,t])
      trans.mat[i,4,6,t] <- 1-mean.phi[3,t-1]
      
      # 5. BREED IN
      trans.mat[i,5,1,t] <- 0 
      trans.mat[i,5,2,t] <- mean.phi[3,t-1]*p.atsea[i,t]
      trans.mat[i,5,3,t] <- mean.phi[3,t-1]*(1-p.atsea[i,t])*(1-p.propensity[i,t])
      trans.mat[i,5,4,t] <- 0
      trans.mat[i,5,5,t] <- mean.phi[3,t-1]*(1-p.atsea[i,t])*(p.propensity[i,t])
      trans.mat[i,5,6,t] <- 1-mean.phi[3,t-1]
      
      # 6. DEAD
      trans.mat[i,6,1,t] <- 0 
      trans.mat[i,6,2,t] <- 0
      trans.mat[i,6,3,t] <- 0
      trans.mat[i,6,4,t] <- 0
      trans.mat[i,6,5,t] <- 0
      trans.mat[i,6,6,t] <- 1
      
      # detection!!!!!!!
      
      # alive at sea, never bred
      obs.mat[i,1, 1,t] <- 1
      obs.mat[i,1, 2,t] <- 0
      
      # alive at sea, last bred in
      obs.mat[i,2, 1,t] <- 1
      obs.mat[i,2, 2,t] <- 0
      
      # loaf in - has bred
      obs.mat[i,3, 1,t] <- 1-p.det[t-1]
      obs.mat[i,3, 2,t] <- p.det[t-1]
      
      # recruit in colony
      obs.mat[i,4, 1,t] <- 1-p.det[t-1]
      obs.mat[i,4, 2,t] <- p.det[t-1]
      
      # breed in colony
      obs.mat[i,5, 1,t] <- 1-p.det[t-1]
      obs.mat[i,5, 2,t] <- p.det[t-1]
      
      # dead
      obs.mat[i,6, 1,t] <- 1
      obs.mat[i,6, 2,t] <- 0
    }
  }
  
})

#### PARAMETERS TO MONITOR ####

params <- c(
  # FECUNDITY
  "mean.fec",
  "sigma.fec",
  "eps.fec",
  
  # SURVIVAL
  "mean.p.propensity", 
  "mean.p.atsea",
  "mean.p.det",
  "mu.p.juv", 
  "agebeta", 
  "sigma.p", 
  "eps.p",
  "mean.phi.juv", 
  "mean.phi.im",
  "mean.phi.ad", 
  "inflation.factor",
  "sigma.phi", 
  "eps.phi",
  
  # ABUNDANCE
  "IM", 
  "N.recruits", 
  "N.ad.surv", 
  "N.breed.ready",
  "Ntot.breed", 
  "N.atsea", 
  "N.loaf",
  "nestlings", 
  "JUV", 
  "Ntot",
  
  "sigma.beta1",
  #"sigma.beta2",
  #"sigma.beta3",
  #"sigma.beta4",
  #"sigma.beta5",
  
  "w1", "w2", "w3", "w4",
  "beta.mit"
  # "beta.ICCAT.ll.e",
  # "beta.ICCAT.ll.mit",
  # "beta.Nam.ll.mit",
  # "beta.SA.ll.mit",
  # "beta.Uru.ll.mit"
)

#### MCMC SETTINGS ####
nb <- 0 #burn-in
ni <- 100000 #total iterations
nt <- 1 # thin
nc <- 1  #chains
adaptInterval = 200
#maxContractions = 1000
#scale = 1
#sliceWidth = 1

#### COMPILE CONFIGURE AND BUILD ####

t.start <- Sys.time()
Rmodel <- nimbleModel(code = code, 
                      constants = const_marginal, 
                      data = dat_marginal, 
                      inits = inits_marginal,
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
out1 <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = nc, inits = inits_marginal,
                setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
#sink()
t.end <- Sys.time()
(runTime <- t.end - t.start)

save(out1, file = "samples_statespace_marginal_loaf_reduced_incolony_COVARIATES_chain1.Rdata")

