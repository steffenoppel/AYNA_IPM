#### LOAD LIBRARIES ####
library(nimble)
library(jagsUI)
library(coda)

#### LOAD DATA ####
load("IPM_AEB_dat.RData")

#### LOAD INITS #####
load("IPM_AEB_inits.Rdata")

#### MODEL CODE ####

sink("IPM_AEB_jags.txt")
cat("
model { 

  #-------------------------------------------------  
  # 1. PRIORS FOR ALL DATA SETS
  #-------------------------------------------------
  
  # -------------------------------------------------        
  # 1.1. Priors and constraints FOR FECUNDITY
  # -------------------------------------------------
  
  for (t in 1:n.years.fec){  
    ann.fec[t] ~ dbeta(32,68)
  } # t
  
  # -------------------------------------------------        
  # 1.2. Priors and constraints FOR POPULATION COUNTS
  # -------------------------------------------------
  for (s in 1:n.sites.count){		### start loop over every study area
    for (t in 1:n.years.fec){		### start loop over every year
      sigma.obs[s,t] ~ dexp(1)	#Prior for SD of observation process (variation in detectability)
    } # t
  } # s
  
  # -------------------------------------------------        
  # 1.3. Priors and constraints FOR SURVIVAL
  # -------------------------------------------------
  
  ### RECAPTURE PROBABILITY
  
  mean.p.ad[1] ~ dunif(0.05, 0.5)	      
  mean.p.ad[2] ~ dunif(0.2, 1)	       
  
  for (gy in 1:2){  
    mu.p.juv[gy] ~ dnorm(-4, 1/0.25) 
    mu.p.ad[gy] <- log(mean.p.ad[gy] / (1-mean.p.ad[gy])) 
  }
  
  agebeta ~ dnorm(1, 1/0.01) # Prior for shape of increase in juvenile recapture probability with age
  
  sigma.p ~ dexp(1) 
  
  ## RANDOM TIME EFFECT ON RESIGHTING PROBABILITY...
  for (t in 1:(n.occasions-1)){
    
    # ...FOR JUVENILES
    for (j in 1:t){ 
      p.juv[t,j] <- 0
    }
    for (j in (t+1):(n.occasions-1)){
      logit(p.juv[t,j])  <- mu.p.juv[goodyear[j]] + agebeta*(j - t)/2 + eps.p[j]
    }
    
    # ...FOR ADULTS
    logit(p.ad[t])  <- mu.p.ad[goodyear[t]]  + eps.p[t]  
    eps.p[t] ~ dnorm(0, 1/sigma.p)
  }
  
  ### SURVIVAL PROBABILITY
  mean.phi.juv ~ dbeta(75.7,24.3)            
  mean.phi.ad ~ dbeta(91,9)              
  
  mu.juv <- log(mean.phi.juv / (1-mean.phi.juv)) # Logit transformation
  mu.ad <- log(mean.phi.ad / (1-mean.phi.ad)) # Logit transformation
  
  sigma.phi ~ dexp(1) 
  
  ## RANDOM TIME EFFECT ON SURVIVAL
  for (j in 1:(n.occasions-1)){
    logit(phi.juv[j]) <- mu.juv + eps.phi[j]*juv.poss[j]
    logit(phi.ad[j]) <- mu.ad + eps.phi[j] 
    eps.phi[j] ~ dnorm(0, 1/sigma.phi) 
  }
  
  #-------------------------------------------------  
  # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
  #-------------------------------------------------
  
  # -------------------------------------------------        
  # 2.1. System process: female based matrix model
  # -------------------------------------------------
  
  ### RECRUIT PROBABILITY ###
  for (age in 1:30) {
    logit(p.juv.recruit.f[age])<-mu.p.juv[2] + (agebeta * age/2)
  }
  
  IM[1,1,1] ~ dnorm(263/2,1/20) T(0, )
  IM[1,1,2] <- 0
  IM[1,1,3] <- round(IM[1,1,1]) - IM[1,1,2]
  
  IM[1,2,1] ~ dnorm(275/2,1/20) T(0, )
  IM[1,2,2] ~ dbin(p.juv.recruit.f[2], round(IM[1,2,1]))
  IM[1,2,3] <- round(IM[1,2,1]) - IM[1,2,2]
  
  IM[1,3,1] ~ dnorm(264/2,1/20) T(0, )
  IM[1,3,2] ~ dbin(p.juv.recruit.f[3], round(IM[1,3,1]))
  IM[1,3,3] <- round(IM[1,3,1]) - IM[1,3,2]
  
  IM[1,4,1] ~ dnorm(177/2,1/20) T(0, )
  IM[1,4,2] ~ dbin(p.juv.recruit.f[4], round(IM[1,4,1]))
  IM[1,4,3] <- round(IM[1,4,1]) - IM[1,4,2]
  
  IM[1,5,1] ~ dnorm(290/2,1/20) T(0, )
  IM[1,5,2] ~ dbin(p.juv.recruit.f[5], round(IM[1,5,1]))
  IM[1,5,3] <- round(IM[1,5,1]) - IM[1,5,2]
  
  IM[1,6,1] ~ dnorm(90/2,1/20) T(0, )
  IM[1,6,2] ~ dbin(p.juv.recruit.f[6], round(IM[1,6,1]))
  IM[1,6,3] <- round(IM[1,6,1]) - IM[1,6,2]
  
  IM[1,7,1] ~ dnorm(158/2,1/20) T(0, )
  IM[1,7,2] ~ dbin(p.juv.recruit.f[7], round(IM[1,7,1]))
  IM[1,7,3] <- round(IM[1,7,1]) - IM[1,7,2]
  
  for(age in 8:30) {
    IM[1,age,1] ~ dbin(mean.phi.ad, IM[1,age-1,3])
    IM[1,age,2] ~ dbin(p.juv.recruit.f[age], IM[1,age,1])
    IM[1,age,3] <- IM[1,age,1] - IM[1,age,2]
  }
  
  N.recruits[1] <- sum(IM[1,1:30,2]) 
  N.ad.surv[1] <- 0
  N.breed.ready[1] <- 0
  Ntot.breed[1] ~ dnorm(640/2,1/20) T(0, )
  N.atsea[1] ~ dnorm(224/2,1/20) T(0, ) 
  nestlings[1] <- 0
  JUV[1] ~ dnorm(232/2, 1/20) T(0, )     
  
  Ntot[1]<-sum(IM[1,1:30,3]) + round(Ntot.breed[1])+round(N.atsea[1]) 
  
  ### FOR EVERY SUBSEQUENT YEAR POPULATION PROCESS
  
  for (tt in 2:n.years.fec){
    for (age in 1:30) {
      logit(p.juv.recruit[age,tt])<- mu.p.juv[2] + eps.p[tt+offset-1] + (agebeta / 2 * age) 
    }
    
    ## IMMATURE MATRIX WITH 3 columns:
    # 1: survivors from previous year
    # 2: recruits in current year
    # 3: unrecruited in current year (available for recruitment next year)
    
    IM[tt,1,1] ~ dbin(phi.juv[tt+offset-1], round(JUV[tt-1]))                                  
    IM[tt,1,2] <- 0 
    IM[tt,1,3] <- IM[tt,1,1] - IM[tt,1,2]
    for(age in 2:30) {
      IM[tt,age,1] ~ dbin(phi.ad[tt+offset-1], IM[tt-1,age-1,3])
      IM[tt,age,2] ~ dbin(p.juv.recruit[age,tt], IM[tt,age,1])
      IM[tt,age,3] <- IM[tt,age,1] - IM[tt,age,2]
    }
    
    N.recruits[tt] <- sum(IM[tt,1:30,2])  ### number of this years recruiters
    N.ad.surv[tt] ~ dbin(phi.ad[tt+offset-1], round(Ntot.breed[tt-1])+round(N.atsea[tt-1]))  ### previous year's adults that survive
    N.breed.ready[tt] ~ dbin(p.ad[tt+offset-1], N.ad.surv[tt]) ### number of available breeders is proportion of survivors that returns
    Ntot.breed[tt]<- N.breed.ready[tt]+N.recruits[tt]  ### number of counted breeders is sum of old breeders returning and first recruits
    N.atsea[tt] <- N.ad.surv[tt]-N.breed.ready[tt] ### potential breeders that remain at sea
    nestlings[tt] <- round(ann.fec[tt] * 0.5 * Ntot.breed[tt]) ### nestlings produced                                                  
    JUV[tt] ~ dpois(nestlings[tt]) ### juveniles produced
    
    ### THE TOTAL POPULATION ###
    Ntot[tt]<-sum(IM[tt,1:30,3]) + Ntot.breed[tt]+N.atsea[tt]  ## total population size is all the immatures plus adult breeders and adults at sea
    
  } # tt
  
  # -------------------------------------------------
  # 2.2. Observation process for population counts: state-space model of annual counts
  # -------------------------------------------------
  
  for (s in 1:n.sites.count){	
    for (t in 1:n.years.fec){
      y.count[t,s] ~ dnorm(Ntot.breed[t]*prop.sites[t,s], 1/sigma.obs[s,t])
    }	# t
  }	# s
  
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
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,1:n.occasions], r.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,1:n.occasions], r.a[t])
  }
  
  # Define the cell probabilities of the m-arrays
  # Main diagonal
  for (t in 1:(n.occasions-1)){
    q.ad[t] <- 1-p.ad[t]            # Probability of non-recapture
    
    for(j in 1:(n.occasions-1)){
      q.juv[t,j] <- 1 - p.juv[t,j]
    }
    
    pr.j[t,t] <- 0
    pr.a[t,t] <- phi.ad[t]*p.ad[t]
    
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- phi.juv[t]*prod(phi.ad[(t+1):j])*prod(q.juv[t,t:(j-1)])*p.juv[t,j]
      pr.a[t,j] <- prod(phi.ad[t:j])*prod(q.ad[t:(j-1)])*p.ad[j]
    } #j
    
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  } #t

}
", fill = TRUE)
sink()

#### DATA ####
dat.jags <- c(dat, const) # jags bundles data and constants

### INITS ####

inits.jags <- head(inits, 12) # jags doesn't let you initialize any deterministic nodes
inits.jags <- inits.jags[-6] 
# could be smarter about this but this is easy, for now
inits.jags <- list(inits.jags, inits.jags, inits.jags)

#### PARAMETERS TO MONITOR ####
params <- c(
  "Ntot", "Ntot.breed", "N.atsea", "N.ad.surv", 
  "N.breed.ready", "N.recruits", "nestlings", "JUV", "IM"
)

#### MCMC SETTINGS ####
nb <- 0 #burn-in
ni <- 200000 #total iterations
nt <- 1 # thin
nc <- 3  #chains
adaptInterval = 10000 # jags does adaptation different than nimble
maxContractions = 1000
scale = 1
sliceWidth = 1

# inits
t.start <- Sys.time()
out1 <- jagsUI::jags(data = dat.jags, inits = inits.jags, parameters.to.save = params, 
                     model.file = "IPM_AEB_jags.txt", n.chains = nc, n.thin = nt, 
                     n.iter = ni, n.burnin = nb, n.adapt = adaptInterval#, parallel = TRUE
)
t.end <- Sys.time()
(runTime <- t.end - t.start)

# no inits
t.start <- Sys.time()
out2 <- jagsUI::jags(data = dat.jags, parameters.to.save = params, 
                     model.file = "IPM_AEB_jags.txt", n.chains = nc, n.thin = nt, 
                     n.iter = ni, n.burnin = nb, n.adapt = adaptInterval#, parallel = TRUE
)
t.end <- Sys.time()
(runTime <- t.end - t.start)

t.start <- Sys.time()
out3 <- jagsUI::jags(data = dat.jags, inits = inits.jags, parameters.to.save = params, 
                     model.file = "IPM_AEB_jags.txt", n.chains = nc, n.thin = nt, 
                     n.iter = ni, n.burnin = nb, n.adapt = adaptInterval, parallel = TRUE
)
t.end <- Sys.time()
(runTime <- t.end - t.start)

# no inits
t.start <- Sys.time()
out4 <- jagsUI::jags(data = dat.jags, parameters.to.save = params, 
                     model.file = "IPM_AEB_jags.txt", n.chains = nc, n.thin = nt, 
                     n.iter = ni, n.burnin = nb, n.adapt = adaptInterval, parallel = TRUE
)
t.end <- Sys.time()
(runTime <- t.end - t.start)