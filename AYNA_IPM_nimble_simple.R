#### COMMENTS SECTION ####

# TODO
# remove hardcoded sections

#### LOAD LIBRARIES #####
library(nimble)
library(here)
library(coda)
library(tidyverse)
library(lubridate)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select

#########################################################################
# LOAD PRE-PREPARED DATA ON COUNTS AND BREEDING SUCCESS
#########################################################################
### see 'IPM_DATA_PREPARATION_AYNA.R' for details on how data are aggregated

## LOAD PREPARED M-ARRAY FOR SURVIVAL ESTIMATION
load("AYNA_IPM_input.marray.RData")

## BOTH ARRAYS MUST HAVE EXACT SAME DIMENSIONS
dim(chick.marray)
dim(adult.marray)

### COUNT DATA FOR POPULATION TREND ######
head(POPSIZE)
names(POPSIZE)

POP<- as.matrix(POPSIZE[,2:12])
n.years.count<-nrow(POP)
n.sites.count<-ncol(POP)

#### BREEDING SUCCESS DATA FOR FECUNDITY ######
J<- as.matrix(CHICKCOUNT[,2:5])
R<- as.matrix(ADCOUNT[,2:5])

### specify constants for JAGS
n.years.fec<-dim(R)[1]		## defines the number of years
n.sites.fec<-dim(R)[2]    ## defines the number of study areas

### reduce R and J to vectors of sum across the study areas for which we have data
## will ensure appropriate weighting of breeding success by n pairs in each study area
# Area 10 has twice as many pairs as other areas

Jlong<-CHICKCOUNT %>% gather(key='Site', value="chicks",-Year)
PROD.DAT<-ADCOUNT %>% gather(key='Site', value="adults",-Year) %>%
  left_join(Jlong, by=c("Year","Site")) %>%
  mutate(include=ifelse(is.na(adults+chicks),0,1)) %>%
  filter(include==1) %>%
  group_by(Year) %>%
  summarise(J=sum(chicks),R=sum(adults))

### DIMENSION MISMATCH IN DATA
# IPM runs from 2008-2021 
# survival analysis runs from 1985-2020, but recapture index refers to columns, which represent year 1986-2021 plus the ones never recaptured (last column)
names(AYNA_CHICK)
POPSIZE$Year

OFFSET<-min(which(!is.na(match(as.numeric(substr(names(AYNA_CHICK)[2:44],1,4)),POPSIZE$Year))))
substr(names(AYNA_CHICK),1,4)[OFFSET+1]

#### MODEL CODE ####
code <- nimbleCode({

  #-------------------------------------------------  
  # 1. PRIORS FOR ALL DATA SETS
  #-------------------------------------------------
  
  # -------------------------------------------------        
  # 1.1. Priors and constraints FOR FECUNDITY
  # -------------------------------------------------
  
  for (t in 1:n.years.fec){  
    ann.fec[t] ~ dbeta(32,68)         ## Informative Priors on fecundity based on Wanless et al 2009
  } #t
  
  # -------------------------------------------------        
  # 1.2. Priors and constraints FOR POPULATION COUNTS
  # -------------------------------------------------
  for (s in 1:n.sites.count){			### start loop over every study area
    for (t in 1:n.years.fec){			### start loop over every year
      sigma.obs[s,t] ~ dexp(0.1)	#Prior for SD of observation process (variation in detectability)
    }
  }
  
  # -------------------------------------------------        
  # 1.3. Priors and constraints FOR SURVIVAL
  # -------------------------------------------------
  
  ### RECAPTURE PROBABILITY
  mean.p.ad[1] ~ dunif(0.05, 0.5)	           # Prior for mean adult recapture - should be higher than 5% but less than 50%
  mean.p.ad[2] ~ dunif(0.2, 1)	           # Prior for mean adult recapture - should be higher than 20%
  
  for (gy in 1:2){    ## for good and poor monitoring years
    mu.p.juv[gy] ~ dnorm(-4, sd = 0.25) # Logit scale prior for mean juvenile recapture - should be almost 0 at age 1 and increase with age/2
    mu.p.ad[gy] <- log(mean.p.ad[gy] / (1-mean.p.ad[gy])) # Logit transformation
  }
  
  agebeta ~ dnorm(1, sd = 0.001)    # Prior for shape osf increase in juvenile recapture probability with age
  
  ## RANDOM TIME EFFECT ON RESIGHTING PROBABILITY OF JUVENILES
  
  # TODO
  # this indexing seems wrong
  for (t in 1:(n.occasions-1)){
    for (j in 1:t){ ## zero by definition (these are never actually used)
      p.juv[t,j] <- 0
    }
    for (j in (t+1):(n.occasions-1)){
      logit(p.juv[t,j])  <- mu.p.juv[goodyear[j]] + agebeta*(j - t)/2 + eps.p[j]
    }
  }
  
  ## PRIORS FOR RANDOM EFFECTSf
  sigma.p ~ dexp(1)                # Prior for standard deviation
  
  ### SURVIVAL PROBABILITY
  mean.phi.juv ~ dbeta(75.7,24.3)             # Prior for mean juvenile survival first year 0.757, second year 0.973 in Laysan albatross
  mean.phi.ad ~ dbeta(91,9)              # Prior for mean adult survival - should be higher than 70%
  mu.juv <- log(mean.phi.juv / (1-mean.phi.juv)) # Logit transformation
  mu.ad <- log(mean.phi.ad / (1-mean.phi.ad)) # Logit transformation
  
  ## PRIORS FOR RANDOM EFFECTS
  sigma.phi ~ dexp(1)                # Prior for standard deviation
  
  ## RANDOM TIME EFFECT ON SURVIVAL AND ADULT RECAPTURE
  
  for (j in 1:(n.occasions-1)){
    logit(phi.juv[j]) <- mu.juv + eps.phi[j]*juv.poss[j] #+ beta.ICCAT.ll.e*ICCAT.ll.e[j] + beta.ICCAT.ll.mit*ICCAT.ll.mit[j] + beta.Nam.ll.mit*Nam.ll.mit[j] + beta.SA.ll.mit*SA.ll.mit[j] + beta.Uru.ll.mit*Uru.ll.mit[j]
    logit(phi.ad[j]) <- mu.ad + eps.phi[j] #+ beta.ICCAT.ll.e*ICCAT.ll.e[j] + beta.ICCAT.ll.mit*ICCAT.ll.mit[j] + beta.Nam.ll.mit*Nam.ll.mit[j] + beta.SA.ll.mit*SA.ll.mit[j] + beta.Uru.ll.mit*Uru.ll.mit[j]
    eps.phi[j] ~ dnorm(0, sd = sigma.phi) 
    logit(p.ad[j])  <- mu.p.ad[goodyear[j]] + eps.p[j]    #### CAT HORSWILL SUGGESTED TO HAVE A CONTINUOUS EFFORT CORRECTION: mu.p.ad + beta.p.eff*goodyear[j] + eps.p[j]
    eps.p[j] ~ dnorm(0, sd = sigma.p)
  }
  
  #-------------------------------------------------  
  # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
  #-------------------------------------------------
  
  # -------------------------------------------------        
  # 2.1. System process: female based matrix model
  # -------------------------------------------------
  
  ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1 - based on deterministic multiplications
  ## ADJUSTED BASED ON PAST POPULATION SIZES + BREEDING SUCCESS IN AREAS WITH COUNTS SINCE 2003
  ## BASED ON WANLESS PAPER, JUVENILES SURVIVE ON AVERAGE WITH RATE 0.757, ADULTS 0.973
  ## GIVES ROUGH ESTIMATE OF EXPECTED NUMBER IN EACH AGE CLASS
  ## CALCULATIONS ARE EST POP SIZE * EST BREEDING SUCCESS  * EST JUV SURVIVAL * EST ADULT SURVIVAL^N.YEARS
  ## TODO - recall it may be useful to start these really high (but same proportional relationships)
  
  # TODO - should it be a binomial function on the twos
  
  IMnonround[1] ~ T(dnorm(263/2,sd = 20), 0, Inf)   ### number of 1-year old survivors in 2007 (700*0.5*0.75) - CAN BE MANIPULATED
  IM[1,1,1] <- round(IMnonround[1])
  IM[1,1,2] <- 0
  IM[1,1,3] <- IM[1,1,1] - IM[1,1,2]
  
  IMnonround[2] ~ T(dnorm(275/2,sd = 20), 0, Inf)   ### number of 2-year old survivors in 2006 (680*0.6*0.75*0.9) CAN BE MANIPULATED
  IM[1,2,1] <- round(IMnonround[2])
  IM[1,2,2] ~ dbin(p.juv.recruit.f[2], IM[1,2,1])
  #IM[1,2,2] <- round(IM[1,2,1]*p.juv.recruit.f[2])
  IM[1,2,3] <- IM[1,2,1] - IM[1,2,2]
  
  IMnonround[3] ~ T(dnorm(264/2,sd = 20), 0, Inf)   ### number of 3-year old survivors in 2005 (680*0.64*0.75*0.9^2) - CAN BE MANIPULATED
  IM[1,3,1] <- round(IMnonround[3])
  #IM[1,3,2] <- round(IM[1,3,1]*p.juv.recruit.f[3])
  IM[1,3,2] ~ dbin(p.juv.recruit.f[3], IM[1,3,1])
  IM[1,3,3] <- IM[1,3,1] - IM[1,3,2]
  
  IMnonround[4] ~ T(dnorm(177/2,sd = 20), 0, Inf)   ### number of 4-year old survivors in 2004 (540*0.6*0.75*0.9^3) - CAN BE MANIPULATED
  IM[1,4,1] <- round(IMnonround[4])
  #IM[1,4,2] <- round(IM[1,4,1]*p.juv.recruit.f[4])
  IM[1,4,2] ~ dbin(p.juv.recruit.f[4], IM[1,4,1])
  IM[1,4,3] <- IM[1,4,1] - IM[1,4,2]
  
  IMnonround[5] ~ T(dnorm(290/2,sd = 20), 0, Inf)   ### number of 5-year old survivors in 2003 (709*0.83*0.75*0.9^4) - CAN BE MANIPULATED
  IM[1,5,1] <- round(IMnonround[5])
  #IM[1,5,2] <- round(IM[1,5,1]*p.juv.recruit.f[5])
  IM[1,5,2] ~ dbin(p.juv.recruit.f[5], IM[1,5,1])
  IM[1,5,3] <- IM[1,5,1] - IM[1,5,2]
  
  IMnonround[6] ~ T(dnorm(90/2,sd = 20), 0, Inf)    ### number of 6-year old survivors in 2002 (600*0.34*0.75*0.9^5) - CAN BE MANIPULATED
  IM[1,6,1] <- round(IMnonround[6])
  #IM[1,6,2] <- round(IM[1,6,1]*p.juv.recruit.f[6])
  IM[1,6,2] ~ dbin(p.juv.recruit.f[6], IM[1,6,1])
  IM[1,6,3] <- IM[1,6,1] - IM[1,6,2]
  
  IMnonround[7] ~ T(dnorm(158/2,sd = 20), 0, Inf)   ### number of 7-year old survivors in 2001 (650*0.61*0.75*0.9^6) - CAN BE MANIPULATED
  IM[1,7,1] <- round(IMnonround[7])
  #IM[1,7,2] <- round(IM[1,7,1]*p.juv.recruit.f[7])
  IM[1,7,2] ~ dbin(p.juv.recruit.f[7], IM[1,7,1])
  IM[1,7,3] <- IM[1,7,1] - IM[1,7,2]
  
  for(age in 8:30) {
    IM[1,age,1] ~ dbin(pow(mean.phi.ad,(age-1)), IM[1,age-1,3])
    #IM[1,age,2] <- round(IM[1,age,1]*p.juv.recruit.f[age])
    IM[1,age,2] ~ dbin(p.juv.recruit.f[age], IM[1,age,1])
    IM[1,age,3] <- IM[1,age,1] - IM[1,age,2]
  }
  N.recruits[1] <- sum(IM[1,1:30,2])  ### number of this years recruiters - irrelevant in year 1 as already included in Ntot.breed prior
  
  # TODO
  # are the correct things divided by 2
  Ntot.breednonround ~ T(dnorm(640/2,sd = 20), 0, Inf)  ### sum of counts is 640 ( sum(POP[1, ]) <- across 11 study areas )
  Ntot.breed[1] <- round(Ntot.breednonround)
  JUVnonround ~ T(dnorm(232/2, sd = 20), 0, Inf)          ### sum of chicks is 232 ( sum(mean.props[1, c(1,2,4,5)]) <- only 4 study areas counted, so correct w proportion )
  JUV[1] <- round(JUVnonround)
  N.atseanonround ~ T(dnorm(224/2,sd = 20), 0, Inf)    ### unknown number, but assume about 65% breeding each year per Cuthbert 2003 (sum(POP[1, ]) * (1-0.65)) - CAN BE MANIPULATED
  N.atsea[1] <- round(N.atseanonround)
  Ntot[1]<-sum(IM[1,1:30,3]) + Ntot.breed[1]+N.atsea[1]  ## total population size is all the immatures plus adult breeders and adults at sea - does not include recruits in Year 1
  
  
  ### FOR EVERY SUBSEQUENT YEAR POPULATION PROCESS
  
  ## recruit probability
  for (age in 1:30) {
    logit(p.juv.recruit.f[age])<-mu.p.juv[2] + (agebeta * age/2)
  }
  
  for (tt in 2:n.years.fec){
    
    ## THE PRE-BREEDING YEARS ##
    
    ## define recruit probability for various ages ##
    for (age in 1:30) {
      #logit(p.juv.recruit[age,tt])<-mu.p.juv[2] + eps.p[tt+offset-1] + (agebeta * age)
      logit(p.juv.recruit[age,tt])<- mu.p.juv[2] + eps.p[tt+offset-1] + (agebeta / 2 * age) # changed to be divided by two
    }
    
    ## IMMATURE MATRIX WITH 3 columns:
    # 1: survivors from previous year
    # 2: recruits in current year
    # 3: unrecruited in current year (available for recruitment next year)
    
    nestlings[tt] <- round(ann.fec[tt] * 0.5 * Ntot.breed[tt])                                                    ### number of locally produced FEMALE chicks
    ### CHECK HERE
    JUV[tt] ~ dpois(nestlings[tt])                                                                     ### need a discrete number otherwise dbin will fail, dpois must be >0
    IM[tt,1,1] ~ dbin(phi.juv[tt+offset-1], JUV[tt-1])                                  ### number of 1-year old survivors
    ### END CHECK
    
    IM[tt,1,2] <- 0 
    IM[tt,1,3] <- IM[tt,1,1] - IM[tt,1,2]
    
    for(age in 2:30) {
      IM[tt,age,1] ~ dbin(phi.ad[tt+offset-1], IM[tt-1,age-1,3])
      # TODO
      # what is going on here
      #IM[tt,age,2] <- min(round(IM[tt,age-1,3]),IM[tt,age,1])*p.juv.recruit[age,tt]
      #IM[tt,age,2] <- IM[tt,age,1]*p.juv.recruit[age,tt]
      IM[tt,age,2] ~ dbin(p.juv.recruit[age,tt], IM[tt,age,1])
      IM[tt,age,3] <- IM[tt,age,1] - IM[tt,age,2]
    }
    N.recruits[tt] <- sum(IM[tt,1:30,2])  ### number of this years recruiters
    
    ## THE BREEDING POPULATION ##
    # Ntot.breed comprised of first-time breeders, previous skippers, and previous unsuccessful breeders
    # simplified in simplified_v2 to just adult survivors with p.ad as proportion returning
    
    ## CAREFUL HERE TO ADD OFFSET SUCH THAT SURVIVAL YEARS ALIGN WITH COUNT YEARS
    
    N.ad.surv[tt] ~ dbin(phi.ad[tt+offset-1], Ntot.breed[tt-1]+N.atsea[tt-1])           ### previous year's adults that survive
    N.breed.ready[tt] ~ dbin(p.ad[tt+offset-1], N.ad.surv[tt])                  ### number of available breeders is proportion of survivors that returns
    Ntot.breed[tt]<- N.breed.ready[tt]+N.recruits[tt]             ### number of counted breeders is sum of old breeders returning and first recruits
    N.atsea[tt] <- N.ad.surv[tt]-N.breed.ready[tt]                    ### potential breeders that remain at sea
    
    ### THE TOTAL AYNA POPULATION ###
    Ntot[tt]<-sum(IM[tt,1:30,3]) + Ntot.breed[tt]+N.atsea[tt]  ## total population size is all the immatures plus adult breeders and adults at sea
    
  } # tt
  
  # -------------------------------------------------
  # 2.2. Observation process for population counts: state-space model of annual counts
  # -------------------------------------------------
  
  for (s in 1:n.sites.count){			### start loop over every study area
    
    ## Observation process
    
    for (t in 1:n.years.fec){
      # TODO - consider lognormal here because of low counts in some site-years
      # TODO - could also eliminate site loop here instead
      y.count[t,s] ~ dnorm(Ntot.breed[t]*prop.sites[t,s], sd = sigma.obs[s,t])								# Distribution for random error in observed numbers (counts)
    }														# run this loop over t= nyears
  }		## end site loop
  
  
  # -------------------------------------------------        
  # 2.3. Likelihood for fecundity: Logistic regression from the number of surveyed broods
  # -------------------------------------------------
  #for (s in 1:n.sites.fec){			### start loop over every study area
  # for (t in 1:(n.years.fec)){ # changed from minus 1 AEB
  #   J[t] ~ dbin(ann.fec[t], R[t])
  # } #	close loop over every year in which we have fecundity data
  #}
  
  # -------------------------------------------------        
  # 2.4. Likelihood for adult and juvenile survival from CMR
  # -------------------------------------------------
  # 
  # # Define the multinomial likelihood
  # for (t in 1:(n.occasions-1)){
  #   marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,1:n.occasions], r.j[t])
  #   marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,1:n.occasions], r.a[t])
  # }
  # 
  # # Define the cell probabilities of the m-arrays
  # # Main diagonal
  # for (t in 1:(n.occasions-1)){
  #   q.ad[t] <- 1-p.ad[t]            # Probability of non-recapture
  #   
  #   for(j in 1:(n.occasions-1)){
  #     q.juv[t,j] <- 1 - p.juv[t,j]
  #   }
  #   
  #   pr.j[t,t] <- 0
  #   pr.a[t,t] <- phi.ad[t]*p.ad[t]
  #   
  #   # Above main diagonal
  #   for (j in (t+1):(n.occasions-1)){
  #     pr.j[t,j] <- phi.juv[t]*prod(phi.ad[(t+1):j])*prod(q.juv[t,t:(j-1)])*p.juv[t,j]
  #     pr.a[t,j] <- prod(phi.ad[t:j])*prod(q.ad[t:(j-1)])*p.ad[j]
  #   } #j
  #   
  #   # Below main diagonal
  #   for (j in 1:(t-1)){
  #     pr.j[t,j] <- 0
  #     pr.a[t,j] <- 0
  #   } #j
  # } #t
  # 
  # # Last column: probability of non-recapture
  # for (t in 1:(n.occasions-1)){
  #   pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
  #   pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  # } #t

})

#### DATA ####
dat <- list(#marr.j = chick.marray,
            #marr.a = adult.marray,
            
            y.count=POP    ### use log(R) here if using the logscale model
            
            ### breeding success data
            #J=PROD.DAT$J,
            #R=PROD.DAT$R
            
            ## future fecundity change - vector with one element for each scenario
) 

#### CONSTANTS ####
const <- list(n.occasions = length(start:2021),
              offset = OFFSET, # difference in start times between population process (2008) and cmr data (1982)
              #r.j=apply(chick.marray,1,sum),
              #r.a=apply(adult.marray,1,sum),
              goodyear=goodyears$p.sel,
              #goodyear=goodyears$prop.seen,   ### if using a continuous effort correction
              juv.poss=phi.juv.possible$JuvSurv, ### sets the annual survival of juveniles to the mean if <70 were ringed
              
              ### count data
              n.sites.count=n.sites.count,
              n.years.count= n.years.count,
              prop.sites=mean.props,  ### need to calculate
              
              n.sites.fec=n.sites.fec,
              n.years.fec= n.years.fec
)

#### INITIAL VALUES ####

# TODO check that all the rounding is correct

p.juv.recruit.f.inits <- plogis(-4 + 1:30/2)

p.juv.recruit.inits <- matrix(rep(plogis(-4 + (1 / 2 * 1:30)), times = n.years.fec), nrow = 30, ncol = n.years.fec, byrow = FALSE) # changed to be divided by two

imean.p.ad <- c(runif(1, 0.05, 0.5), runif(1, 0.2, 1))
iann.fec <- rep(0.4, n.years.fec)

# TODO 
# this seems wrong
# update - yea, the indexing is wrong. there should only be 35 rows because those
# are the only years you would have a nonzero prob of detecting in year 36
ip.juv <- matrix(NA, nrow = length(start:2021)-1, ncol = length(start:2021)-1)
for (t in 1:(length(start:2021)-1)){
  for (j in 1:t){ ## zero by definition (these are never actually used)
    ip.juv[t,j] <- 0
  }
  for (j in (t+1):(length(start:2021)-1)){
    print(paste(j))
    ip.juv[t,j]  <- plogis(-4 + 1*(j - t)/2 + 0)
  }
}

ip.ad <- numeric(length(start:2020))
for (j in 1:length(start:2020)) {
  tmp <- c(0.4, 0.8)
  ip.ad[j] <- tmp[goodyears$p.sel[j]] 
}

iN.ad.surv <- rep(999, n.years.fec)
iN.breed.ready <- rep(999, n.years.fec)
inestlings <- rep(999, n.years.fec)
iN.recruits <- rep(NA, n.years.fec)
iNtot.breed <- rep(NA, n.years.fec)
iJUV <- rep(NA, n.years.fec)
iN.atsea <- rep(NA, n.years.fec)
iNtot <- rep(NA, n.years.fec)
IMinits <- array(NA, dim = c(n.years.fec, 30, 3))

IMinits[1,1,1] = max(c(rnorm(1, 263*2, 20)), 0) %>% round() # TODO change sd???
IMinits[1,2,1] = max(c(rnorm(1, 275*2, 20)), 0)%>% round()
IMinits[1,3,1] = max(c(rnorm(1, 264*2, 20)), 0)%>% round()
IMinits[1,4,1] = max(c(rnorm(1, 177*2, 20)), 0)%>% round()
IMinits[1,5,1] = max(c(rnorm(1, 290*2, 20)), 0)%>% round()
IMinits[1,6,1] = max(c(rnorm(1, 90*2, 20)), 0)%>% round()
IMinits[1,7,1] = max(c(rnorm(1, 158*2, 20)), 0)%>% round()

IMinits[1, 1, 2] <- 0
IMinits[1, 1, 3] <- IMinits[1,1,1] - IMinits[1,1,2]

for (age in 2:7) {
  IMinits[1, age, 2] <- rbinom(1, IMinits[1, age, 1], p.juv.recruit.f.inits[age])
  IMinits[1, age, 3] <- IMinits[1,age,1] - IMinits[1,age,2]
}

for (age in 8:30) {
  IMinits[1, age, 1] = rbinom(1, IMinits[1, age-1, 3], pow(0.9,(age-1)))
  IMinits[1, age, 2] <- rbinom(1, IMinits[1, age, 1], p.juv.recruit.f.inits[age])
  print(paste(age, IMinits[1,age,1] - IMinits[1,age,2], sep = " "))
  IMinits[1, age, 3] <- IMinits[1,age,1] - IMinits[1,age,2]
}

IMinits[1, , ]
IMinits[,, 1]

iN.recruits[1] <- sum(IMinits[1,1:30,2]) 
iNtot.breed[1] <- rnorm(1, 640*2,sd = 20)  %>% round()#change here

iJUV[1] <- rnorm(1, 232*2, sd = 20) %>% round()
iN.atsea[1] <- rnorm(1, 224*2,sd = 20)%>% round() #change here
iNtot[1]<-sum(IMinits[1,1:30,3]) + iNtot.breed[1]+iN.atsea[1]  

for (tt in 2:n.years.fec) {
  iN.ad.surv[tt] <- rbinom(1, iNtot.breed[tt-1]+iN.atsea[tt-1], 0.9)
  iN.breed.ready[tt] <- rbinom(1, iN.ad.surv[tt], ip.ad[1])
  
  IMinits[tt,1,1] <- rbinom(1, iJUV[tt-1], 0.75)
  IMinits[tt, 1, 2] <- 0
  IMinits[tt, 1, 3] <- IMinits[tt,1,1] - IMinits[tt,1,2]
  
  for (age in 2:30) {
    IMinits[tt, age, 1] <- rbinom(1, IMinits[tt-1, age-1, 3], 0.9)
    # TODO wtf is this
    #IMinits[tt, age, 2] <- min(IMinits[tt, age, 1], IMinits[tt, age-1, 3]) * p.juv.recruit.inits[age, tt] %>% round()
    IMinits[tt, age, 2] <- rbinom(1, IMinits[tt, age, 1], p.juv.recruit.inits[age, tt])
    IMinits[tt, age, 3] <- IMinits[tt,age,1] - IMinits[tt,age,2]
  }
  
  iN.recruits[tt] <- sum(IMinits[tt,1:30,2])
  iNtot.breed[tt]<- iN.breed.ready[tt]+iN.recruits[tt]
  inestlings[tt] <- (iann.fec[tt] * 0.5 * iNtot.breed[tt]) %>% round()
  iJUV[tt] <- rpois(1, inestlings[tt])
  
  iN.atsea[tt] <- iN.ad.surv[tt]-iN.breed.ready[tt]
  iNtot[tt]<-sum(IMinits[tt,1:30,3]) + iNtot.breed[tt]+iN.atsea[tt]
}

IMinits[,, 1]
IMinits[1, , ]

# checking if the inits are reasonable - this should be close to 0, on average
iNtot.breed * mean.props[1:13, ] - POP[1:13, ]
mean(iNtot.breed * mean.props[1:13, ] - POP[1:13, ])

# logic check
# there should not be negatives here
  for (year in 2:n.years.fec) {
    for(age in 2:30) {
    print(paste(year, age, IMinits[year-1, age-1, 3] - IMinits[year, age, 1] , sep = ' '))
  }
}

inits <- list(
  ann.fec = iann.fec,
  
  sigma.obs=matrix(rexp(n.sites.count*n.years.count, 0.1),ncol=n.years.count),
  
  p.juv.recruit.f = p.juv.recruit.f.inits,
  p.juv.recruit = p.juv.recruit.inits,
  
  mean.p.ad = imean.p.ad, 
  # 
  mu.p.juv = rnorm(2, -4, 0), # fixed because of the IM
  mu.p.ad = c(log(0.4/(1-0.4)), log(0.8/(1-0.8))),
   
  agebeta = rnorm(1, 1, 0), # fixed because of the IM
  
  p.juv = ip.juv,

  sigma.p = 0, # fixed because of the IM 
   
  mean.phi.ad = 0.9 , # fixed because of the IM
  mean.phi.juv = 0.75, # fixed because of the IM
  mu.juv = logit(0.75),
  mu.ad = logit(0.9),
   
  sigma.phi = 0, # fixed because of the IM
   
  phi.juv = rep(0.75, length(start:2021)-1),
  phi.ad = rep(0.9, length(start:2021)-1),
   
  eps.phi = rep(0, length(start:2021)-1),
   
  p.ad = ip.ad,
  # 
  eps.p = rep(0, length(start:2021)-1),
  # 
  # maybe better to not do this below
  IMnonround = IMinits[1, 1:7, 1],
  Ntot.breednonround = iNtot.breed[1],
  JUVnonround = iJUV[1],
  N.atseanonround = iN.atsea[1],
  IM = IMinits, # TODO
  N.ad.surv = iN.ad.surv,
  N.breed.ready = iN.breed.ready,
  nestlings = inestlings,
  N.recruits = iN.recruits,
  Ntot.breed = iNtot.breed,
  JUV = iJUV,
  N.atsea = iN.atsea,
  Ntot = iNtot
  
  #####
)


#### PARAMETERS TO MONITOR ####
params <- c(
  # SURVIVAL
  "mean.phi.juv", "mean.phi.ad", "sigma.phi", 
  "mu.p.juv", "mean.p.ad", "agebeta", "sigma.p",
  # FECUNDITY
  "ann.fec",
  # ABUNDANCE
  "Ntot", "Ntot.breed", "N.atsea", 
  "sigma.obs"
  # FUTURE
  #"fut.growth.rate", "fut.lambda", 
  #"Nobs.f"
)

#### MCMC SETTINGS ####
nb <- 1 #burn-in
ni <- 5000 + nb #total iterations
nt <- 1  #thin
nc <- 1  #chains
adaptInterval = 100
maxContractions = 1000

#### COMPILE CONFIGURE AND BUILD ####
Rmodel <- nimbleModel(code = code, constants = const, data = dat, 
                      check = FALSE, calculate = FALSE, inits = inits)
Rmodel$initializeInfo()
conf <- configureMCMC(Rmodel, monitors = params, thin = nt, 
                      control = list(maxContractions = maxContractions, 
                                     adaptInterval = adaptInterval)) # SLOWW
conf$printSamplers(type = "conjugate")
conf$printSamplers(type = "posterior") # check sampler defaults

# conf$removeSamplers("ann.fec")
# conf$addSampler(target = "ann.fec[1]", type="conjugate")
# conf$addSampler(target = "ann.fec[1:13]", type="AF_slice")
# conf$addSampler(target = "ann.fec[13]", type="posterior_predictive_branch")
# conf$printSamplers("ann.fec"

# TODO
# could block wrt to time
# using RW block samplers
# or AF slice samplers which are generally faster and mix better

# need to figure out why there is a conjugate sampler on one of the nodes

# Ntot breed N at see highly correlated
# juveniles and annual fecundity

Rmcmc <- buildMCMC(conf)  
Cmodel <- compileNimble(Rmodel, showCompilerOutput = FALSE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
library(beepr)
beep(sound = 1)

#### RUN MCMC ####
t.start <- Sys.time()
sink("somanyerrors.txt")
out <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = nc, inits = inits,
               setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE) 
sink()
t.end <- Sys.time()
(runTime <- t.end - t.start)

error.vec <- read_lines("somanyerrors.txt")
error.vec <- error.vec[str_detect(error.vec, "") & 
  !str_detect(error.vec, "lifted") &
  !str_detect(error.vec, "slice")] %>% unique() #%>% sort()
write_lines(error.vec, "somanyerrors.txt")
