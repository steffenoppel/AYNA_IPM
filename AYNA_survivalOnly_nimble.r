#### COMMENTS SECTION ####


#### LOAD LIBRARIES #####
library(nimble)
library(here)
library(coda)
library(tidyverse)
library(lubridate)
library(data.table)
#library(jagsUI)
#library(runjags)   ## added by Beth in July 2021 because jagsUI would not converge
filter<-dplyr::filter
select<-dplyr::select





#########################################################################
# LOAD PRE-PREPARED DATA ON COUNTS AND BREEDING SUCCESS
#########################################################################
### see 'IPM_DATA_PREPARATION_AYNA.R' for details on how data are aggregated

## LOAD PREPARED M-ARRAY FOR SURVIVAL ESTIMATION
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM")
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


### PLOT TO SPOT ANY OUTLIERS OF BCOUNTS
#ggplot(AYNA.pop, aes(x=Year,y=tot)) +geom_point(size=2, color='darkred')+geom_smooth(method='lm') 


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
# survival analysis runs from 1978-2021, but recapture index refers to columns, which represent year 1979-2021 plus the ones never recaptured (last column)
# very difficult
names(AYNA_CHICK)
POPSIZE$Year

OFFSET<-min(which(!is.na(match(as.numeric(substr(names(AYNA_CHICK)[2:44],1,4)),POPSIZE$Year))))
substr(names(AYNA_CHICK),1,4)[OFFSET+1]

### SCALE NUMBER OF HOOKS

longline <- longline %>% mutate(n_hooks = scale(n_hooks)) 
ave.since.2010 <- longline %>% filter(Year > 2010) %>% select(2) 
ave.since.2010 <- mean(ave.since.2010$n_hooks, na.rm = T)
longline <- longline %>% 
  mutate(n_hooks = if_else(Year %in% c(2020, 2021), ave.since.2010, n_hooks))
longline

#########################################################################
# SPECIFY FUTURE DECREASE IN SURVIVAL
#########################################################################

dec.surv=0.9  ## we assume that adult survival will decrease by 10%
lag.time=10    ## the decrease will take 10 years to materialise
PROJECTION.years<-seq(1,30,1)  ## we specify the relative survival decrease for all 30 years in the projection

fut.surv.change<- expand.grid(PROJECTION.years,dec.surv,lag.time) %>%
  rename(Year=Var1,SURV3=Var2,LAG=Var3) %>%
  mutate(ann.offset=(SURV3-1)/LAG) %>%
  mutate(SURV3=ifelse(Year<LAG,1+(Year*ann.offset),SURV3)) %>%
  mutate(SURV1=1,SURV2=1) %>%
  select(Year, SURV1,SURV2,SURV3)

#### MODEL CODE ####
code <- nimbleCode({
  #-------------------------------------------------
  # integrated population model for the Gough AYNA population
  # - age structured model with 30 age classes 
  # - adult survival based on CMR ringing data
  # - pre breeding census, female-based assuming equal sex ratio & survival
  # - productivity based on all areas incu and chick counts
  # - linked population process with SUM OF count data
  # - v4 includes 3 scenarios of future projection: no change, improved fecundity, reduced adult survival
  # - marray_v1 uses marray for survival estimation to speed up computation time
  # -------------------------------------------------
  
  # changes from jags syntax
  # - trunction notation
  # - everything has to be indexed within sums over age classes
  # - easier to use sd instead of tau for normally distributed stuff
  
  # TODO AEB
  # tidy this up
  # remove some of the hardcoding
  # add covariates
  
  #-------------------------------------------------  
  # 1. PRIORS FOR ALL DATA SETS
  #-------------------------------------------------
  
  
  # -------------------------------------------------        
  # 1.3. Priors and constraints FOR SURVIVAL
  # -------------------------------------------------
  
  ### RECAPTURE PROBABILITY
  for (gy in 1:2){    ## for good and poor monitoring years
    # TODO - could put more informative priors here
    #mean.p.juv[gy] ~ dunif(0, 0.01)	         # Prior for mean juvenile recapture - should be higher than 20% if they survive!
    mean.p.ad[gy] ~ dunif(0.2, 1)	           # Prior for mean adult recapture - should be higher than 20%
    #mu.p.juv[gy] <- log(mean.p.juv[gy] / (1-mean.p.juv[gy])) # Logit transformation
    mu.p.ad[gy] <- log(mean.p.ad[gy] / (1-mean.p.ad[gy])) # Logit transformation
  }
  agebeta ~ dunif(0,1)    # Prior for shape of increase in juvenile recapture probability with age
  beta.fe ~ dnorm(0, sd = 1)  # TODO - change precison?
  
  ## RANDOM TIME EFFECT ON RESIGHTING PROBABILITY OF JUVENILES
  for (t in 1:(n.occasions-1)){
    for (j in 1:t){ ## zero by definition (these are never actually used)
      p.juv[t,j] <- 0
    }
    for (j in (t+1):(n.occasions-1)){
      #logit(p.juv[t,j])  <- mu.p.juv[goodyear[j]] + agebeta*(j - t) + eps.p[j]
      logit(p.juv[t,j])  <-  agebeta*(j - t) + eps.p[j]
      
    }
  }
  
  ## PRIORS FOR RANDOM EFFECTS
  sigma.p ~ dexp(1)                # Prior for standard deviation
  tau.p <- pow(sigma.p, -2)
  
  
  ### SURVIVAL PROBABILITY
  mean.phi.juv ~ dbeta(75.7,24.3)             # Prior for mean juvenile survival first year 0.757, second year 0.973 in Laysan albatross
  mean.phi.ad ~ dbeta(91,9)              # Prior for mean adult survival - should be higher than 70%
  mu.juv <- log(mean.phi.juv / (1-mean.phi.juv)) # Logit transformation
  mu.ad <- log(mean.phi.ad / (1-mean.phi.ad)) # Logit transformation
  
  ## PRIORS FOR RANDOM EFFECTS
  sigma.phi ~ dexp(1)                # Prior for standard deviation
  tau.phi <- pow(sigma.phi, -2)
  
  ## RANDOM TIME EFFECT ON SURVIVAL AND ADULT RECAPTURE
  
  for (j in 1:(n.occasions-1)){
    logit(phi.juv[j]) <- mu.juv + eps.phi[j]*juv.poss[j] #+ beta.fe*longline[j] 
    logit(phi.ad[j]) <- mu.ad + eps.phi[j] #+ beta.fe*longline[j]
    eps.phi[j] ~ dnorm(0, sd = sigma.phi) 
    logit(p.ad[j])  <- mu.p.ad[goodyear[j]] + eps.p[j]    #### CAT HORSWILL SUGGESTED TO HAVE A CONTINUOUS EFFORT CORRECTION: mu.p.ad + beta.p.eff*goodyear[j] + eps.p[j]
    eps.p[j] ~ dnorm(0, sd = sigma.p)
  }
  
  
  
  #-------------------------------------------------  
  # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
  #-------------------------------------------------
  

  # -------------------------------------------------        
  # 2.4. Likelihood for adult and juvenile survival from CMR
  # -------------------------------------------------
  
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
  
})

#### DATA ####
dat <- list(marr.j = chick.marray,
            marr.a = adult.marray#,
            
            #y.count=POP,    ### use log(R) here if using the logscale model
            
            ### breeding success data
            #J=PROD.DAT$J,
            #R=PROD.DAT$R
            
            ## future fecundity change - vector with one element for each scenario
) 

#### CONSTANTS ####
const <- list(n.occasions = length(start:2021),
              #offset = OFFSET, # difference in start times between population process (2008) and cmr data (1982)
              r.j=apply(chick.marray,1,sum),
              r.a=apply(adult.marray,1,sum),
              goodyear=goodyears$p.sel,
              #goodyear=goodyears$prop.seen,   ### if using a continuous effort correction
              juv.poss=phi.juv.possible$JuvSurv#, ### sets the annual survival of juveniles to the mean if <70 were ringed
              
              ### count data
              # n.sites.count=n.sites.count,
              # n.years.count= n.years.count,
              # prop.sites=mean.props,  ### need to calculate
              # 
              # n.sites.fec=n.sites.fec,
              # n.years.fec= n.years.fec,
              # 
              # ### longline effort data
              # longline=longline$n_hooks %>% as.numeric(),
              # 
              # # ### FUTURE PROJECTION
              # FUT.YEAR=30,  ### for different scenarios future starts at 1
              # n.scenarios=1,
              # fut.surv.change=as.matrix(fut.surv.change[,2]),  ## future survival rate change - matrix that adjusts gradual decrease in survival
              # fut.fec.change=rep(1, 30)  
)

#### INITIAL VALUES ####
inits <- list(sigma.phi = rexp(1, 1),
              mean.phi.ad = rbeta(1, 91,9) ,
              mean.p.ad = runif(2, 0.2, 1),
              sigma.p = rexp(1, 1), 
              #Ntot.breed= c(rnorm(1, 640, sd = 20),rep(NA,n.years.fec-1)), # TODO change this
              #JUV= c(rnorm(1, 232, sd = 20),rep(NA,n.years.fec-1)), # TODO change this
              #N.atsea= c(rnorm(1, 224, sd = 20),rep(NA,n.years.fec-1)), # TODO change this
              beta.fe = rnorm(1, 0, 1),
              agebeta = runif(1, 0, 1)#,    
              # maybe better to not do this below
              # IM[1:n.years.count,1,1] = c(rnorm(1, 263, 20),rep(NA,n.years.count-1)), # TODO change sd???
              # IM[1:n.years.count,2,1] = c(rnorm(1, 275, 20),rep(NA,n.years.count-1)),
              # IM[1:n.years.count,3,1] = c(rnorm(1, 264, 20),rep(NA,n.years.count-1)),
              # IM[1:n.years.count,4,1] = c(rnorm(1, 177, 20),rep(NA,n.years.count-1)),
              # IM[1:n.years.count,5,1] = c(rnorm(1, 290, 20),rep(NA,n.years.count-1)),
              # IM[1:n.years.count,6,1] = c(rnorm(1, 90, 20),rep(NA,n.years.count-1)),
              # IM[1:n.years.count,7,1] = c(rnorm(1, 158, 20),rep(NA,n.years.count-1)),
              #####
              #ann.fec = rbeta(n.years.fec, 32,68),
              #sigma.obs=matrix(rexp(n.sites.count*n.years.count, 0.1),ncol=n.years.count)
)

#### PARAMETERS TO MONITOR ####
# TODO check that this has everything 
params <- c("mean.phi.ad","mean.phi.juv",#"mean.fec",
            #"pop.growth.rate" ,#"fut.growth.rate", 
            "agebeta",#"Ntot","Ntot.breed",
            "phi.ad","phi.juv",    
            #"ann.fec", "sigma.obs",
            "mean.p.ad"
            )

#### MCMC SETTINGS ####
nb <- 1 #burn-in
ni <- 50 + nb #total iterations
nt <- 1  #thin
nc <- 3  #chains
adaptInterval = 100
maxContractions = 1000

#### COMPILE CONFIGURE AND BUILD ####
Rmodel <- nimbleModel(code = code, constants = const, data = dat, 
                      check = TRUE, calculate = TRUE, inits = inits)
Rmodel$simulate()
Rmodel$calculate()
write_lines(Rmodel$getCode(), "AYNAipm_nimble.txt")


conf <- configureMCMC(Rmodel, monitors = params, thin = nt, 
                      control = list(maxContractions = maxContractions, 
                                     adaptInterval = adaptInterval)) # SLOWW
# lots of initial model checking you can do by exploring conf
# if you wanted to change samplers this is where you would do that

conf$printSamplers(type = "conjugate")
conf$printSamplers(type = "posterior") # check sampler defaults
conf

# TODO
# could block wrt to time
# using RW block samplers
# or AF slice samplers which are generally faster and mix better

# need to figure out why there is a conjugate sampler on one of the nodes

# Ntot breed N at see highly correlated
# juveniles and annual fecundity

## example
# conf$removeSamplers("ann.fec")
# conf$addSampler(target = "ann.fec[1]", type="conjugate")
# conf$addSampler(target = "ann.fec[1:13]", type="AF_slice")
# conf$addSampler(target = "ann.fec[13]", type="posterior_predictive_branch")
# conf$printSamplers("ann.fec")

Rmcmc <- buildMCMC(conf)  
Cmodel <- compileNimble(Rmodel, showCompilerOutput = FALSE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
# library(beepr)
# beep(sound = 8)

#### RUN MCMC ####
t.start <- Sys.time()
sink("somanyerrors.txt")
out <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = nc, inits = inits,
               setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE) 
sink()
t.end <- Sys.time()
(runTime <- t.end - t.start)

error.vec <- read_lines("somanyerrors.txt")
error.vec <- error.vec[!(str_detect(error.vec, "initializing") & 
                           #str_detect(error.vec, "IM\\[") & 
                           str_detect(error.vec, "Inf") |
                           str_detect(error.vec, "lifted") | str_detect(error.vec, "slice") )
] %>% unique() %>% sort()
write_lines(error.vec, "somanyerrors.txt")

#### MAKE BEAUTIFUL PLOTS AND STUFF ####
summary(out)
traceplot(out)