##########################################################################
#
# ATLANTIC YELLOW-NOSED ALBATROSS INTEGRATED POPULATION MODEL 2008-2050
#
##########################################################################
# based on https://github.com/steffenoppel/TRAL_IPM
# modified for AYNA on 9 Nov 2021


library(tidyverse)
library(lubridate)
library(data.table)
library(jagsUI)
library(runjags)   ## added by Beth in July 2021 because jagsUI would not converge
filter<-dplyr::filter
select<-dplyr::select



#########################################################################
# LOAD PRE-PREPARED DATA ON COUNTS AND BREEDING SUCCESS
#########################################################################
### see 'IPM_DATA_PREPARATION_AYNA.R' for details on how data are aggregated

### NOTE THAT SORT ORDER OF GONYDALE AND GP VALLEY HAS SHIFTED ON 15 Jan 2021 (due to switch to if_else on R4.0.2)

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

### SCALE NUMBER OF HOOKS

longline <- longline %>% mutate(n_hooks = scale(n_hooks)) 
ave.since.2010 <- longline %>% filter(Year > 2009 & Year < 2020) %>% select(2) %>% unlist() %>% mean() 
longline <- longline %>% 
  mutate(n_hooks = if_else(Year == 2020 | Year == 2021, ave.since.2010, n_hooks))
longline


#########################################################################
# SPECIFY MODEL IN JAGS
#########################################################################
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM")
sink("AYNA_survivalOnly.jags")
cat("

model {
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
    
#-------------------------------------------------  
# 1. PRIORS FOR ALL DATA SETS
#-------------------------------------------------
    
    
    # -------------------------------------------------        
    # 1.3. Priors and constraints FOR SURVIVAL
    # -------------------------------------------------
    
    ### RECAPTURE PROBABILITY
    mean.p.ad[1] ~ dunif(0.0, 0.5)	           # Prior for mean adult recapture - should be higher than 5% but less than 50%
    mean.p.ad[2] ~ dunif(0.2, 1)	           # Prior for mean adult recapture - should be higher than 20%

    for (gy in 1:2){    ## for good and poor monitoring years
      # TODO - could put more informative priors here
      # but also note that the uniform prior on the logit scale is informative
      mean.p.juv[gy] ~ dunif(0, 1)	         # Prior for mean juvenile recapture - should be higher than 20% if they survive!
      mu.p.juv[gy] <- log(mean.p.juv[gy] / (1-mean.p.juv[gy])) # Logit transformation
      mu.p.ad[gy] <- log(mean.p.ad[gy] / (1-mean.p.ad[gy])) # Logit transformation
    }
    agebeta ~ dunif(0,1)    # Prior for shape of increase in juvenile recapture probability with age
    beta.ICCAT.ll.e ~ dnorm(0, 1)  # TODO - change precison?
    beta.ICCAT.ll.mit ~ dnorm(0, 1)  # TODO - change precison?
    beta.Nam.ll.mit ~ dnorm(0, 1) # TODO - change precison?
    beta.SA.ll.mit ~ dnorm(0, 1) # TODO - change precison?
    beta.Uru.ll.mit ~ dnorm(0, 1) # TODO - change precison?

    ## RANDOM TIME EFFECT ON RESIGHTING PROBABILITY OF JUVENILES
    for (t in 1:(n.occasions-1)){
      for (j in 1:t){ ## zero by definition (these are never actually used)
        p.juv[t,j] <- 0
      }
      for (j in (t+1):(n.occasions-1)){
        logit(p.juv[t,j])  <- mu.p.juv[goodyear[j]] + agebeta*(j - t) + eps.p[j]
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
    # TODO - add additional covariates wrt to fishing effory and bycatch mitigation
    for (j in 1:(n.occasions-1)){
      logit(phi.juv[j]) <- mu.juv + eps.phi[j]*juv.poss[j] + beta.ICCAT.ll.e*ICCAT.ll.e[j] + beta.ICCAT.ll.mit*ICCAT.ll.mit[j] + beta.Nam.ll.mit*Nam.ll.mit[j] + beta.SA.ll.mit*SA.ll.mit[j] + beta.Uru.ll.mit*Uru.ll.mit[j]
      logit(phi.ad[j]) <- mu.ad + eps.phi[j] + beta.ICCAT.ll.e*ICCAT.ll.e[j] + beta.ICCAT.ll.mit*ICCAT.ll.mit[j] + beta.Nam.ll.mit*Nam.ll.mit[j] + beta.SA.ll.mit*SA.ll.mit[j] + beta.Uru.ll.mit*Uru.ll.mit[j]
      eps.phi[j] ~ dnorm(0, tau.phi) 
      logit(p.ad[j])  <- mu.p.ad[goodyear[j]] + eps.p[j]    #### CAT HORSWILL SUGGESTED TO HAVE A CONTINUOUS EFFORT CORRECTION: mu.p.ad + beta.p.eff*goodyear[j] + eps.p[j]
      eps.p[j] ~ dnorm(0, tau.p)
    }
    
    
    
#-------------------------------------------------  
# 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
#-------------------------------------------------
    
    
    # -------------------------------------------------        
    # 2.4. Likelihood for adult and juvenile survival from CMR
    # -------------------------------------------------
    
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
      marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], r.j[t])
      marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], r.a[t])
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
    
}  ## end model loop

",fill = TRUE)
sink()





#########################################################################
# PREPARE DATA FOR MODEL
#########################################################################

# Bundle data
jags.data <- list(marr.j = chick.marray,
                  marr.a = adult.marray,
                  n.occasions = dim(chick.marray)[2],
                  r.j=apply(chick.marray,1,sum),
                  r.a=apply(adult.marray,1,sum),
                  goodyear=goodyears$p.sel,
                  #goodyear=goodyears$prop.seen,   ### if using a continuous effort correction
                  juv.poss=phi.juv.possible$JuvSurv, ### sets the annual survival of juveniles to the mean if <70 were ringed
                  
                  ### count data
                  #n.sites.count=n.sites.count,
                  #n.years.count= n.years.count,
                  #prop.sites=mean.props,  ### need to calculate
                  #y.count=POP,    ### use log(R) here if using the logscale model
                  
                  ### breeding success data
                  #J=PROD.DAT$J,
                  #R=PROD.DAT$R,
                  #n.sites.fec=n.sites.fec,
                  #n.years.fec= n.years.fec,
                  
                  ### longline effort data
                  ICCAT.ll.e = longline$n_hooks %>% as.numeric(), 
                  ICCAT.ll.mit = longline$mit.ICCAT %>% as.numeric(), 
                  Nam.ll.mit = longline$mit.NAM %>% as.numeric(), 
                  SA.ll.mit = longline$mit.RSA %>% as.numeric(),
                  Uru.ll.mit = longline$mit.URU %>% as.numeric()
                  
                  # ### FUTURE PROJECTION
                  #FUT.YEAR=30,  ### for different scenarios future starts at 1
                  #n.scenarios=1,
                  #fut.surv.change=as.matrix(fut.surv.change[,2]),  ## future survival rate change - matrix that adjusts gradual decrease in survival
                  #fut.fec.change=c(1)     ## future fecundity change - vector with one element for each scenario
)


# Initial values 
inits <- function(){list(mean.phi.ad = runif(1, 0.7, 0.97),
                         mean.phi.juv = runif(1, 0.5, 0.9),
                         mean.p.ad = c(runif(1, 0.05, 0.5), runif(1, 0.2, 1)),
                         mean.p.juv = runif(2, 0, 1),
                         beta.ICCAT.ll.e = rnorm(1, 0, 1),
                         beta.ICCAT.ll.mit = rnorm(1, 0, 1),
                         beta.Nam.ll.mit = rnorm(1, 0, 1), 
                         beta.SA.ll.mit = rnorm(1, 0, 1), 
                         beta.Uru.ll.mit = rnorm(1, 0, 1)
                        
                         #Ntot.breed= c(runif(1, 4950, 5050),rep(NA,n.years.fec-1)), # TODO change this
                         #JUV= c(rnorm(1, 246, 0.1),rep(NA,n.years.fec-1)), # TODO change this
                         #N.atsea= c(rnorm(1, 530, 0.1),rep(NA,n.years.fec-1)), # TODO change this
                         # IM[,1,1]= c(rnorm(1, 324, 0.1),rep(NA,n.years-1)),
                         # IM[,2,1]= c(rnorm(1, 257, 0.1),rep(NA,n.years-1)),
                         # IM[,3,1]= c(rnorm(1, 462, 0.1),rep(NA,n.years-1)),
                         # IM[,4,1]= c(rnorm(1, 207, 0.1),rep(NA,n.years-1)),
                         # IM[,5,1]= c(rnorm(1, 700, 0.1),rep(NA,n.years-1)),
                         # IM[,6,1]= c(runif(1, 150, 300),rep(NA,n.years-1)),
                         #sigma.obs=matrix(runif(n.sites.count*n.years.count,1,20),ncol=n.years.count))
                         )}



# Parameters monitored
parameters <- c("mean.phi.ad","mean.phi.juv", "mean.p.ad", 'mean.p.juv', "phi.ad", "phi.juv", 
                "beta.ICCAT.ll.e", "beta.ICCAT.ll.mit", "beta.Nam.ll.mit", "beta.SA.ll.mit", "beta.Uru.ll.mit")

# MCMC settings
nt <- 1#0
nb <- 25000
nad <- 2000
nc <- 3
ns <- 20000#0 #longest

# run the model in run jags
start.time <- Sys.time()
AYNAipm <- run.jags(data=jags.data, inits=inits, parameters, 
                    model="AYNA_survivalOnly.jags",
                    n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                    method = "rjparallel") 
end.time <- Sys.time()
(run.time <- end.time - start.time)


#########################################################################
# SAVE OUTPUT - RESULT PROCESSING in AYNA_IPM_result_summaries.r
#########################################################################
### DO NOT UPLOAD THIS TO GITHUB - IT WILL CORRUPT THE REPOSITORY

## updated script for 'runjags' output
summary_AYNAipm <- summary(AYNAipm)
library(coda)
plot(AYNAipm)
gelman.diag(AYNAipm, multivariate = FALSE, autoburnin = TRUE)
summary(AYNAipm)

# monitor annual survival values and plot against whether it's a good or bad year
goodyears$p.sel
goodyears$prop.seen
library(stringr)
survival_posteriors <- AYNAipm$mcmc[, str_detect(colnames(AYNAipm$mcmc[[1]]),"phi.ad\\[")][[1]]
plot(goodyears$prop.seen[1:43],apply(survival_posteriors, 2, median))
boxplot(apply(survival_posteriors, 2, median) ~ goodyears$p.sel[1:43])

summary_AYNAipm_df <- as.data.frame(summary_AYNAipm)
View(summary_AYNAipm_df)
head(summary_AYNAipm_df)
min(summary_AYNAipm_df$SSeff) #Ntot[1]
max(summary_AYNAipm_df$psrf) #Ntot[1]

addsummary_AYNAipm <- add.summary(AYNAipm,plots = runjags.getOption("predraw.plots"))
addsummary_AYNAipm #18 min

plot(addsummary_AYNAipm, layout=c(2,2))

predictions <- data.frame(summary(addsummary_AYNAipm),
                          parameter = row.names(summary(addsummary_AYNAipm)))
head(predictions)
row.names(predictions) <- 1:nrow(predictions)

predictions <- predictions[1:218,]   ### 200 cuts off ann.fec
#predictions[1:5,]

predictions$Mode <- NULL
np <- names(predictions) 
names(predictions) <- c("lcl",np[2],"ucl",np[4:9],"Rhat",np[11])

max(predictions$Rhat)



setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM")
save.image("AYNA_IPM_output_FINAL.RData")