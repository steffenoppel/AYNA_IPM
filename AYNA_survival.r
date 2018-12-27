##########################################################################
#
# ATLANTIC YELLOW-NOSED ALBATROSS SURVIVAL ANALYSIS 2000-2018
#
##########################################################################
# based on Kery and Schaub 2012, Chapter 7
# modified by Steffen oppel, December 2018

library(tidyverse)
library(jagsUI)
library(data.table)
library(nimble)



#########################################################################
# LOAD PRE-PREPARED DATA
#########################################################################

try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\AYNA_IPM"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Gough\\ANALYSIS\\AYNA_IPM"), silent=T)
AYNA<-fread("AYNA_simple_encounter_history_1982_2018.csv")
names(AYNA)
CH<-as.matrix(AYNA[,3:39], dimnames=F)
AYNA$AGE[is.na(AYNA$AGE)]<-1    ## set all NA as 'adult'

### check that there are contacts in every season
apply(CH,2,sum)



#########################################################################
# CREATE MATRIX OF AGE FOR EACH OCCASION AND INDIVIDUAL
#########################################################################

## this matrix will relate to the parameter estimates chosen in the model
## simple model only has 2 survival parameters:
## 1 - juvenile and immature survival (years 1-5)
## 2 - adult survival (birds >5 years old)

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x==1))
f <- apply(CH, 1, get.first)


## REMOVE BIRDS THAT ARE TOO YOUNG TO HAVE HAD A CHANCE TO RETURN
tooyoung<-ifelse(f>(dim(CH)[2]-5),ifelse(AYNA$AGE==0,1,0),0)
CH<-CH[tooyoung==0,]  ## removes individuals that were ringed as chicks <5 years before end of time series
ages<-AYNA$AGE[tooyoung==0]

## CREATE BLANK AGE MATRIX
AGEMAT<-CH
AGEMAT[,]<-2 ### set default to adult survival, and then insert the juvenile years
n.occ<-dim(AGEMAT)[2]

## LOOP OVER EACH BIRD RINGED AND SET PRE-CAPTURE DATA TO NA AND ADJUST AGE
for (l in 1:nrow(AGEMAT)){
  firstocc<-get.first(CH[l,])
  lastjuv<-firstocc+4
  lastjuv<-ifelse(lastjuv>n.occ,n.occ,lastjuv)
  young<-ages[l]
  if(firstocc>1){AGEMAT[l,1:(firstocc-1)]<-NA}  ## sets everything before first contact to NA
  if(young==0){AGEMAT[l,firstocc:lastjuv]<-1}  ## sets all juvenile years to 1
}




#########################################################################
# RE-ARRANGE DATA TO REMOVE CONTACTS BEFORE 2000 (and individuals with no contacts after 2000)
#########################################################################

## remove encounter occasions before 2000
rCH<-CH[,c(19:37)]	      ## reduce data set to exclude years before 2000
AGEMAT<-AGEMAT[,c(19:37)]	## reduce data set to exclude years before 2000
dim(rCH)
exclude<- apply(rCH,1,sum)
rCH<-rCH[exclude>0,]        ## removes individuals that were not observed in the last 19 years - could also set >1 to remove transients
AGEMAT<-AGEMAT[exclude>0,]  ## removes individuals that were not observed in the last 19 years
ages<-ages[exclude>0]
dim(rCH)


## PREPARE CONSTANTS
n.ind<-dim(rCH)[1]		## defines the number of individuals
n.years<-dim(rCH)[2]  ## defines the number of years
f <- apply(rCH, 1, get.first)








#########################################################################
# Specify basic CJS model with random time effects
#########################################################################

sink("AYNA_CJS_simple.jags")
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phi: survival probability, grouped by age category
    # p: recapture probability when breeding
    # -------------------------------------------------


    ## Priors and constraints
      for (i in 1:nind){
        for (t in f[i]:(n.occasions-1)){
          logit(phi[i,t]) <- mu[AGEMAT[i,t]] + surv.raneff[t]
          logit(p[i,t]) <- mean.p  + capt.raneff[t]
        } #t
      } #i
    
    
    ## AGE-SPECIFIC SURVIVAL 
    for (age in 1:2){
      beta[age] ~ dunif(0, 1)                         # Priors for age-specific survival
      mu[age] <- log(beta[age] / (1-beta[age]))       # Logit transformation
    }

    ## RANDOM TIME EFFECT ON SURVIVAL AND RECAPTURE 
    for (t in 1:(n.occasions-1)){
      surv.raneff[t] ~ dnorm(0, tau.surv)
      capt.raneff[t] ~ dnorm(0, tau.capt)
    }
    

    sigma.surv ~ dunif(0, 10)                     # Prior for standard deviation of survival
    tau.surv <- pow(sigma.surv, -2)
    
    sigma.capt ~ dunif(0, 10)                     # Prior for standard deviation of capture
    tau.capt <- pow(sigma.capt, -2)
    mean.p ~ dunif(0, 1)                          # Prior for mean recapture
    logit.p <- log(mean.p / (1-mean.p))           # Logit transformation
    


    # Likelihood 
    for (i in 1:nind){
      # Define latent state at first capture
      z[i,f[i]] <- 1
        for (t in (f[i]+1):n.occasions){
          # State process
          z[i,t] ~ dbern(mu1[i,t])
          mu1[i,t] <- phi * z[i,t-1]
          # Observation process
          y[i,t] ~ dbern(mu2[i,t])
          mu2[i,t] <- p * z[i,t]
        } #t
      } #i

    }
    ",fill = TRUE)
sink()





#########################################################################
# PREPARE DATA FOR MODEL
#########################################################################

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = n.years, nind = n.ind, AGEMAT=AGEMAT)

# Initial values 
inits <- function(){list(beta = runif(2, 0, 1),
                         mean.p = runif(1, 0, 1))}
 

# Parameters monitored
parameters <- c("beta", "mean.p")

# MCMC settings
ni <- 3000
nt <- 1
nb <- 1000
nc <- 4

# Call JAGS from R
AYNAsurv <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\AYNA_IPM\\AYNA_CJS_simple.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)




#########################################################################
# PRODUCE OUTPUT TABLE
#########################################################################

out<-as.data.frame(AYNAsurv$summary)
out$parameter<-row.names(AYNAsurv$summary)
out
write.table(out,"AYNA_Gough_Survival_estimates.csv", sep=",", row.names=F)



#########################################################################
# PRODUCE OUTPUT GRAPH
#########################################################################


out[1:11,] %>% select(c(1,5,2,3,7)) %>%
  setNames(c('Mean', 'Median','SD','lcl', 'ucl')) %>%
  mutate(Year=colnames(rCH)[1:11]) %>%
  
  ggplot(aes(y=Median, x=Year)) + geom_point(size=2.5)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  ylab("Annual adult survival probability") +
  scale_y_continuous(breaks=seq(0.5,1,0.1), limits=c(0.5,1))+
  #scale_x_continuous(breaks=seq(2006,2017,1))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()



