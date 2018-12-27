#####################################################################################################
##### ATLANTIC YELLOW-NOSED ALBATROSS POPULATION TREND 2000 - 2018 #################
#####################################################################################################
# code from SHPL recovery paper (Leo et al. 2018)
# data prepared using script IPM_DATA_PREPARATION.R
# first modified 27 Dec 2018
# removed multi-site to run model in NIMBLE


################################## SET WORKING DIRECTORY  ######################################

library(jagsUI)
library(nimble)
library(tidyverse)
library(data.table)
library(ggplot2)

try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\AYNA_IPM"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Gough\\ANALYSIS\\AYNA_IPM"), silent=T)




################################## LOAD DATA  ######################################

AYNA<-fread("AYNA_pop_counts_1982_2018.csv")
AYNA<-subset(AYNA,Year>1999)	## reduce data set to remove NA in 4 years
n.years<-dim(AYNA)[1]		## defines the number of years
n.sites<-dim(AYNA)[2]-1 ## defines the number of study areas
str(AYNA)
names(AYNA)



################################## SIMPLE STATE SPACE MODEL  ######################################

AYNA_trend<-nimbleCode({

# Priors and constraints

for (s in 1:n.sites){			### start loop over every study area
  N.est[1,s] ~ dunif(0,200)   ## draw random value from a uniform distribution between 0 and 200 for initial population size
  mean.lambda[s] ~ dunif(0,10)	#Prior for mean growth rate
  sigma.proc[s] ~ dunif(0,10)	#Prior for SD of state process (annual variation in pop size)
  sigma2.proc[s]<-pow(sigma.proc[s],2)
  tau.proc[s]<-pow(sigma.proc[s],-2)
  sigma.obs[s] ~ dunif(0,100)	#Prior for SD of observation process (variation in detectability)
  sigma2.obs[s]<-pow(sigma.obs[s],2)
  tau.obs[s]<-pow(sigma.obs[s],-2)
  }


##### Likelihood function

for (s in 1:n.sites){			### start loop over every study area
## State process for entire time series

  for (t in 1:(T-1)){
    lambda[t,s] ~ dnorm(mean.lambda[s], tau.proc[s])								# Distribution for random error of growth rate
    N.est[t+1,s]<-N.est[t,s]*lambda[t,s]										# Linear predictor (population size based on past pop size and change rate)
  }														# run this loop over nyears


## Observation process

for (t in 1:T){
  y[t,s] ~ dnorm(N.est[t,s], tau.obs[s])								# Distribution for random error in observed numbers (counts)
  }														# run this loop over t= nyears
}		## end site loop

## Derived parameters
for (t in 1:T){
  pop.size[t]<-sum(N.est[t,1:n.sites])
  }

  mlam <- mean(lambda[1:(T-1),1:n.sites])  				# Arithmetic mean for whole time series
site.lam<-mean(mean.lambda[1:n.sites])

})													# close the model loop






#################################################################
#	PREPARE DATA AND SET UP MODEL RUN
#################################################################

## Initial values
N.init=matrix(NA, nrow=n.years,ncol=n.sites)
N.init[1,]<-as.matrix(AYNA[1,2:12])
inits<- function() {list(sigma.proc=runif(n.sites,0,5),
    mean.lambda=runif(n.sites,0.1,2),
    sigma.obs=runif(n.sites,0,10),
    N.est=N.init)}

AYNAinits<-inits()

## Parameters to be estimated ('monitored') by NIMBLE
params<-c("pop.size", "site.lam")

## MCMC settings
ni<-60000  	## number of iterations (draws per chain)
nt<-1		## thinning rate
nb<-35000	## length of burn-in (number of iterations discarded at the start of each chain until convergence is reached)
nc<-4		## number of chains


## call JAGS from R   ###
#AYNAtrend<- jags(data=bugs.data, inits=inits, parameters.to.save=params, model.file="A:\\RSPB\\UKOT\\StHelena\\Science\\Birds\\Census_data\\SHPL_trend.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.cores=nc, parallel=T)
#AYNAtrend<- jags(data=bugs.data, inits=inits, parameters.to.save=params, model.file="C:\\STEFFEN\\RSPB\\UKOT\\StHelena\\Science\\Birds\\Census_data\\SHPL_trend.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.cores=nc, parallel=T)

AYNAtrend <- nimbleMCMC(AYNA_trend, 
                     data = list(y=as.matrix(AYNA[,2:12])),			### could also use win.data as constants here
                     constants = list(T=n.years, n.sites=n.sites),
                     inits = AYNAinits,
                     monitors=params,
                     nchains=nc, thin=nt, niter=ni, nburnin=nb,
                     samplesAsCodaMCMC = TRUE, summary=TRUE)

str(AYNAtrend)




#################################################################
#	SUMMARISE OUTPUT AND PLOT POPULATION TREND FOR GLOBAL POPULATION
#################################################################


## COMPILE DATA FROM NIMBLE OUTPUT
OUT<-as.data.frame(AYNAtrend$summary$all.chains)
write.table(OUT,"AYNA_trend_estimates2018_nimble.csv", sep=",")



################ SET UP PLOT WINDOW AND PLOT THE GRAPH FOR THE WHOLE ISLAND ####################
library(data.table)
library(tidyverse)


pdf("AYNA_pop_trend_Gough_2000_2018.pdf", width=11, height=8)
OUT %>% filter(Mean>10) %>%
  mutate(Year=AYNA$Year) %>%
  setNames(c('Mean', 'Median','SD','lcl', 'ucl','Year')) %>%
  
ggplot(aes(y=Median, x=Year)) + geom_point(size=2.5)+ geom_line()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  #geom_point(data=data.frame(N=N,year=AYNA$Year), aes(x=year,y=N), colour='red',pch=4,size=2.5)+ 
  ylab("Number of AYNA pairs") +
  scale_y_continuous(breaks=seq(0,1000,100), limits=c(0,1000))+
  scale_x_continuous(breaks=seq(2000,2018,2))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
dev.off()



