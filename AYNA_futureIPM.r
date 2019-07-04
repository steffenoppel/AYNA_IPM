##########################################################################
#
# ATLANTIC YELLOW-NOSED ALBATROSS INTEGRATED POPULATION MODEL 2000-2018
#
##########################################################################
# based on Kery and Schaub 2012, Chapter 11
# modified by Steffen oppel, December 2018
# code for survival and trend SSM in separate files
# this IPM is much simpler than previous attempts by Horswill, Converse etc.

# modified on 2 January 2019 to include projection into future
# fixed errors based on Martyn Plummers comments on 3 Jan 2019
# persisting 'invalid parent' error on ann.recruits - presumably because pool of immature birds drops to 0?
# 4 Jan 2019 - reduced recruitment rate prior from 0.05,0.95 to 0.15,0.75
# completed first future run (v3), but pop growth rate in future is <<1, so exploring whether changing parameters will stabilise pop trajectory
# 6 January 2019 - v4 includes mean demographic rates rather than single-year realisations - this always led to invalid parent error
# 7 January 2019 - reverted back to v3 but fixed the future growth rate


# 20 May 2019 - first attempt to include fishing effort data (provided by Nina DaRocha)
# 28 May 2019 - corrected ICCAT data after realising they provide lats/longs in quadrants

# revised 10 June 2019 to include the new longline effort data sent by Ana Carneiro - removed 17 June because data from 2017 are questionable
# revised IPM to base future projection on past 3 years rather than total series - to simulate projection with improved mitigation

# revised 1 July 2019 to include Namibian demersal longline data provided by Nina daRocha (ATF)
# v3 of IPM includes adjustments to project future pop growth on recent surv and average fecundity

# revised 2 July after exploring if survival of potential recruiters should be imm or ad survival - future trajectory changes if using imm survival (negative trend) vs. ad survival (stable)
# decided to use adult survival for N6 upwards, as that matches the age-matrix of the CJS model

# revised 2 July : build in 4 scenarios - mouse eradication, bycatch reduction, both, neither
# mouse eradication: fec goes from 0.56 - 0.69
# bycatch reduction: use surv from last 4 years rather than mean across earlier years

# output processing outsourced to AYNA_IPM_result_summaries.r



library(tidyverse)
library(jagsUI)
library(data.table)
#library(nimble)
filter<-dplyr::filter
select<-dplyr::select


#########################################################################
# LOAD FISHERY DATA PROVIDED BY JOEL RICE
#########################################################################
try(setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch"), silent=T)
nhooks<-fread("SH_longline_effort_30MAY_2019.csv")
head(nhooks)
nhooksSummary<-nhooks %>% 
  filter(lat5>(-40.1)) %>% filter(lat5<(-20.1)) %>%
  filter(lon5<20.1) %>% filter(lon5>(-15.1)) %>%
  filter(yy>1999) %>%
  group_by(yy) %>%
  summarise(N=sum(hooks))

## scale
longlineJR<- (nhooksSummary$N-mean(nhooksSummary$N))/sd(nhooksSummary$N)


# #########################################################################
# # LOAD FISHERY DATA FROM ICCAT (n hooks 2000 - 2017)
# #########################################################################
# ### more detailed data from Namibia do not cover the same time horizon
# 
# ## this CSV file is based on a query that uses QuadID==2 for only SE Atlantic!
# ## see https://www.iccat.int/Data/t2ce-ENG.pdf
# 
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\BycatchData"), silent=T)
nhooks<-fread("N_Hooks2000_2017.csv")
head(nhooks)


#### FIXED WITH QuadID==2
### format coordinates - these are unbelievably not specified as N or S but 'QuadID' indicates hemisphere
### extract data from 20-40 S and 20E to 10W

nhooksSummary<-nhooks %>% #mutate(Lat=ifelse(grepl('n',DSetTypeID)==T,Lat,Lat*-1)) %>%
  #mutate(Lon=ifelse(grepl('w',DSetTypeID)==T,Lon*-1,Lon)) %>%
  filter(Lat<(40.1)) %>% filter(Lat>(19.9)) %>%
  filter(Lon<20.1) %>%
  group_by(YearC) %>%
  summarise(N=sum(SumOfEff1))

## scale
longlineICCAT<- (nhooksSummary$N-mean(nhooksSummary$N))/sd(nhooksSummary$N)
# 
# 
# ### load long-line fishing effort from Namibia
# ### this only goes back to 2009, but in some years is greater than what ICCAT report for SE Atlantic!!
# 
# NamLLeff<-fread("Longline_effort_Namibia.csv")
# head(NamLLeff)
# 
# NamLLeffSummary<-NamLLeff %>% 
#   group_by(Year) %>%
#   summarise(N=sum(`Number of HOOKS_SET`))
# 
# NamLLeffSummary
# 
# 
# 
# ### load trawling effort from Namibia
# ### this combines both wet and frozen fish
# 
# NamTRwet<-fread("Trawl_effort_Namibia_wet.csv")
# NamTRfrozen<-fread("Trawl_effort_Namibia_frozen.csv")
# head(NamTRfrozen)
# head(NamTRwet)
# names(NamTRfrozen)<-names(NamTRwet)
# 
# NamTrawlSummary<-NamTRwet %>% bind_rows(NamTRfrozen) %>% 
#   group_by(YEAR) %>%
#   summarise(N=sum(`DURATION(HOURS)`))
# NamTrawlSummary
# 
# 
# 
# 



#########################################################################
# LOAD FISHERY DATA PROVIDED BY NAMIBIAN FISHERIES DEPARTMENT
#########################################################################
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\BycatchData"), silent=T)
nhooks<-fread("Namibia_DemersalLongLine.csv")
head(nhooks)
nhooksSummary<-nhooks %>% 
  filter(Year>1999 & Year<2018) %>%
  group_by(Year) %>%
  summarise(N=sum(HOOKS_SET))

## scale
longlineNAM<- (nhooksSummary$N-mean(nhooksSummary$N))/sd(nhooksSummary$N)






### PLOT THE THREE DATASETS TOGETHER ON ONE GRAPH
plotdat<- data.frame(Year=rep(seq(2000,2017), 3), N=c(longlineJR, longlineICCAT, longlineNAM), Source=rep(c("Joel Rice","ICCAT","Namibia"), each=18))

ggplot(data=plotdat) +
  geom_line(aes(x=Year, y=N, col=Source), size=2) +
  scale_x_continuous(name="Year", limits=c(2000,2018), breaks=seq(2000,2018,2)) +
  scale_y_continuous(name="Longline fishing effort (standardized)", limits=c(-3,3), breaks=seq(-3,3,0.5)) +
  theme(panel.background=element_rect(fill="white", colour="black"),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black", angle=45, vjust=0.5),
        axis.title=element_text(size=20),
        strip.text.x=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())



# ### ESTIMATE CORRELATION WITH SURVIVAL ESTIMATES
# ## THERE IS A REASONABLE POSITIVE CORRELATION
# 
# try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM"), silent=T)
# AYNA<-fread("AYNA_Gough_IPM_estimates.csv")
# head(AYNA)
# surv<-AYNA %>% filter(parameter=="adult.survival") %>%
#   filter(Year!='mean') %>%
#   select(Year,Mean,Median)
# 
# cor.test(surv$Median,log(nhooksSummary$N))
# plot(surv$Median~log(nhooksSummary$N))
# #cor.test(surv$Median,lag(nhooksSummary$N,2))
# #cor.test(surv$Median,lead(nhooksSummary$N,2))
# #cor.test(surv$Median[8:17],NamLLeffSummary$N)
# 
# 
# # surv<-AYNA %>% filter(parameter=="adult.survival") %>%
# #   filter(Year %in% c('2009','2010','2016','2017')) %>%
# #   select(Year,Mean,Median)
# # 
# # cor.test(surv$Mean,NamTrawlSummary$N)




#########################################################################
# INCORPORATE BYCATCH MITIGATION ADOPTION
#########################################################################

mitigation=c(rep(1,13),0.8,0.6,0.4,0.2,0.1)
nhooksSummary$Neff<-nhooksSummary$N*mitigation

## scale
longlineNAM<- (nhooksSummary$Neff-mean(nhooksSummary$Neff))/sd(nhooksSummary$Neff)


#########################################################################
# LOAD PRE-PREPARED DATA
#########################################################################
### see 'IPM_DATA_PREPARATION.R' for details on how data are aggregated


#### CMR SURVIVAL DATA ######

try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM"), silent=T)
AYNA<-fread("AYNA_simple_encounter_history_1982_2018.csv")
names(AYNA)
CH<-as.matrix(AYNA[,3:39], dimnames=F)
AYNA$AGE[is.na(AYNA$AGE)]<-1    ## set all NA as 'adult'

### check that there are contacts in every season
apply(CH,2,sum)




#### COUNT DATA FOR POPULATION TREND ######

AYNA.pop<-fread("AYNA_pop_counts_1982_2018.csv")
AYNA.pop<-subset(AYNA.pop,Year>1999)	## reduce data set to remove NA in 4 years
n.years<-dim(AYNA.pop)[1]		## defines the number of years
n.sites<-dim(AYNA.pop)[2]-1 ## defines the number of study areas
str(AYNA.pop)
names(AYNA.pop)



#### BREEDING SUCCESS DATA FOR FECUNDITY ######

AYNA.bs<-fread("AYNA_breed_success_1982_2017.csv")
AYNA.bs<-subset(AYNA.bs,Year>1999)	## reduce data set to specified time period
J<-as.integer(AYNA.bs$n_nests*AYNA.bs$BREED_SUCC)
R<-AYNA.bs$n_nests




#########################################################################
# MANIPULATE DATA: CREATE MATRIX OF AGE FOR EACH OCCASION AND INDIVIDUAL
#########################################################################

## this matrix will relate to the survival parameter estimates chosen in the model
## simple model only has 2 survival parameters:
## 1 - juvenile and immature survival (years 1-5)
## 2 - adult survival (birds >5 years old)

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x==1))
f <- apply(CH, 1, get.first)


## REMOVE BIRDS THAT ARE TOO YOUNG TO HAVE HAD A CHANCE TO RETURN
tooyoung<-ifelse(f>(dim(CH)[2]-5),ifelse(AYNA$AGE==0,1,0),0)
CH<-CH[tooyoung==0,]  ## removes individuals that were ringed as chicks <5 years before end of time series
f <- apply(CH, 1, get.first)
toolate<-ifelse(f==dim(CH)[2],1,0)
CH<-CH[toolate==0,]  ## removes individuals ringed in last occasion end of time series
ages<-AYNA$AGE[tooyoung==0]

## CREATE BLANK AGE MATRIX
AGEMAT<-matrix(2,nrow=nrow(CH),ncol=ncol(CH))
n.occ<-ncol(CH)

## LOOP OVER EACH BIRD RINGED AND SET PRE-CAPTURE DATA TO NA AND ADJUST AGE
for (l in 1:nrow(AGEMAT)){
  firstocc<-get.first(CH[l,])
  lastjuv<-firstocc+4
  lastjuv<-ifelse(lastjuv>n.occ,n.occ,lastjuv)
  young<-ages[l]
  if(firstocc>1){AGEMAT[l,1:(firstocc-1)]<-NA}  ## sets everything before first contact to NA
  if(young==0){AGEMAT[l,firstocc:lastjuv]<-1}  ## sets all juvenile years to 1
}

### CHECK WHETHER IT LOOKS OK ###
head(AGEMAT)
head(CH)


#########################################################################
# MANIPULATE DATA: RE-ARRANGE DATA TO REMOVE CONTACTS BEFORE 2000 (and individuals with no contacts after 2000)
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
dim(AGEMAT)
head(rCH)
head(AGEMAT)


## PREPARE CONSTANTS
n.ind<-dim(rCH)[1]		## defines the number of individuals
n.years<-dim(rCH)[2]  ## defines the number of years
f <- apply(rCH, 1, get.first)




#########################################################################
# MANIPULATE DATA: INITIAL VALUES
#########################################################################


## CREATE MATRIX for INITIAL STATE Z FOR SURVIVAL MODEL
zinit<-rCH
for (l in 1:nrow(zinit)){
  firstocc<-get.first(zinit[l,])
  zinit[l,1:firstocc]<-NA  ## sets everything up to first contact to NA
  zinit[l,(firstocc+1):n.years]<-1  ## alive after first contact
}
dim(zinit)


## CREATE MATRIX for COUNTS

N.init=matrix(NA, nrow=n.years,ncol=n.sites)
N.init[1,]<-as.matrix(AYNA.pop[1,2:12])






#########################################################################
# SPECIFY MODEL IN JAGS
#########################################################################
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM")
sink("AYNA_IPM_projection_scenario_0.jags")
cat("

  
    model {
    #-------------------------------------------------
    # integrated population model for the Gough AYNA population
    # - age structured model with 6 age classes 
    # - adult survival based on CMR ringing data
    # - pre breeding census, female-based assuming equal sex ratio & survival
    # - productivity based on Area 1 nest monitoring data
    # - simplified population process with informed prior for adults skipping breeding and uninformed immatures recruiting
    # - FOUR future scenarios to project population growth after eradication, bycatch mitigation or no management
    # -------------------------------------------------
    
    #-------------------------------------------------  
    # 1. PRIORS FOR ALL DATA SETS
    #-------------------------------------------------
    
    
    # -------------------------------------------------        
    # 1.1. Priors and constraints FOR FECUNDITY
    # -------------------------------------------------
    
    for (t in 1:T){  
      ann.fec[t] ~ dunif(0,1)           # Priors on fecundity can range from 0-1 chicks per pair (constrained based on our data)
      imm.rec[t]~dunif(0,1)                ## RECRUITMENT PROBABILITY COULD SET MORE INFORMATIVE PRIOR HERE
      skip.prob[t]~dunif(0,1)              ## PRIOR FOR ADULT BREEDER SKIPPING PROBABILITY from Cuthbert paper that reported breeding propensity of 0.66
    } #t
    
    
    
    
    # -------------------------------------------------        
    # 1.2. Priors and constraints FOR POPULATION COUNTS
    # -------------------------------------------------
    
    
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
    
    
    # -------------------------------------------------        
    # 1.3. Priors and constraints FOR SURVIVAL
    # -------------------------------------------------
    
    ### RECAPTURE PROBABILITY
    mean.p ~ dunif(0, 1)                          # Prior for mean recapture
    logit.p <- log(mean.p / (1-mean.p))           # Logit transformation
    
    for (t in 1:T){
      logit(p[t]) <- logit.p  + capt.raneff[t]
      capt.raneff[t] ~ dnorm(0, tau.capt)
    }
    
    ### SURVIVAL PROBABILITY
    for (i in 1:nind){
      for (t in f[i]:(T-1)){
        logit(phi[i,t]) <- mu[AGEMAT[i,t]] + surv.raneff[t] + bycatch*longline[t]
      } #t
    } #i
    
    
    ## AGE-SPECIFIC SURVIVAL 
    for (age in 1:2){
      beta[age] ~ dunif(0, 1)                         # Priors for age-specific survival
      mu[age] <- log(beta[age] / (1-beta[age]))       # Logit transformation
    }
    
    ## RANDOM TIME EFFECT ON SURVIVAL 
    for (t in 1:(T-1)){
      surv.raneff[t] ~ dnorm(0, tau.surv)
    }
    
    ### PRIORS FOR RANDOM EFFECTS
    sigma.surv ~ dunif(0, 10)                     # Prior for standard deviation of survival
    tau.surv <- pow(sigma.surv, -2)
    
    sigma.capt ~ dunif(0, 10)                     # Prior for standard deviation of capture
    tau.capt <- pow(sigma.capt, -2)
    
    
    ### PRIOR FOR BYCATCH EFFECTS
    bycatch ~ dnorm(0,tau.byc)
    sigma.byc ~ dunif(0, 10)                     # Prior for standard deviation of capture    
    tau.byc <- pow(sigma.byc, -2)
    
    
    #-------------------------------------------------  
    # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
    #-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------
    
    for (tt in 2:T){
    
      ## THE PRE-BREEDING YEARS ##
    
      nestlings[tt] <- ann.fec[tt] * 0.5 * Ntot.breed[tt]                                                     ### number of locally produced FEMALE chicks
      JUV[tt] ~ dpois(nestlings[tt])                                                                     ### need a discrete number otherwise dbin will fail, dpois must be >0
      N1[tt]  ~ dbin(ann.surv[1,tt-1], round(JUV[tt-1]))                                                    ### number of 1-year old survivors 
      N2[tt] ~ dbin(ann.surv[1,tt-1], round(N1[tt-1]))                                                      ### number of 2-year old survivors
      N3[tt] ~ dbin(ann.surv[1,tt-1], round(N2[tt-1]))                                                       ### number of 3-year old survivors
      N4[tt] ~ dbin(ann.surv[1,tt-1], round(N3[tt-1]))                                                       ### number of 4-year old survivors
      N5[tt] ~ dbin(ann.surv[1,tt-1], round(N4[tt-1]))                                                       ### number of 5-year old survivors
    
    
      ## THE POTENTIAL RECRUITING YEARS ##
    
      N6[tt] ~ dbin(ann.surv[2,tt-1], round(N5[tt-1]))                                     ### number of 6-year old survivors that are ready for recruitment - using adult survival
      N.notrecruited[tt] ~ dbin(ann.surv[2,tt-1], round(max(10,non.recruits[tt-1])))       ### number of not-yet-recruited birds surviving from previous year
      non.recruits[tt]<-(N6[tt]+N.notrecruited[tt])-ann.recruits[tt]                      ## number of birds that do not recruit is the sum of all available minus the ones that do recruit
    
    
      ## THE BREEDING YEARS ##
    
      Ntot.breed[tt] ~ dpois(pop.size[tt])                                           ### the annual number of breeding birds is the estimate from the count SSM
      ann.recruits[tt] ~ dbin(imm.rec[tt],round(N6[tt]+N.notrecruited[tt]))          ### Ntot.breed[tt]-Nold.breed[tt]+1))           ### this total number comprises a bunch of new recruits, which is the number of total breeders that are not old breeders
      Nold.breed[tt]<- N.pot.breed[tt]-N.non.breed[tt]                              ### number of old breeders is survivors from previous year minus those that skip a year of breeding
      N.pot.breed[tt] ~ dbin(ann.surv[2,tt-1], round(sum(Ntot.breed[tt-1],N.non.breed[tt-1])))   ### number of potential old breeders is the number of survivors from previous year breeders and nonbreeders
      N.non.breed[tt] ~ dbin(skip.prob[tt], round(N.pot.breed[tt]))                             ### number of old nonbreeders (birds that have bred before and skip breeding) 
    
    } # tt
    
    
    
    ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1 - based on stable stage distribution from previous model
    
    JUV[1]<-round(Ntot.breed[1]*0.5*ann.fec[1])
    N1[1]<-round(Ntot.breed[1]*0.17574058)
    N2[1]<-round(Ntot.breed[1]*0.11926872)
    N3[1]<-round(Ntot.breed[1]*0.10201077)
    N4[1]<-round(Ntot.breed[1]*0.08725001)
    N5[1]<-round(Ntot.breed[1]*0.07462511)
    non.recruits[1]<-round(Ntot.breed[1]*0.3147774)
    Ntot.breed[1]<-sum(y.count[1,])
    N.non.breed[1]<- round(Ntot.breed[1]*0.12632740)
    
    
    
    # -------------------------------------------------        
    # 2.2. Observation process for population counts: state-space model of annual counts
    # -------------------------------------------------
    
    for (s in 1:n.sites){			### start loop over every study area
    
      ## State process for entire time series
    
      for (t in 1:(T-1)){
        lambda[t,s] ~ dnorm(mean.lambda[s], tau.proc[s])								# Distribution for random error of growth rate
        N.est[t+1,s]<-N.est[t,s]*lambda[t,s]										        # Linear predictor (population size based on past pop size and change rate)
      }														# run this loop over nyears
    
    
      ## Observation process
    
      for (t in 1:T){
        y.count[t,s] ~ dnorm(N.est[t,s], tau.obs[s])								# Distribution for random error in observed numbers (counts)
      }														# run this loop over t= nyears
    }		## end site loop
    
    
    
    # -------------------------------------------------        
    # 2.3. Likelihood for fecundity: Poisson regression from the number of surveyed broods
    # -------------------------------------------------
    for (t in 1:(T-1)){
      J[t] ~ dpois(rho.fec[t])
      rho.fec[t] <- R[t]*ann.fec[t]
    } #	close loop over every year in which we have fecundity data
    
    
    
    
    # -------------------------------------------------        
    # 2.4. Likelihood for adult and juvenile survival from CMR
    # -------------------------------------------------
    
    # Likelihood 
    for (i in 1:nind){

      # Define latent state at first capture
      z[i,f[i]] <- 1

      for (t in (f[i]+1):T){
    
        # State process
        z[i,t] ~ dbern(mu1[i,t])
        mu1[i,t] <- phi[i,t-1] * z[i,t-1]
    
        # Observation process
        y[i,t] ~ dbern(mu2[i,t])
        mu2[i,t] <- p[t] * z[i,t]
      } #t
    } #i
    
    
    
    
    #-------------------------------------------------  
    # 3. DERIVED PARAMETERS FOR OUTPUT REPORTING
    #-------------------------------------------------
    
    ## DERIVED SURVIVAL PROBABILITIES PER YEAR 
    for (t in 1:(T-1)){
      for (age in 1:2){
        logit(ann.surv[age,t]) <- mu[age] + surv.raneff[t]
      }
    }
    
    
    ## DERIVED POPULATION SIZE PER YEAR 
    for (t in 1:T){
      pop.size[t]<-max(10,sum(N.est[t,1:n.sites]))               ## introduced max to prevent this number from being 0 which leads to invalid parent error on Ntot.breed
    }
    
    
    ## DERIVED OVERALL POPULATION GROWTH RATE 
    pop.growth.rate <- mean(lambda[1:(T-1),1:n.sites])  				# Arithmetic mean for whole time series
    
    ## DERIVED MEAN FECUNDITY 
    mean.fec <- mean(ann.fec)
    sd.fec <- sd(ann.fec)
    tau.fec <- pow(max(sd.fec,0.01),-2)
    
    #-------------------------------------------------  
    # 4. PROJECTION INTO FUTURE
    #-------------------------------------------------


    for (tt in (T+1):FUT.YEAR){
    
      ## RANDOMLY DRAW DEMOGRAPHIC RATES FROM PREVIOUS YEARS WHILE AVOIDING THAT INDEX BECOMES 0
    
      FUT[tt] ~ dunif(1.5,15.5)           ### CHANGE FROM 1.5 to 15.5 to only sample from last three years when survival was high
      FUT.int[tt]<-round(FUT[tt])
      fut.fec[tt] ~ dnorm(mean.fec,tau.fec)   ### CHANGE FROM mean.fec to 0.69 + 0.16 from Caravaggi et al. 2018 for eradication scenario
    
    
    
      # -------------------------------------------------        
      # 4.1. System process for future
      # -------------------------------------------------
    
      ## THE PRE-BREEDING YEARS ##
    
      nestlings[tt] <- round(fut.fec[tt]* 0.5 * Ntot.breed[tt])                                             ### number of locally produced FEMALE chicks based on average fecundity - to use just one take ann.fec[FUT.int[tt]] 
      N1[tt]  ~ dbin(ann.surv[1,FUT.int[tt]-1], max(1,round(nestlings[tt-1])))                                                    ### number of 1-year old survivors 
      N2[tt] ~ dbin(ann.surv[1,FUT.int[tt]-1], round(N1[tt-1]))                                                      ### number of 2-year old survivors
      N3[tt] ~ dbin(ann.surv[1,FUT.int[tt]-1], round(N2[tt-1]))                                                       ### number of 3-year old survivors
      N4[tt] ~ dbin(ann.surv[1,FUT.int[tt]-1], round(N3[tt-1]))                                                       ### number of 4-year old survivors
      N5[tt] ~ dbin(ann.surv[1,FUT.int[tt]-1], round(N4[tt-1]))                                                       ### number of 5-year old survivors
    
    
      ## THE POTENTIAL RECRUITING YEARS ##
    
      N6[tt] ~ dbin(ann.surv[2,FUT.int[tt]-1], round(N5[tt-1]))                                     ### number of 6-year old survivors that are ready for recruitment - using adult survival
      N.notrecruited[tt] ~ dbin(ann.surv[2,FUT.int[tt]-1], round(max(10,non.recruits[tt-1])))       ### number of not-yet-recruited birds surviving from previous year
      non.recruits[tt]<-(N6[tt]+N.notrecruited[tt])-ann.recruits[tt]                                ### number of birds that do not recruit is the sum of all available minus the ones that do recruit
      ann.recruits[tt] ~ dbin(imm.rec[FUT.int[tt]],round(N6[tt]+N.notrecruited[tt]))                       ### new recruits
    
    
      ## THE BREEDING YEARS ##
    
      Ntot.breed[tt] <- Nold.breed[tt] + ann.recruits[tt]                                         ### the annual number of breeding birds is the estimate from the count SSM
      Nold.breed[tt]<- N.pot.breed[tt]-N.non.breed[tt]                                            ### number of old breeders is survivors from previous year minus those that skip a year of breeding
      N.pot.breed[tt] ~ dbin(ann.surv[2,FUT.int[tt]-1], round(sum(Ntot.breed[tt-1],N.non.breed[tt-1])))   ### number of potential old breeders is the number of survivors from previous year breeders and nonbreeders
      N.non.breed[tt] ~ dbin(skip.prob[FUT.int[tt]], round(N.pot.breed[tt]))                             ### number of old nonbreeders (birds that have bred before and skip breeding) 
    
    
      ## CALCULATE ANNUAL POP GROWTH RATE ##
      fut.lambda[tt-19] <- Ntot.breed[tt]/max(1,Ntot.breed[tt-1])                                 ### inserted safety to prevent denominator being 0
    
    } # tt
    
    # -------------------------------------------------        
    # 4.2. DERIVED POPULATION GROWTH RATE FOR FUTURE
    # -------------------------------------------------
    
    ## DERIVED OVERALL POPULATION GROWTH RATE 
    future.growth.rate <- mean(fut.lambda[1:10])  				# projected ANNUAL growth rate in the future 
    
    }
    
    
    
    ",fill = TRUE)
sink()





#########################################################################
# PREPARE DATA FOR MODEL
#########################################################################

# Bundle data
jags.data <- list(y = rCH,
                  f = f,
                  T = n.years,
                  nind = n.ind,
                  AGEMAT=AGEMAT,
                  
                  ### count data
                  n.sites=n.sites,
                  y.count=as.matrix(AYNA.pop[,2:12]),
                  
                  ### breeding success data
                  J=J,
                  R=R,
                  
                  ### longline effort data
                  longline=longlineICCAT,
                  
                  ### FUTURE PROJECTION
                  FUT.YEAR=n.years+10,
                  FUT.int=c(seq(1,(n.years-1),1),rep(NA,11)),
                  fut.fec=c(rep(0.5,(n.years)),rep(NA,10))     ## blank vector to hold index for future demographic rates
                  )


# Initial values 
inits <- function(){list(beta = runif(2, 0, 1),
                         z = zinit,
                         mean.p = runif(1, 0, 1),
                         bycatch = rnorm(1,0,0.01),
                         #hookpod = rnorm(1,0,0.01),
                         
                         ### count data
                         sigma.proc=runif(n.sites,0,5),
                         mean.lambda=runif(n.sites,0.1,2),
                         sigma.obs=runif(n.sites,0,10),
                         N.est=N.init)}
 

# Parameters monitored
parameters <- c("Ntot.breed","ann.fec","skip.prob","imm.rec","ann.surv","beta","pop.growth.rate","future.growth.rate","mean.fec","bycatch")  #,"hookpod"

# MCMC settings
ni <- 30000
nt <- 4
nb <- 10000
nc <- 4






# RUN THE FOUR SCENARIOS
AYNAscenario0 <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\AYNA_IPM_projection_scenario_0.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)
#AYNAscenarioM <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\AYNA_IPM_projection_scenarioM.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)
#AYNAscenarioB <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\AYNA_IPM_projection_scenarioB.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)
#AYNAscenarioMB <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\AYNA_IPM_projection_scenarioMB.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)


# COMPARE WITH BYCATCH MITIGATION PROPORTION AS INPUT
jags.data$longline<-mitigation
AYNAscenario0byc <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\AYNA_IPM_projection_scenario_0.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)


#########################################################################
# SAVE OUTPUT - RESULT PROCESSING in AYNA_IPM_result_summaries.r
#########################################################################
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM")
save.image("AYNA_IPM_output_4scenarios.RData")





















#########################################################################
# THIS MODEL DID NOT WORK - ALWAYS INVALID PARENT ERROR SOMEWHERE
#########################################################################

sink("AYNA_IPM_projection_v4.jags")
cat("
    model {
    #-------------------------------------------------
    # integrated population model for the Gough AYNA population
    # - age structured model with 6 age classes 
    # - adult survival based on CMR ringing data
    # - pre breeding census, female-based assuming equal sex ratio & survival
    # - productivity based on Area 1 nest monitoring data
    # - simplified population process with informed prior for adults skipping breeding and uninformed immatures recruiting
    # -------------------------------------------------
    
    #-------------------------------------------------  
    # 1. PRIORS FOR ALL DATA SETS
    #-------------------------------------------------
    
    
    # -------------------------------------------------        
    # 1.1. Priors and constraints FOR FECUNDITY
    # -------------------------------------------------
    
    for (t in 1:T){  
    ann.fec[t] ~ dunif(0.2,0.8)           # Priors on fecundity can range from 0-1 chicks per pair (constrained based on our data)
    imm.rec[t]~dunif(0.15,0.75)                ## RECRUITMENT PROBABILITY COULD SET MORE INFORMATIVE PRIOR HERE
    skip.prob[t]~dunif(0.15,0.45)              ## PRIOR FOR ADULT BREEDER SKIPPING PROBABILITY from Cuthbert paper that reported breeding propensity of 0.66
    } #t
    
    
    
    
    # -------------------------------------------------        
    # 1.2. Priors and constraints FOR POPULATION COUNTS
    # -------------------------------------------------
    
    
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
    
    
    # -------------------------------------------------        
    # 1.3. Priors and constraints FOR SURVIVAL
    # -------------------------------------------------
    
    ### RECAPTURE PROBABILITY
    mean.p ~ dunif(0, 1)                          # Prior for mean recapture
    logit.p <- log(mean.p / (1-mean.p))           # Logit transformation
    
    for (t in 1:T){
    logit(p[t]) <- logit.p  + capt.raneff[t]
    capt.raneff[t] ~ dnorm(0, tau.capt)
    }
    
    ### SURVIVAL PROBABILITY
    for (i in 1:nind){
    for (t in f[i]:(T-1)){
    logit(phi[i,t]) <- mu[AGEMAT[i,t]] + surv.raneff[t]
    } #t
    } #i
    
    
    ## AGE-SPECIFIC SURVIVAL 
    for (age in 1:2){
    beta[age] ~ dunif(0, 1)                         # Priors for age-specific survival
    mu[age] <- log(beta[age] / (1-beta[age]))       # Logit transformation
    }
    
    ## RANDOM TIME EFFECT ON SURVIVAL 
    for (t in 1:(T-1)){
    surv.raneff[t] ~ dnorm(0, tau.surv)
    }
    
    ### PRIORS FOR RANDOM EFFECTS
    sigma.surv ~ dunif(0, 10)                     # Prior for standard deviation of survival
    tau.surv <- pow(sigma.surv, -2)
    
    sigma.capt ~ dunif(0, 10)                     # Prior for standard deviation of capture
    tau.capt <- pow(sigma.capt, -2)
    
    
    
    
    
    
    #-------------------------------------------------  
    # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
    #-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------
    
    for (tt in 2:T){
    
    ## THE PRE-BREEDING YEARS ##
    
    nestlings[tt] <- ann.fec[tt] * 0.5 * Ntot.breed[tt]                                                     ### number of locally produced FEMALE chicks
    JUV[tt] ~ dpois(nestlings[tt])                                                                     ### need a discrete number otherwise dbin will fail, dpois must be >0
    N1[tt]  ~ dbin(ann.surv[1,tt-1], round(JUV[tt-1]))                                                    ### number of 1-year old survivors 
    N2[tt] ~ dbin(ann.surv[1,tt-1], round(N1[tt-1]))                                                      ### number of 2-year old survivors
    N3[tt] ~ dbin(ann.surv[1,tt-1], round(N2[tt-1]))                                                       ### number of 3-year old survivors
    N4[tt] ~ dbin(ann.surv[1,tt-1], round(N3[tt-1]))                                                       ### number of 4-year old survivors
    N5[tt] ~ dbin(ann.surv[1,tt-1], round(N4[tt-1]))                                                       ### number of 5-year old survivors
    
    
    ## THE POTENTIAL RECRUITING YEARS ##
    
    N6[tt] ~ dbin(ann.surv[1,tt-1], round(N5[tt-1]))                                     ### number of 6-year old survivors that are ready for recruitment
    N.notrecruited[tt] ~ dbin(ann.surv[2,tt-1], round(max(10,non.recruits[tt-1])))       ### number of not-yet-recruited birds surviving from previous year
    non.recruits[tt]<-(N6[tt]+N.notrecruited[tt])-ann.recruits[tt]                      ## number of birds that do not recruit is the sum of all available minus the ones that do recruit
    
    
    ## THE BREEDING YEARS ##
    
    Ntot.breed[tt] ~ dpois(pop.size[tt])                                           ### the annual number of breeding birds is the estimate from the count SSM
    ann.recruits[tt] ~ dbin(imm.rec[tt],round(N6[tt]+N.notrecruited[tt]))          ### Ntot.breed[tt]-Nold.breed[tt]+1))           ### this total number comprises a bunch of new recruits, which is the number of total breeders that are not old breeders
    Nold.breed[tt]<- N.pot.breed[tt]-N.non.breed[tt]                              ### number of old breeders is survivors from previous year minus those that skip a year of breeding
    N.pot.breed[tt] ~ dbin(ann.surv[2,tt-1], round(sum(Ntot.breed[tt-1],N.non.breed[tt-1])))   ### number of potential old breeders is the number of survivors from previous year breeders and nonbreeders
    N.non.breed[tt] ~ dbin(skip.prob[tt], round(N.pot.breed[tt]))                             ### number of old nonbreeders (birds that have bred before and skip breeding) 
    
    
    
    } # tt
    
    
    
    ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1 - based on stable stage distribution from previous model
    
    JUV[1]<-round(Ntot.breed[1]*0.5*ann.fec[1])
    N1[1]<-round(Ntot.breed[1]*0.17574058)
    N2[1]<-round(Ntot.breed[1]*0.11926872)
    N3[1]<-round(Ntot.breed[1]*0.10201077)
    N4[1]<-round(Ntot.breed[1]*0.08725001)
    N5[1]<-round(Ntot.breed[1]*0.07462511)
    non.recruits[1]<-round(Ntot.breed[1]*0.3147774)
    Ntot.breed[1]<-sum(y.count[1,])
    N.non.breed[1]<- round(Ntot.breed[1]*0.12632740)
    
    
    
    # -------------------------------------------------        
    # 2.2. Observation process for population counts: state-space model of annual counts
    # -------------------------------------------------
    
    for (s in 1:n.sites){			### start loop over every study area
    
    ## State process for entire time series
    
    for (t in 1:(T-1)){
    lambda[t,s] ~ dnorm(mean.lambda[s], tau.proc[s])								# Distribution for random error of growth rate
    N.est[t+1,s]<-N.est[t,s]*lambda[t,s]										        # Linear predictor (population size based on past pop size and change rate)
    }														# run this loop over nyears
    
    
    ## Observation process
    
    for (t in 1:T){
    y.count[t,s] ~ dnorm(N.est[t,s], tau.obs[s])								# Distribution for random error in observed numbers (counts)
    }														# run this loop over t= nyears
    }		## end site loop
    
    
    
    # -------------------------------------------------        
    # 2.3. Likelihood for fecundity: Poisson regression from the number of surveyed broods
    # -------------------------------------------------
    for (t in 1:(T-1)){
    J[t] ~ dpois(rho.fec[t])
    rho.fec[t] <- R[t]*ann.fec[t]
    } #	close loop over every year in which we have fecundity data
    
    
    
    
    # -------------------------------------------------        
    # 2.4. Likelihood for adult and juvenile survival from CMR
    # -------------------------------------------------
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):T){
    # State process
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t] <- phi[i,t-1] * z[i,t-1]
    
    # Observation process
    y[i,t] ~ dbern(mu2[i,t])
    mu2[i,t] <- p[t] * z[i,t]
    } #t
    } #i
    
    
    
    
    #-------------------------------------------------  
    # 3. DERIVED PARAMETERS FOR OUTPUT REPORTING
    #-------------------------------------------------
    
    ## DERIVED SURVIVAL PROBABILITIES PER YEAR 
    for (t in 1:(T-1)){
      for (age in 1:2){
        logit(ann.surv[age,t]) <- mu[age] + surv.raneff[t]
      }
    }
    
    
    ## DERIVED POPULATION SIZE PER YEAR 
    for (t in 1:T){
      pop.size[t]<-sum(N.est[t,1:n.sites])
    }
    
    
    ## DERIVED OVERALL POPULATION GROWTH RATE 
    pop.growth.rate <- mean(lambda[1:(T-1),1:n.sites])  				# Arithmetic mean for whole time series
    
    
    
    #-------------------------------------------------  
    # 4. PROJECTION INTO FUTURE
    #-------------------------------------------------
    
    
    for (tt in (T+1):FUT.YEAR){
    
    ## RANDOMLY DRAW DEMOGRAPHIC RATES FROM MEAN AND SD OF PREVIOUS YEARS 
    
    #ann.fec[tt] ~ dnorm(mean(ann.fec[1:(T-1)]),pow(sd(ann.fec[1:(T-1)]),-2))     ### this always causes an invalid parent error
    #ann.surv[1,tt-1] ~ dnorm(mean(ann.surv[1,1:(T-1)]),pow(sd(ann.surv[1,1:(T-1)]),-2))
    #ann.surv[2,tt-1] ~ dnorm(mean(ann.surv[2,1:(T-1)]),pow(sd(ann.surv[2,1:(T-1)]),-2))
    
    ann.fec[tt] ~ dunif(0.3,0.75)            ### the model has no value for ann.fec[19] because at the time of writing the chicks had not fledged (March 2019)
    ann.surv[1,tt-1] ~ dnorm(mu[1],0.01)
    ann.surv[2,tt-1] ~ dnorm(mu[2],0.05)
    imm.rec[tt] ~ dnorm(mean(imm.rec[1:(T-1)]),pow(sd(imm.rec[1:(T-1)]),-2))
    #imm.rec[tt]~dunif(0.25,0.75)     
    skip.prob[tt] ~ dnorm(mean(skip.prob[1:(T-1)]),pow(sd(skip.prob[1:(T-1)]),-2))
    
    
    
    # -------------------------------------------------        
    # 4.1. System process for future
    # -------------------------------------------------
    
    ## THE PRE-BREEDING YEARS ##
    
    nestlings[tt] <- round(ann.fec[tt] * 0.5 * Ntot.breed[tt])                                             ### number of locally produced FEMALE chicks
    N1[tt]  ~ dbin(ann.surv[1,tt-1], round(nestlings[tt-1]))                                                    ### number of 1-year old survivors 
    N2[tt] ~ dbin(ann.surv[1,tt-1], round(N1[tt-1]))                                                      ### number of 2-year old survivors
    N3[tt] ~ dbin(ann.surv[1,tt-1], round(N2[tt-1]))                                                       ### number of 3-year old survivors
    N4[tt] ~ dbin(ann.surv[1,tt-1], round(N3[tt-1]))                                                       ### number of 4-year old survivors
    N5[tt] ~ dbin(ann.surv[1,tt-1], round(N4[tt-1]))                                                       ### number of 5-year old survivors
    
    
    ## THE POTENTIAL RECRUITING YEARS ##
    
    N6[tt] ~ dbin(ann.surv[1,tt-1], round(N5[tt-1]))                                     ### number of 6-year old survivors that are ready for recruitment
    N.notrecruited[tt] ~ dbin(ann.surv[2,tt-1], round(max(10,non.recruits[tt-1])))       ### number of not-yet-recruited birds surviving from previous year
    non.recruits[tt]<-(N6[tt]+N.notrecruited[tt])-ann.recruits[tt]                                ### number of birds that do not recruit is the sum of all available minus the ones that do recruit
    ann.recruits[tt] ~ dbin(imm.rec[tt],round(N6[tt]+N.notrecruited[tt]))                       ### new recruits
    
    
    ## THE BREEDING YEARS ##
    
    Ntot.breed[tt] <- Nold.breed[tt] + ann.recruits[tt]                                         ### the annual number of breeding birds is the estimate from the count SSM
    Nold.breed[tt]<- N.pot.breed[tt]-N.non.breed[tt]                                            ### number of old breeders is survivors from previous year minus those that skip a year of breeding
    N.pot.breed[tt] ~ dbin(ann.surv[2,tt-1], round(sum(Ntot.breed[tt-1],N.non.breed[tt-1])))   ### number of potential old breeders is the number of survivors from previous year breeders and nonbreeders
    N.non.breed[tt] ~ dbin(skip.prob[tt], round(N.pot.breed[tt]))                             ### number of old nonbreeders (birds that have bred before and skip breeding) 
    
    } # tt
    
    # -------------------------------------------------        
    # 4.2. DERIVED POPULATION GROWTH RATE FOR FUTURE
    # -------------------------------------------------
    
    ## DERIVED OVERALL POPULATION GROWTH RATE 
    future.growth.rate <- pow(Ntot.breed[FUT.YEAR]/Ntot.breed[T],-10)  				# projected ANNUAL growth rate in the future 
    
    }
    ",fill = TRUE)
sink()




