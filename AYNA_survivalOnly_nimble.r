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
#setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM")
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
  
  agebeta ~ dnorm(1, sd = 0.001)    # Prior for shape of increase in juvenile recapture probability with age

  # TODO 
  # fix this
  # beta.ICCAT.ll.e ~ dnorm(0, 1)  # TODO - change precison?
  # beta.ICCAT.ll.mit ~ dnorm(0, 1)  # TODO - change precison?
  # beta.Nam.ll.mit ~ dnorm(0, 1) # TODO - change precison?
  # beta.SA.ll.mit ~ dnorm(0, 1) # TODO - change precison?
  # beta.Uru.ll.mit ~ dnorm(0, 1) # TODO - change precison?
  
  ## RANDOM TIME EFFECT ON RESIGHTING PROBABILITY OF JUVENILES
  for (t in 1:(n.occasions-1)){
    for (j in 1:t){ ## zero by definition (these are never actually used)
      p.juv[t,j] <- 0
    }
    for (j in (t+1):(n.occasions-1)){
      logit(p.juv[t,j])  <- mu.p.juv[goodyear[j]] + agebeta*(j - t)/2 + eps.p[j]
    }
  }
  
  ## PRIORS FOR RANDOM EFFECTS
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
            marr.a = adult.marray) 

#### CONSTANTS ####
const <- list(n.occasions = length(start:2021),
              r.j=apply(chick.marray,1,sum),
              r.a=apply(adult.marray,1,sum),
              goodyear=goodyears$p.sel,
              juv.poss=phi.juv.possible$JuvSurv#, ### sets the annual survival of juveniles to the mean if <70 were ringed
)

#### INITIAL VALUES ####

inits <- list(sigma.phi = rexp(1, 1),
              mean.phi.ad = rbeta(1, 91,9) ,
              mean.phi.juv = rbeta(1, 75.7, 24.3),
              sigma.p = rexp(1, 1), 
              agebeta = rnorm(1, 0, 0.001),
              mean.p.ad = c(runif(1, 0.05, 0.5), runif(1, 0.2, 1)), 
              mu.p.juv = rnorm(2, -4, 0.25)
)

#### PARAMETERS TO MONITOR ####
# TODO check that this has everything 
params <- c("mean.phi.ad","mean.phi.juv",
            "sigma.phi", "sigma.p",
            "agebeta",
            "phi.ad","phi.juv",    
            "mean.p.ad", "mu.p.juv"
            )

#### MCMC SETTINGS ####
nb <- 25000 #burn-in
ni <- 20000 + nb #total iterations
nt <- 1  #thin
nc <- 3  #chains
adaptInterval = 100
maxContractions = 1000

#### COMPILE CONFIGURE AND BUILD ####
Rmodel <- nimbleModel(code = code, constants = const, data = dat, 
                      check = TRUE, calculate = TRUE, inits = inits)
Rmodel$simulate()
Rmodel$calculate()
#write_lines(Rmodel$getCode(), "AYNAipm_nimble.txt")


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
library(beepr)
beep(sound = 8)

#### RUN MCMC ####
t.start <- Sys.time()
#sink("somanyerrors.txt")
out <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = nc, inits = inits,
               setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE) 
#sink()
t.end <- Sys.time()
(runTime <- t.end - t.start)

# error.vec <- read_lines("somanyerrors.txt")
# error.vec <- error.vec[!(str_detect(error.vec, "initializing") & 
#                            #str_detect(error.vec, "IM\\[") & 
#                            str_detect(error.vec, "Inf") |
#                            str_detect(error.vec, "lifted") | str_detect(error.vec, "slice") )
# ] %>% unique() %>% sort()
# write_lines(error.vec, "somanyerrors.txt")

#### MAKE BEAUTIFUL PLOTS AND STUFF ####
pdf("survivalplots.pdf")
plot(out)
dev.off()

geldiag <- gelman.diag(out, multivariate=FALSE)
geldiag <- geldiag$psrf
View(geldiag)

summ <- summary(out) 
View(summ$statistics)

colnames(out$chain1)


#####  Steffen plotting code, modified by abby

library(tidyverse)
library(cowplot)
library(patchwork)
library(ggpomological)
library(rphylopic)
library(png)
library(magick)
library(mapproj)
library(sf)
library(rgdal)
library(here)
library(maps)
library(spData)
library(RColorBrewer)
library(bayesplot)
library(rgeos)
library(ggsn)
library(ggspatial)
library(tidybayes)
library(ggdist)
library(wesanderson)
library(strex)

pal <- wes_palette("Zissou1", 5, type = "continuous")
subadultcol <- pal[2]
adultcol <- pal[4]
age.pal <- c(subadultcol, adultcol)

plotdat <- rbind(out$chain1, out$chain2, out$chain3) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  as_tibble() %>% 
  filter(str_detect(name, "phi") & !str_detect(name, "mean") & !str_detect(name, "sigma")) %>% 
  mutate(Age = if_else(str_detect(name, "ad"), "Adult", "Juvenile")) %>% 
  mutate(Year = c(1985:2020)[str_first_number(name)])

png("survival_only_violin.png", width = 7, height = 5, units = "in", res = 300)
ggplot(plotdat, 
       mapping = aes(x = Year, y = value, color = Age, fill = Age)) +
  stat_eye(alpha = 0.7, .width = c(0.5, 0.95)) +
  theme_minimal() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text = element_text(size = 12),
        legend.position = "bottom", 
        plot.title.position = "plot",
        axis.title=element_text(size=12)) +
  scale_fill_manual(values = c(subadultcol, adultcol)) +
  scale_color_manual(values = c(subadultcol, adultcol)) +
  scale_y_continuous(expand = c(0, 0.05), 
                     limits = c(0, 1)) +  #no expansion below or above min/max
  xlab("Year") +
  ylab("Annual survival probability")
dev.off()

plotdat2 <- cbind(Count =rowSums(POP), Year = 1:length(rowSums(POP))) %>% as.data.frame()
plotdat2$Year <- 2007 + plotdat2$Year
plotdat2 <- left_join(plotdat2, PROD.DAT, by = "Year") %>% 
  mutate(coef = 4)

png("prod_and_counts.png", width = 7, height = 5, units = "in", res = 300)
ggplot(plotdat2) +
  geom_point(mapping = aes(Year, y = Count), alpha = 0.9, size = 2, color = pal[5]) +
  geom_smooth(mapping = aes(Year, y = Count),
              #method = "lm", 
              color = pal[5], fill = pal[5]) +
  geom_point(mapping = aes(Year, y = J*coef), alpha = 0.9, size = 2, color = pal[1]) +
  geom_smooth(mapping = aes(Year, y = J*coef), 
              #method = "lm", 
              color = pal[1], fill = pal[1]) +
  theme_minimal() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text = element_text(size = 12),
        legend.position = "bottom", 
        plot.title.position = "plot",
        axis.title=element_text(size=12)) +
  scale_y_continuous(
    
    # Features of the first axis
    name = "Population Counts",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./4, name="Chick Counts")
  ) +
  xlab("Year")
dev.off()

survival_posteriors <- out$chain1[, str_detect(colnames(out$chain1),"phi.ad\\[")]
ggplot() +
  geom_point(data = as.data.frame(survival_posteriors[1:36, ]), 
             aes(x=1985:2020,y=apply(survival_posteriors,2, median)), 
             size=2, color='darkred')+
  geom_smooth(method='lm') +
  xlab("Year") +
  ylab("Annual adult survival probability") +
  ylim(c(0, 1)) +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=14, color="black"), 
        axis.title=element_text(size=16), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank()) +
  geom_point(data = as.data.frame(survival_posteriors.juvs[1:36, ]),
             aes(x=1985:2020,y=apply(survival_posteriors.juvs,2, median)), 
             size=2, color='green') 
  

