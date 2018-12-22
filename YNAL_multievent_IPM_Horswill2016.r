##########################################################################
#
# YELLOW-NOSED ALBATROSS MULTI-STATE INTEGRATED POPULATION MODEL
#
##########################################################################
## originally written by Steffen Oppel in August 2013
## multievent survival model contributed by Sarah Converse September 2015
## initial development in YNAL_multistate_IPM.r, migrated to here on 15 Sept 2015
## added fecundity data into model to accommodate annual variation in fecundity on 19 Sept 2015
## MODIFIED BY CAT HORSWILL on 22 MARCH 2016
## updated 29 March 2016 to include lambda.t and rec_prop
## updated 4 April to include data from 2015/2016
## 2015 was an OUT year!?
# UPDATED 1 April 2016: SO: used mean B - possible extension would be to use status-specific B if we split successful/failed and nonbreeders
# UPDATED 4 April 2016: CH: 1.  increased precision on the observation model for popualiton size and fecundity following conversation with AB about census data
#                           2.  Change assignment of priors to deterministic for initial number of breeding pairs and fledglings in intial popualiton structure following conversation with AB about census data
#                           3.  Add in parameter to estimate annual breeding propensity based on population size and mean(B[])
#                           4.  Change population structure to 14 age classes following conversation between SC and PR; survival model was replaced with code from SC for same structure.      
#                           5.  Removal of emigration term

# UPDATED 8 April 2016: re-inserted emigration term but set prior to dunif(-10,10) to allow both emi- and immigration
# divided nB by 2 to get to account for sex-ratio

# UPDATED 15 SEPT 2016 to incorporate ICCAT FISHING EFFORT and split survival into two periods (<2006, 2006-2014)
# UPDATED 15 SEPT 2016 to extend projection into future - not yet functional
# NEEDED: Consider 1995 as outlier year with very low effort! (Rich Cuthbert, phone 15 Sept 2016)  


# Load necessary library
library(jagsUI)


#############################################################
#
# LOAD AND MANIPULATE DATA
#
##############################################################

#read in COUNT DATA 
setwd("G:\\STEFFEN\\RSPB\\UKOT\\Gough\\Raw_data")
setwd("S:\\ConSci\\DptShare\\steffenoppel\\RSPB\\UKOT\\Gough\\Raw_data")

#Number of breeding pairs and chicks fledged from census data
#extra year added at start to enable recruitment from t-1 in the population model
counts <- read.csv("YNAL_counts.csv")
nP<-c(NA, counts[,2])
nB<-c(NA, counts[,3])
nB<-as.integer(nB/2)		### this is necessary because the model does not include a 0.5 term to split the juveniles into males and females
#nP<-c(NA, 64, 38, 38, 59, 70, 45, 53, 53, 52, 38, 47, 46, 27, 29, 35, 45, 52, 45, 45, 55, 47, 60, 73, 55, 57, 69, 65, 63, 47, 31, 40, 44, 52,45)
#nB<-c(NA, 20, 12, 11, 16, 26, 15, 18, 18, 20, 12, 18, 16, 10, NA, 14, 16, 18,  8, 12, 17,  8, 23, 22, 18, 17, 20, 12, 22, 14,  8, NA, NA, NA,NA)
Fm<-nB/nP  #half measured fecundity to model only females



######### READ IN FISHING EFFORT ############################

fisheff<-read.table("AYNA_fish_effort_covariate.csv", header=T, sep=",")
f.eff<-fisheff$EFF_stand



######### READ IN ENCOUNTER HISTORY ############################

setwd("G:\\STEFFEN\\RSPB\\UKOT\\Gough\\AYNA_IPM")
setwd("S:\\ConSci\\DptShare\\steffenoppel\\RSPB\\UKOT\\Gough\\AYNA_IPM")

eh <- read.csv("YNdata_thru15_working AB 20160404.csv")
eh <- as.matrix(eh)
eh <- eh[,-(1:6)]

YNstate <- eh
colnames(YNstate) = colnames(eh)

# Events are
#0 = unobserved                        -> #6  Unobs
#1 = loaf in colony,                   -> #2  Loaf In
#2 = loaf out of colony,               -> #6  Unobs
#3 = successful breed in colony,       -> #3  Breed-in-S
#4 = breed outside colony,             -> #5  Breed-Out
#5 = failed breed in colony,           -> #4  Breed-in-F
#7 = confirmed dead,                   -> #6  Unobs
#8 = juvenile hatched in colony,       -> #1  Juv
#9 = juvenile hatched out of colony    -> #6  Unobs

#Condition on being captured first in the study area as a breeder or a juvenile
condcap <- rep(NA,nrow(YNstate))
for(i in 1:nrow(YNstate)){
  condcap[i] <- min(c(which(YNstate[i,] == 3),which(YNstate[i,] == 5),which(YNstate[i,] == 8),(ncol(YNstate)+1)))
  if(condcap[i] == '1'){
    YNstate[i,] <- YNstate[i,]
  }else YNstate[i,(1:(condcap[i]-1))] <- 0 
}    

#Reassign events as
#1 = Juvenile  
#2 = Loaf In
#3 = Breed In S 
#4 = Breed In F
#5 = Breed Out
#6 = no observation

#reassign events according to matrices in model        
for(i in 1:nrow(YNstate)){
  for(j in 1:ncol(YNstate)){
    if(YNstate[i,j] == '0'){
      YNstate[i,j] <- '6'
    }else if(YNstate[i,j] == '1'){
      YNstate[i,j] <- '2'
    }else if(YNstate[i,j] == '2'){
      YNstate[i,j] <- '6'
    }else if(YNstate[i,j] == '3'){
      YNstate[i,j] <- '3'
    }else if(YNstate[i,j] == '4'){
      YNstate[i,j] <- '5'
    }else if(YNstate[i,j] == '5'){
      YNstate[i,j] <- '4'
    }else if(YNstate[i,j] == '7'){
      YNstate[i,j] <- '6'
    }else if(YNstate[i,j] == '8'){
      YNstate[i,j] <- '1'
    }else if(YNstate[i,j] == '9'){
      YNstate[i,j] <- '6'
    }
  }
}

class(YNstate) <- 'numeric'


#pull out birds never observed in the admissable events
admit <- function(x) length(which(x!=6))
not.admit <- apply(YNstate,1,admit)
Y <- YNstate[-c(which(not.admit==0)),]

#number of individuals and number of years
nind <- nrow(Y)
nyear <- ncol(Y)

#determine first capture occassion for bounding likelihood
get.first <- function(x) min(which(x<6))
first <- apply(Y,1,get.first)

#Determine event at first release for all birds (should be 1, 2 or 3) 
first.event <- rep(NA,nrow(Y))
for(i in 1:nrow(Y)){
  first.event[i] <- Y[i,min(which(Y[i,]<6))]
}
#get birds that were captured as juveniles
juv <- Y[which(first.event==1),]
#get when first bred
first.breed <- rep(NA,nrow(juv))
for(i in 1:nrow(juv)){
  first.breed[i] <- min(which(juv[i,]==3 | juv[i,]==4 | juv[i,] ==5),(nyear+1))
}
#get those indviduals that recruited
juv.rec <- juv[which(first.breed <(nyear+1)),]
first.breed <- first.cap <- rep(NA,nrow(juv.rec))
#get the age at observed recruitment for those individuals that recruited
for(i in 1:nrow(juv.rec)){
  first.cap[i] <- which(juv.rec[i,]==1)
  first.breed[i] <- min(which(juv.rec[i,]==3 | juv.rec[i,]==4 | juv.rec[i,] ==5),(nyear+1))
}
recruit.obs <- first.breed-first.cap

#age for known age birds 
age <- matrix(data=NA,nrow=nrow(Y),ncol=ncol(Y))
condcap3 <- rep(NA,nrow(Y))
for(i in 1:nrow(Y)){
  condcap3[i] <- min(c(which(Y[i,] == 1),(ncol(Y)+1)))
  for(j in 1:ncol(age)){
    if(j < condcap3[i]){
      age[i,j] <- 'NA'
    }else if(j == condcap3[i]){
      age[i,j] <- 'NA'
    }else age[i,j] <- j-condcap3[i]
  }
}

class(age) <- 'numeric'

age[is.na(age)]<- nyear			### changed from 33


#Deal with variable years for captures 

#ps which are 0 out of colony
yrs.out <- c(18,22,23,24,25,26,27,32,34)
yrs.in <- c(1:(nyear-1)); yrs.in <- yrs.in[-yrs.out]
nyear.out <- length(yrs.out)
nyear.in <- length(yrs.in)

for(i in 1:nrow(Y)){
  for(j in c(1,yrs.in+1)){
    if(Y[i,j] == 5){
      Y[i,j] <- 6
    }
  }
}

#give this as data to initialize the first time step
z.first <- rep(NA,nind)
for(i in 1:nind){
  if(first.event[i]==1){
    z.first[i] <- 1
  }else if(first.event[i]==3){
    z.first[i] <- 3
  }else if(first.event[i]==4){
    z.first[i] <- 4
  }  
}

#Initialize process matrix
z.start <- Y
for(i in 1:nind){
  if(first[i]>1){
    z.start[i,1:(first[i])] <- NA
  }
}
first.breed <- rep(NA,nind)
for(i in 1:nind){
  first.breed[i] <- min(c(which(Y[i,] == 3),which(Y[i,] == 4),which(Y[i,] == 5)),(nyear+1)) 
}
for(i in 1:nind){
  for(t in first[i]:nyear){
    if(Y[i,t] == 1){
      z.start[i,t] <- NA
    }else if (Y[i,t] == 2){
      if(first.breed[i]<t){
        z.start[i,t] <- 5
      }else z.start[i,t] <- 2
    }else if (Y[i,t] == 3){
      z.start[i,t] <- NA
    }else if (Y[i,t] == 4){
      z.start[i,t] <- NA
    }else if (Y[i,t] == 5){
      z.start[i,t] <- 4
    }else if (Y[i,t] == 6){
      if(first.breed[i]<t){
        z.start[i,t] <- 4
      }else z.start[i,t] <- 2 
    }
  }
}
for(i in 1:nind){
  for(t in first[i]:nyear){
    if(t < first.breed[i]){
      if(age[i,t] > 14 & age[i,t]<nyear){
        z.start[i,t] <- 4
      }
    }
  }
}
for(i in 1:nind){
  z.start[i,first[i]] <- NA
}

z.data <- z.start
for(i in 1:nind){
  for(t in first[i]:nyear){
    if(Y[i,t] == 1){
      z.data[i,t] <- 1
    }else if (Y[i,t] == 2){
      z.data[i,t] <- NA
    }else if (Y[i,t] == 3){
      z.data[i,t] <- 3
    }else if (Y[i,t] == 4){
      z.data[i,t] <- 4
    }else z.data[i,t] <- NA 
  }
}
for(i in 1:nind){
  z.data[i,first[i]] <- NA
}


S.rand.st <- matrix(runif(2*(nyear-1)),nrow=2,ncol=nyear-1)




######### CREATE TWO PERIODS FOR SURVIVAL ####################

surv.period<-c(rep(1,24),rep(2,9))







#############################################################
#
# SPECIFY THE MODEL
#
##############################################################


# Specify model in BUGS language
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\AYNA_IPM")
setwd("S:\\ConSci\\DptShare\\steffenoppel\\RSPB\\UKOT\\Gough\\AYNA_IPM")
setwd("C:\\Users\\steffenoppel\\Desktop")

sink("YNAL.ipm.v10.txt")
cat("

    model {
    
    
    #observed 
    #            |--------------------------------------- Observed event ---------------------------------------|
    #true state     Juv   Loaf-In   Breed-In-S  Breed-In-F  Breed-Out  Unobs;
    #Juvenile      
    #Pre-breed
    #Breed-S
    #Breed-F
    #Skip
    #Dead
    
    
    #OBSERVATION MATRIX
    
    for(t in 1:(nyear-1)){
    
    pi[1,t,1]<-0;       pi[1,t,2]<-0;              pi[1,t,3]<-0;              pi[1,t,4]<-0;             pi[1,t,5]<-0;                  pi[1,t,6]<-0;
    
    pi[2,t,1]<-0;       pi[2,t,2]<-p[1,1,t];       pi[2,t,3]<-0;              pi[2,t,4]<-0;             pi[2,t,5]<-0;                  pi[2,t,6]<-(1-p[1,1,t]);   
    
    pi[3,t,1]<-0;       pi[3,t,2]<-0;              pi[3,t,3]<-p[2,1,t];       pi[3,t,4]<-0;             pi[3,t,5]<-p[2,2,t];           pi[3,t,6]<-(1-p[2,1,t]-p[2,2,t]); 
    
    pi[4,t,1]<-0;       pi[4,t,2]<-p[3,1,t];       pi[4,t,3]<-0;              pi[4,t,4]<-p[3,2,t];      pi[4,t,5]<-p[3,3,t];           pi[4,t,6]<-(1-p[3,1,t]-p[3,2,t]-p[3,3,t]);    
    
    pi[5,t,1]<-0;       pi[5,t,2]<-p[4,1,t];       pi[5,t,3]<-0;              pi[5,t,4]<-0;             pi[5,t,5]<-0;                  pi[5,t,6]<-(1-p[4,1,t]);
    
    pi[6,t,1]<-0;       pi[6,t,2]<-0;              pi[6,t,3]<-0;              pi[6,t,4]<-0;             pi[6,t,5]<-0;                  pi[6,t,6]<-1;
    
    
    p[1,2,t] <- 0;
    p[1,3,t] <- 0;
    p[2,3,t] <- 0;
    p[4,2,t] <- 0;
    p[4,3,t] <- 0;
    
    }
    
    #OBSERVATION MODELS AND PRIORS 
    
    for(t in 1:(nyear-1)){
    p[1,1,t] <- 1/(1+exp(-p.link[1,t])) 
    p.link[1,t] ~ dnorm(int.p[1],tau.p[1])
    
    p[4,1,t] <- 1/(1+exp(-p.link[2,t])) 
    p.link[2,t] ~ dnorm(int.p[2],tau.p[2])
    
    }
    
    for(t in 1:nyear.in){
    
    p[2,1,yrs.in[t]] <- 1/(1+exp(-p.link[3,t])) 
    p.link[3,t] ~ dnorm(int.p[3],tau.p[3])
    
    p[2,2,yrs.in[t]] <- 0
    
    
    p[3,1,yrs.in[t]] <- exp(p.link[4,t])/(1+exp(p.link[4,t])+exp(p.link[5,t]))
    p.link[4,t] ~ dnorm(int.p[4],tau.p[4])
    
    p[3,2,yrs.in[t]] <- exp(p.link[5,t])/(1+exp(p.link[4,t])+exp(p.link[5,t]))
    p.link[5,t] ~ dnorm(int.p[5],tau.p[5])
    
    p[3,3,yrs.in[t]] <- 0
    
    }
    
    for(t in 1:nyear.out){
    
    p[2,1,yrs.out[t]] <- exp(p.link[6,t])/(1+exp(p.link[6,t])+exp(p.link[7,t]))
    p.link[6,t] ~ dnorm(int.p[6],tau.p[6])
    
    p[2,2,yrs.out[t]] <- exp(p.link[7,t])/(1+exp(p.link[6,t])+exp(p.link[7,t]))
    p.link[7,t] ~ dnorm(int.p[7],tau.p[7])
    
    
    p[3,1,yrs.out[t]] <- exp(p.link[8,t])/(1+exp(p.link[8,t])+exp(p.link[9,t])+exp(p.link[10,t]))
    p.link[8,t] ~ dnorm(int.p[8],tau.p[8])
    
    p[3,2,yrs.out[t]] <- exp(p.link[9,t])/(1+exp(p.link[8,t])+exp(p.link[9,t])+exp(p.link[10,t]))
    p.link[9,t] ~ dnorm(int.p[9],tau.p[9])
    
    p[3,3,yrs.out[t]] <- exp(p.link[10,t])/(1+exp(p.link[8,t])+exp(p.link[9,t])+exp(p.link[10,t]))
    p.link[10,t] ~ dnorm(int.p[10],tau.p[10])
    
    }
    
    for(g in 1:10){
    int.p[g] ~ dunif(-5,5)
    
    tau.p[g] <- pow(sigma.p[g],-2)
    sigma.p[g] ~ dunif(0,10)
    }
    
    #STATE TRANSITION MATRIX
    
    #S = survive
    #R = recruit (breed for first time)
    #B = breed (breed again)
    #F = fledge
    
    for(i in 1:nind){
    for(t in 1:(nyear-1)){
    
    S.t[1,t,i,1]<-0;  S.t[1,t,i,2]<-S[1,t];             S.t[1,t,i,3]<-0;                        S.t[1,t,i,4]<-0;                           S.t[1,t,i,5]<-0;                S.t[1,t,i,6]<-1-S[1,t];         
    
    S.t[2,t,i,1]<-0;  S.t[2,t,i,2]<-S[1,t]*(1-R[i,t]);  S.t[2,t,i,3]<-S[1,t]*R[i,t]*F[t];       S.t[2,t,i,4]<-S[1,t]*R[i,t]*(1-F[t]);      S.t[2,t,i,5]<-0;                S.t[2,t,i,6]<-1-S[1,t];  
    
    S.t[3,t,i,1]<-0;  S.t[3,t,i,2]<-0;                  S.t[3,t,i,3]<-S[2,t]*B[1]*F[t];         S.t[3,t,i,4]<-S[2,t]*B[1]*(1-F[t]);        S.t[3,t,i,5]<-S[2,t]*(1-B[1]);  S.t[3,t,i,6]<-1-S[2,t];        
    
    S.t[4,t,i,1]<-0;  S.t[4,t,i,2]<-0;                  S.t[4,t,i,3]<-S[2,t]*B[2]*F[t];         S.t[4,t,i,4]<-S[2,t]*B[2]*(1-F[t]);        S.t[4,t,i,5]<-S[2,t]*(1-B[2]);  S.t[4,t,i,6]<-1-S[2,t];         
    
    S.t[5,t,i,1]<-0;  S.t[5,t,i,2]<-0;                  S.t[5,t,i,3]<-S[2,t]*B[3]*F[t];         S.t[5,t,i,4]<-S[2,t]*B[3]*(1-F[t]);        S.t[5,t,i,5]<-S[2,t]*(1-B[3]);  S.t[5,t,i,6]<-1-S[2,t];        
    
    S.t[6,t,i,1]<-0;  S.t[6,t,i,2]<-0;                  S.t[6,t,i,3]<-0;                        S.t[6,t,i,4]<-0;                           S.t[6,t,i,5]<-0;                S.t[6,t,i,6]<-1;        
    
    }
    }
    
    #PROCESS MODELS AND PRIORS
    beta.fish~dnorm(0,0.001)
    surv.intcpt~dnorm(0,0.001)    
    for(t in 1:(nyear-1)){
    for(m in 1:2){
    #S[m,t] <- 1/(1+exp(-S.rand[m,t]))
    logit(S[m,t]) <- surv.intcpt + beta.fish*f.eff[t] + S.rand[m,t]				### updated 15 Sept: included fishing effort effect on survival
    S.rand[m,t] ~ dnorm(int.S[m,surv.period[t]],tau.S[m,surv.period[t]])		### updated 15 Sept: included period-specific survival
    }
    }
    for(m in 1:2){
	for(per in 1:2){						### updated 15 Sept: included period-specific survival
    int.S[m,per] ~ dunif(-5,5)				### 
    tau.S[m,per] <- pow(sigma.S[m,per],-2)
    sigma.S[m,per] ~ dunif(0,10)
	}
    }
    
    for(i in 1:nind){
    for(t in 1:(nyear-1)){
    R[i,t] <- R.age[age[i,t]] 
    }
    }
    R.age[1] <- 0
    R.age[2] <- 0
    R.age[3] <- 0
    for(g in 4:14){
    R.age[g] <- 1/(1+exp(-R.rand[g]))
    R.rand[g] ~ dnorm(int.R,tau.R)
    }
    for(g in 15:34){
    R.age[g] <- 0
    }
    
    int.R ~ dunif(-5,5)
    tau.R <- pow(sigma.R,-2)
    sigma.R ~ dunif(0,10)
    
    for(m in 1:3){
    B[m] ~ dunif(0,1)
    }
    
    for(t in 1:(nyear-1)){
    F[t] <- 1/(1+exp(-F.rand[t]))
    F.rand[t] ~ dnorm(int.F,tau.F)
    }
    int.F ~ dunif(-5,5)
    tau.F <- pow(sigma.F,-2)
    sigma.F ~ dunif(0,10)
    
    
    #LIKELIHOOD
    
    for(i in 1:nind){
    z[i,first[i]] <- z.first[i]
    for(t in (first[i]+1):nyear){
    #state equation
    z[i,t] ~ dcat(S.t[z[i,(t-1)],t-1,i,])
    #observation equation
    Y[i,t] ~ dcat(pi[z[i,t],t-1,])
    }
    }
  
  ###################################################
  # 2 Likelihood for IPM
  ###################################################
  
  #############################
  # 2.1 System process
  #############################
  #six states #1. Juvenile 2. Pre-breed  3. Breed-S  4. Breed-F  5. Skip 6. Dead
  for (r in 2:nyear)   {			#### changed from nyear to nproj to facilitate future calculations
    # Update rules for process model
    #population size
    n[1,r+1] <- ntotf[r]
    n[2,r+1] <- sS[1,r]
    n[3,r+1] <- sS[2,r]
    n[4,r+1] <- sS[3,r]
    n[5,r+1] <- sS[4,r] 
    n[6,r+1] <- sS[5,r] 
    n[7,r+1] <- sS[6,r] 
    n[8,r+1] <- sS[7,r] 
    n[9,r+1] <- sS[8,r] 
    n[10,r+1] <- sS[9,r] 
    n[11,r+1] <- sS[10,r]
    n[12,r+1] <- sS[11,r]
    n[13,r+1] <- sS[12,r]
    n[14,r+1] <- sS[13,r] + sS[14,r]
    
    # Total number of experienced breeders in each age class
    #this separates the birds that have recruited from those that havent; the n[] matrix includes both
    ex_br[1,r+1] <- 0
    ex_br[2,r+1] <- 0
    ex_br[3,r+1] <- 0
    ex_br[4,r+1] <- n_rec2[3,r]
    ex_br[5,r+1] <- ex_br2[4,r] + n_rec2[4,r] 
    ex_br[6,r+1] <- ex_br2[5,r] + n_rec2[5,r] 
    ex_br[7,r+1] <- ex_br2[6,r] + n_rec2[6,r] 
    ex_br[8,r+1] <- ex_br2[7,r] + n_rec2[7,r] 
    ex_br[9,r+1] <- ex_br2[8,r] + n_rec2[8,r] 
    ex_br[10,r+1] <- ex_br2[9,r] + n_rec2[9,r] 
    ex_br[11,r+1] <- ex_br2[10,r] + n_rec2[10,r]
    ex_br[12,r+1] <- ex_br2[11,r] + n_rec2[11,r]
    ex_br[13,r+1] <- ex_br2[12,r] + n_rec2[12,r]
    ex_br[14,r+1] <- ex_br2[13,r] + n_rec2[13,r] +  ex_br2[14,r] + n_rec2[14,r]
    
    #=================================================================
    #observation model for population counts
    prec[r]<-100		#0000 # spike precision because no error on census data
    nP[r]~dnorm(ntotmu[r],prec[r])  # data on number of females in breeding popn; i.e. no of pairs 
    #emm[r] ~ dpois() # allow the census data to inform the population transition matrix
    emm[r] ~ dunif(-10,10) # allow both emigration and immigration to match with the census data
    
    #observation model for fecundity
    prec2[r]<-100		#0000 # spike precision because no error on census data
    Fm[r]~dnorm(F2[r],prec2[r])  #data on census fecundity re-estimates fecundity from MR
    
    
    ################################################
    ##PROCESSS MODEL FOR SURVIVAL
    #demographic parameters (S1,S2,F) taken from the survival model are [t-1] to put on same time scale as
    #population model, which needs to be one year forward to allow recrutiment terms to be estimated
    ### juveniles   
    ndum[1,r]<-round(equals(n[1,r],0)+n[1,r]-equals(n[1,r],0)*n[1,r])     # Check for zero number of trials     
    sdum[1,r]<-S[1,r-1]-equals(n[1,r],0)*S[1,r-1]  # Check for zero number of trials
    sS[1,r] ~ dbin(sdum[1,r],ndum[1,r]) # Number of survivors from fledgling age class
    
    ### FROM AGE 2:14  
    for (i in 2:14) {
      # Background survival probability 
      ndum[i,r]<-round(equals(n[i,r],0)+n[i,r]-equals(n[i,r],0)*n[i,r]) # Check for zero number of trials
      sdum[i,r]<-S[2,r-1]-equals(n[i,r],0)*S[2,r-1]  # Check for zero number of trials
      sS[i,r] ~ dbin(sdum[i,r],ndum[i,r]) # Number of survivors from each age class
    }
    
    #################################################
    ##PROCESS MODEL FOR FECUNDITY
    #First time breeders that recruit in each age class (R.age[i])
    F2[r]<-F[r-1]/2 #half the productivity to get on same scale as breeding pairs; i.e. only modelling females - assuming 50:50 sex ratio
    
    for (i in 1:14)  {
      ndum2[i,r]<-round(equals(n[i,r-1],0)+n[i,r-1]-equals(n[i,r-1],0)*n[i,r-1])  # Check for zero number of trials
      rdum[i,r]<-R.age[i]-equals(n[i,r-1],0)*R.age[i]# Check for zero number of trials
      n_rec[i,r] ~ dbin(rdum[i,r],ndum2[i,r]) # Total recruits per age class
      n_rec2[i,r]<-round(n_rec[i,r]) # integer value for recruits
      b1dum[i,r]<-F2[r]-equals(n[i,r-1],0)*F2[r]  # Check for zero number of trials
      Bi[i,r] ~ dbin(b1dum[i,r],n_rec2[i,r]) # Total births per recruiting age class
      
      #Experienced breeders
      #this is slightly underestimating the transition from failed to breeding; (B[1]) only
      # UPDATED 1 April 2016: used mean B - possible extension would be to use status-specific B if we split successful/failed and nonbreeders
      #**********************************************************
      # Background number of breeders
      ndum3[i,r]<-round(equals(ex_br[i,r],0)+ex_br[i,r]-equals(ex_br[i,r],0)*ex_br[i,r]) # Check for zero number of trials  
      b2dum[i,r]<-mean(B[])-equals(ex_br[i,r],0)*mean(B[])  # Check for zero number of trials, CHANGED B[1] TO mean(B[])
      #how many actually breed
      ndum4[i,r]<-round(ndum3[i,r]*b2dum[i,r])
      #how many succesfully raised chicks
      b3dum[i,r]<-F2[r]-equals(ex_br[i,r],0)*F2[r]  # Check for zero number of trials
      Bi_e[i,r] ~ dbin(b3dum[i,r],ndum4[i,r]) # Total births per experienced age class

      # update rules for experienced breeders using survival probability
      ndum5[i,r]<-round(equals(ex_br[i,r],0)+ex_br[i,r]-equals(ex_br[i,r],0)*ex_br[i,r]) # Check for zero number of trials
      b4dum[i,r]<-S[2,r-1]-equals(ex_br[i,r],0)*S[2,r-1]  # Check for zero number of trials
      ex_br2[i,r] ~ dbin(b4dum[i,r],ndum5[i,r]) # Number of survivors from each age class
    }
    
    ntotf[r]<- sum(Bi_e[4:14,r])   			#total number of chick fledged per year
    ntotmu[r]<- sum(ndum4[4:14,r])-emm[r] 	#total number of breeding pairs per year
    
    ntotr[r]<- sum(n_rec2[4:14,r])  #total number of recruits per year
    ntotn[r]<- sum(n[4:14,r])  #total number individuals >4yr old per year
    ntotbr[r]<- sum(ndum4[4:14,r]) # total individuals that breed
    
    r_prop[r]<-  ntotr[r]/(ntotn[r-1]) #annual recruitment propensity; r-1 as recrutiment is drawn from previous year
    br_prop[r]<- ntotbr[r]/sum(ex_br[4:14,r]) # annual breeding propensity    

    #estimating population growth rate between years
    ndumlambda[r]<-round(equals(ntotmu[r-1],0)+ntotmu[r-1]-equals(ntotmu[r-1],0)*ntotmu[r-1]) # Check for zero number of trials
    lambda.t[r]<- ntotmu[r]/ndumlambda[r]   # CHANGED TO GROWTH RATE rather than difference in number of breeding pairs between years
    loglam.t[r]<-log(lambda.t[r])## for calculating geometric mean of overall population growth rate
  }
  
  #### PRIORS FOR POPULATION MODEL #########  
  #Dummy variables
  ### CHANGED ON 12 AUGUST TO INCLUDE JUV SURVIVAL UP TO 13 YEARS - would need to distinguish between breeder and non-breeder though
    #Survival
    for (i in 1:13) {
    ndum[i,1]<-round(n[i,1])     
    sdum[i,1]<-S[1,1]  
    sS[i,1] ~ dbin(sdum[i,1],ndum[i,1]) 
    }
    
    ndum[14,1]<-round(n[14,1]) 
    sdum[14,1]<-S[2,1]      
    sS[14,1] ~ dbin(sdum[14,1],ndum[14,1]) 
  
  #Fecundity
  for (i in 1:14) {
    ndum2[i,1]<-round(n[i,1])  # dummy variable
    rdum[i,1]<-R.age[i]   # dummy variable
    n_rec[i,1] ~ dbin(rdum[i,1],ndum2[i,1]) # dummy variable
    n_rec2[i,1]<-round(n_rec[i,1]) # dummy variable
    b1dum[i,1]<-F2[2]    # dummy variable
    Bi[i,1] ~ dbin(b1dum[i,1],n_rec2[i,1]) # dummy variable
    
    ndum3[i,1]<-round(ex_br[i,1])  # dummy variable
    b2dum[i,1]<-mean(B[])  # dummy variable, CHANGED B[1] TO mean(B[])
    ndum4[i,1]<-round(ndum3[i,1]*b2dum[i,1]) # dummy variable
    b3dum[i,1]<-F2[1]  # dummy variable
    Bi_e[i,1] ~ dbin(b3dum[i,1],ndum4[i,1]) # dummy variable
    
    ndum5[i,1]<-round(ex_br[i,1])
    b4dum[i,1]<-S[2,1]
    ex_br2[i,1]~ dbin(b4dum[i,1],ndum5[i,1]) 
  }
  
  F2[1]<-F[1]/2
  ntotmu[1]<-br_ind
  ntotf[1]<-round(nBirth)
  emm[1]~ dunif(-10,10)
  ntotr[1]~ dunif(0,10)
  ntotn[1]<- round(ntot*0.6029799) #total population size age >4 yrs at t=1
  ntotbr[1]<-br_ind
  ndumlambda[1]<-1
  lambda.t[1]~ dunif(0,1)
  loglam.t[1]<-log(lambda.t[1])   
  r_prop[1]~ dunif(0,1)
  br_prop[1]~ dunif(0,1)
  
  # Initial total number of experienced breeders in each age class; year 1 is dummy year
  br_ind<-64 #Number of breeding pairs at t=1
  ind<-br_ind/prop # true number of breeding pairs including those on sabbatical
  prop~dunif(0.7,1.0)
  
  for (i in 1:2) {
    ex_br[1,i]<-0
    ex_br[2,i]<-0
    ex_br[3,i]<-0
    ex_br[4,i]<-0
    ex_br[5,i]<-round(ind*0.07142300)
    ex_br[6,i]<-round(ind*0.02433816)
    ex_br[7,i]<-round(ind*0.03263166)
    ex_br[8,i]<-round(ind*0.04375125)
    ex_br[9,i]<-round(ind*0.05865996)
    ex_br[10,i]<-round(ind*0.07864898)
    ex_br[11,i]<-round(ind*0.10544947)
    ex_br[12,i]<-round(ind*0.14138253)
    ex_br[13,i]<-round(ind*0.18956016)
    ex_br[14,i]<-round(ind*0.25415484)
  }
  
  
  #Initial population structure using stable age structure that starts breeding at age 4; 
  #calculated in R
  nBirth<-20 # Number of fledglings at at t=1
  ntot <- nBirth/0.17574058 #give estimate of total population size    
  #Year 1 is dummy year
  for (i in 1:2) {
    n[1,i]<-round(nBirth)
    n[2,i]<-round(ntot*0.11926872)
    n[3,i]<-round(ntot*0.10201077)
    n[4,i]<-round(ntot*0.08725001)
    n[5,i]<-round(ntot*0.07462511)
    n[6,i]<-round(ntot*0.06382700)
    n[7,i]<-round(ntot*0.05459136)
    n[8,i]<-round(ntot*0.04669210)
    n[9,i]<-round(ntot*0.03993584)
    n[10,i]<-round(ntot*0.03415720)
    n[11,i]<-round(ntot*0.02921472)
    n[12,i]<-round(ntot*0.02498741)
    n[13,i]<-round(ntot*0.02137178)
    n[14,i]<-round(ntot*0.12632740)
  }
  #### DERIVED PARAMETER: OVERALL POPULATION GROWTH RATE  #########
  #sumlam<-sum(lambda.t[3:nyear])
  #geometrate<-(1/(nyear-2))
  #mean.lambda<-pow(sumlam,geometrate)## geometric mean growth rate
  mean.lambda<-exp((1/(nyear-2))*sum(loglam.t[3:nyear]))   # Geometric mean 
  
  }  # End Model
",fill=TRUE)
sink()






#####################################################################################################################################################
#
# SET UP SIMULATION AND RUN MODEL
# 
#####################################################################################################################################################

######### CREATE PROJECTION PERIOD FOR FUTURE PVA ####################

#nproj<- nyear+25
#nP<-c(nP,rep(NA,25))
#Fm<-c(Fm,rep(NA,25))


# DEFINE SPECIFICATIONS FOR RUNNING THE MODEL

# 1. MCMC specification
ni <- 5000
nt <- 1
nb <- 1500
nc <- 4



# 2. Combine all data into a list and specify parameters


dataset <- list (Y=Y,z=z.data,z.first=z.first,first=first,nind=nind,nyear=nyear,
                 age=age,nyear.in=nyear.in,yrs.in=yrs.in,nyear.out=nyear.out,
                 yrs.out=yrs.out,nP=nP, Fm=Fm, f.eff=f.eff, surv.period=surv.period)

parameters <- c("int.S","sigma.S","R.age","int.R","sigma.R","B","F","int.F",
                "sigma.F","n_rec2","Bi","ntotmu","ntotf","n","Bi_e","ex_br",
                "ntotr","ntotn","r_prop","lambda.t","F","S","B","mean.lambda","emm","beta.fish") 


# 3. Specify initialisation values

inits <- function(){                                                                                                           
  list (z=z.start,S.rand=S.rand.st,int.S=matrix(runif(4),ncol=2),sigma.S=matrix(runif(4),ncol=2),B=runif(3),F.rand=runif(nyear-1))
}                                                                             




############# RUN THE MODEL ################

beg.time <- Sys.time()

# Call JAGS from R 
ipm.YNAL  <- jags(data=dataset, inits=inits, parameters.to.save=parameters, model.file="C:\\Users\\steffenoppel\\Desktop\\YNAL.ipm.v10.txt", n.chains = nc, n.thin = nt,n.iter = ni, n.burnin = nb, parallel = T, n.cores=8)

Sys.time() - beg.time


setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\AYNA_IPM")
setwd("S:\\ConSci\\DptShare\\steffenoppel\\RSPB\\UKOT\\Gough\\AYNA_IPM")

#ipm.YNAL $summary
write.table(ipm.YNAL$summary,"YNAL.IPM10.output.csv", sep=",",row.names=T)

plot(nP[2:nyear-1],ipm.YNAL$mean$ntotmu[2:nyear-1])
plot(ipm.YNAL$mean$ntotmu[3:nyear-1])


plot(ipm.YNAL $mean$r_prop[2:nyear-1],ipm.YNAL $mean$lambda.t[2:nyear-1]) #recruitment
plot(ipm.YNAL $mean$S[1,1:nyear-1],ipm.YNAL $mean$lambda.t[2:nyear]) #juv survival
plot(ipm.YNAL $mean$S[2,1:nyear-1],ipm.YNAL $mean$lambda.t[2:nyear]) #ad survival
plot(ipm.YNAL $mean$F[1:nyear-1],ipm.YNAL $mean$lambda.t[2:nyear]) #productivity
plot(ipm.YNAL $mean$emm[2:nyear-1],ipm.YNAL $mean$lambda.t[2:nyear-1]) #emmigration

ipm.YNAL $mean$B





#####################################################################################################################################################
#
# INSPECT MODEL OUTPUT
# 
#####################################################################################################################################################
library(Hmisc)

print(ipm.YNAL , 3)



##### PLOT POPULATION TREND ######
pdf("YNAL_population_trend.pdf", width=16, height=8)

fitted<-lower<-upper<-numeric()
year<-c(1982:2015)
n.years<-length(year)
for (i in 1:(n.years)){
fitted[i]<-quantile(ipm.YNAL$sims.list$ntotmu[,i], 0.5)
lower[i]<-quantile(ipm.YNAL$sims.list$ntotmu[,i], 0.025)
upper[i]<-quantile(ipm.YNAL$sims.list$ntotmu[,i], 0.975)}


par(mar=c(4.5,5,2,0), oma=c(0,0,0,0),cex=1.2)
errbar(c(1982:2014),fitted[2:n.years], lower[2:n.years],upper[2:n.years],ylim=c(0,100), xlim=c(1982,2014), ylab="Yellow-nosed Albatross population", xlab="Year", las=1, type='p', pch=16, main="",frame=FALSE, axes=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
axis(2, at=seq(0,100,20),labels=T, las=1, mgp=c(3,0.5,0), cex.axis=1.3)
axis(1, at=c(1982:2015),labels=T, , cex.axis=1.3)
points(x=c(1982:2014),fitted[2:n.years], type='l', col='black', lwd=2, lty=1)
points(x=c(1981:2015),nP, type='p', col='red', pch=4, cex=1.3)
legend(x=1982,y=20,legend=c("raw counts","model estimates"), col=c("red","black"), pch=c(4,16), bty='n',cex=1.3)
dev.off()



##### PLOT SURVIVAL ESTIMATES OVER TIME ######
pdf("YNAL_survival_trend.pdf", width=16, height=8)

fitted<-lower<-upper<-numeric()
year<-c(1982:2014)
n.years<-length(year)
for (i in 1:(n.years)){
fitted[i]<-quantile(ipm.YNAL$sims.list$S[,2,i], 0.5)
lower[i]<-quantile(ipm.YNAL$sims.list$S[,2,i], 0.025)
upper[i]<-quantile(ipm.YNAL$sims.list$S[,2,i], 0.975)}

mean(fitted)
mean(lower)
mean(upper)

par(mar=c(4.5,5,2,0), oma=c(0,0,0,0),cex=1.2)
errbar(c(1982:2014)+0.1,fitted, lower,upper,ylim=c(0,1), xlim=c(1982,2014), ylab="Adult survival probability", xlab="Year", las=1, type='p', pch=16, main="",frame=FALSE, axes=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
axis(2, at=seq(0,1,0.2),labels=T, las=1, mgp=c(3,0.5,0), cex.axis=1.3)
axis(1, at=c(1982:2014),labels=T, , cex.axis=1.3)
points(x=c(1982:2014),fitted, type='l', col='black', lwd=2, lty=1)


par(new=T)


n.years<-length(year)
for (i in 2:(n.years)){
fitted[i]<-quantile(ipm.YNAL$sims.list$S[,1,i], 0.5)
lower[i]<-quantile(ipm.YNAL$sims.list$S[,1,i], 0.025)
upper[i]<-quantile(ipm.YNAL$sims.list$S[,1,i], 0.975)}

mean(fitted, na.rm=T)
mean(fitted)
mean(lower)
mean(upper)

par(mar=c(4.5,5,2,0), oma=c(0,0,0,0),cex=1.2)
errbar(c(1982:2014)-0.1,fitted, lower,upper,ylim=c(0,1), xlim=c(1982,2014), ylab="", xlab="", las=1, type='p', pch=1, main="",frame=FALSE, axes=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)

legend(x=1983,y=0.2,legend=c("adult survival","juvenile survival"), pch=c(16,1), bty='n',cex=1.3)
#legend(x=1983,y=0.2,legend=c("mean adult survival","mean juvenile survival"), col=c("red","green"), lty=c(1,1), bty='n',cex=1.3)

abline(v=2005.5, lty=2, col='red', lwd=2)
text(2010,0.1,"bycatch mitigation launched", col='red')
dev.off()




##### PLOT PRODUCTIVITY OVER TIME ######
pdf("YNAL_productivity_trend.pdf", width=16, height=8)

fitted<-lower<-upper<-numeric()
year<-c(1982:2014)
n.years<-length(year)
for (i in 1:(n.years)){
fitted[i]<-quantile(ipm.YNAL$sims.list$F[,i], 0.5)
lower[i]<-quantile(ipm.YNAL$sims.list$F[,i], 0.025)
upper[i]<-quantile(ipm.YNAL$sims.list$F[,i], 0.975)}

mean(fitted)
mean(lower)
mean(upper)

par(mar=c(4.5,5,2,0), oma=c(0,0,0,0),cex=1.2)
errbar(c(1982:2014),fitted, lower,upper,ylim=c(0,1), xlim=c(1982,2013), ylab="annual productivity", xlab="Year", las=1, type='p', pch=16, main="",frame=FALSE, axes=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
axis(2, at=seq(0,1,0.2),labels=T, las=1, mgp=c(3,0.5,0), cex.axis=1.3)
axis(1, at=c(1982:2014),labels=T, , cex.axis=1.3)
## mean survival
abline(h=mean(ipm.YNAL$mean$F), lty=1, lwd=1,col='red')
abline(h=quantile(ipm.YNAL$sims.list$F, 0.025), lty=2, lwd=1,col='red')
abline(h=quantile(ipm.YNAL$sims.list$F, 0.975), lty=2, lwd=1,col='red')

dev.off()




#######################################################################################################################
################## PLOTTING THE POP GROWTH RATE AGAINST THE DEMOGRAPHIC PARAMETERS ####################################
#######################################################################################################################
### NOTE THAT POP SIZES AND LAMDA HAVE A DUMMY YEAR PRIOR TO FIRST DATA - this is ZERO in the output
### to correlate lambda and survival the indices need to be offset by 1 (lamda[2] ~ S[1])



### calculating annual pop growth ###

year<-c(1981:2014)
n.years<-length(year)
l.fitted<-l.lower<-l.upper<-rep(1,n.years)
#for (i in 3:(n.years)){
#l.fitted[i]<-quantile(ipm.YNAL$sims.list$ntotmu[,i], 0.5)/quantile(ipm.YNAL$sims.list$ntotmu[,i-1], 0.5)
#l.lower[i]<-quantile(ipm.YNAL$sims.list$ntotmu[,i], 0.025)/quantile(ipm.YNAL$sims.list$ntotmu[,i-1], 0.025)
#l.upper[i]<-quantile(ipm.YNAL$sims.list$ntotmu[,i], 0.975)/quantile(ipm.YNAL$sims.list$ntotmu[,i-1], 0.975)}
rowsout<-substr(row.names(ipm.YNAL$summary),1,8)=="lambda.t"
l.fitted<-ipm.YNAL$mean$lambda.t
l.lower<-ipm.YNAL$summary[rowsout,3]
l.upper<-ipm.YNAL$summary[rowsout,7]



### calculating number of recruits per year across all age classes ###

re.fitted<-re.lower<-re.upper<-rep(0,n.years)
rowsout<-substr(row.names(ipm.YNAL$summary),1,6)=="r_prop"
re.fitted<-ipm.YNAL$mean$r_prop
re.lower<-ipm.YNAL$summary[rowsout,3]
re.upper<-ipm.YNAL$summary[rowsout,7]




### compiling output for other demographic parameters ###

ad.fitted<-ad.lower<-ad.upper<-ju.fitted<-ju.lower<-ju.upper<-pr.fitted<-pr.lower<-pr.upper<-em.fitted<-em.lower<-em.upper<-rep(1,n.years)
year<-c(1982:2014)
n.years<-length(year)
for (i in 1:(n.years)){

ad.fitted[i]<-quantile(ipm.YNAL$sims.list$S[,2,i], 0.5)
ad.lower[i]<-quantile(ipm.YNAL$sims.list$S[,2,i], 0.025)
ad.upper[i]<-quantile(ipm.YNAL$sims.list$S[,2,i], 0.975)

ju.fitted[i]<-quantile(ipm.YNAL$sims.list$S[,1,i], 0.5)
ju.lower[i]<-quantile(ipm.YNAL$sims.list$S[,1,i], 0.025)
ju.upper[i]<-quantile(ipm.YNAL$sims.list$S[,1,i], 0.975)

pr.fitted[i]<-quantile(ipm.YNAL$sims.list$F[,i], 0.5)
pr.lower[i]<-quantile(ipm.YNAL$sims.list$F[,i], 0.025)
pr.upper[i]<-quantile(ipm.YNAL$sims.list$F[,i], 0.975)}

for (i in 1:(n.years)+1){
em.fitted[i]<-quantile(ipm.YNAL$sims.list$emm[,i], 0.5)/(quantile(ipm.YNAL$sims.list$ntotmu[,i], 0.5)+quantile(ipm.YNAL$sims.list$emm[,i], 0.5))
em.lower[i]<-quantile(ipm.YNAL$sims.list$emm[,i], 0.025)/(quantile(ipm.YNAL$sims.list$ntotmu[,i], 0.025)+quantile(ipm.YNAL$sims.list$emm[,i], 0.025))
em.upper[i]<-quantile(ipm.YNAL$sims.list$emm[,i], 0.975)/(quantile(ipm.YNAL$sims.list$ntotmu[,i], 0.975)+quantile(ipm.YNAL$sims.list$emm[,i], 0.975))}



##### PLOTTING OUTPUT ####
## note that time series are shifted


pdf("YNAL_lambda_correlations.pdf", width=9, height=9)

par(mfrow=c(2,2))

plot(l.fitted[2:nyear]~ad.fitted[1:nyear-1], xlim=c(0.5,1), ylim=c(0.5,1.5), xlab="Adult survival probability",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(ad.lower[1:nyear-1],l.fitted[2:nyear],ad.upper[1:nyear-1],l.fitted[2:nyear] ,col="gray", lty=1, lwd=0.5)
segments(ad.fitted[1:nyear-1],l.lower[2:nyear],ad.fitted[1:nyear-1],l.upper[2:nyear] ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted[2:nyear],ad.fitted[1:nyear-1],alternative = c("two.sided"),method = "spearman")
text(0.5,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

plot(l.fitted[2:nyear]~ju.fitted[1:nyear-1], xlim=c(0.5,1), ylim=c(0.5,1.5), xlab="Immature survival probability",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(ju.lower[1:nyear-1],l.fitted[2:nyear],ju.upper[1:nyear-1],l.fitted[2:nyear] ,col="gray", lty=1, lwd=0.5)
segments(ju.fitted[1:nyear-1],l.lower[2:nyear],ju.fitted[1:nyear-1],l.upper[2:nyear] ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted[2:nyear],ju.fitted[1:nyear-1],alternative = c("two.sided"),method = "spearman")
text(0.5,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

plot(l.fitted[2:nyear]~pr.fitted[1:nyear-1], xlim=c(0,0.8), ylim=c(0.5,1.5), xlab="Annual productivity",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(pr.lower[1:nyear-1],l.fitted[2:nyear],pr.upper[1:nyear-1],l.fitted[2:nyear] ,col="gray", lty=1, lwd=0.5)
segments(pr.fitted[1:nyear-1],l.lower[2:nyear],pr.fitted[1:nyear-1],l.upper[2:nyear] ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted[2:nyear],pr.fitted[1:nyear-1],alternative = c("two.sided"),method = "spearman")
text(0,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

plot(l.fitted[2:(nyear-1)]~re.fitted[2:(nyear-1)], xlim=c(0,0.6), ylim=c(0.5,1.5), xlab="Annual recruitment",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(re.lower[2:(nyear-1)],l.fitted[2:(nyear-1)],re.upper[2:(nyear-1)],l.fitted[2:(nyear-1)] ,col="gray", lty=1, lwd=0.5)
segments(re.fitted[2:(nyear-1)],l.lower[2:(nyear-1)],re.fitted[2:(nyear-1)],l.upper[2:(nyear-1)] ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted[2:(nyear-1)],re.fitted[2:(nyear-1)],alternative = c("two.sided"),method = "spearman")
text(0,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

mtext("Atlantic Yellow-nosed Albatross demographic correlations", side=3, cex=1.9, outer=T, line=-2)
dev.off()



#### removed from final model on 4 April 2016

pdf("YNAL_lambda_emigration_correlation.pdf", width=6, height=6)

plot(l.fitted[2:(nyear-1)]~em.fitted[2:(nyear-1)], xlim=c(-0.5,0.5), ylim=c(0.5,1.5), xlab="Annual emigration",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(em.lower[2:(nyear-1)],l.fitted[2:(nyear-1)],em.upper[2:(nyear-1)],l.fitted[2:(nyear-1)] ,col="gray", lty=1, lwd=0.5)
segments(em.fitted[2:(nyear-1)],l.lower[2:(nyear-1)],em.fitted[2:(nyear-1)],l.upper[2:(nyear-1)] ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted[2:(nyear-1)],em.fitted[2:(nyear-1)],alternative = c("two.sided"),method = "spearman")
text(0,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

dev.off()





#######################################################################################################################
################## CALCULATING SENSITIVITY AND ELASTICITY ####################################
#######################################################################################################################

### DEFINITIONS
# sensitivity: delta(lambda)/delta(param) or
# elasticity: param/lambda*delta(lambda)/delta(param)

output<-data.frame()

sens.ad.surv<-numeric()
elas.ad.surv<-numeric()
for (y in 2:(nyear-1)){
sens.ad.surv[y-1]<-abs(l.fitted[y+1]-(l.fitted[y]))/abs(ad.fitted[y]-ad.fitted[y-1])
elas.ad.surv[y-1]<-(abs(l.fitted[y+1]-l.fitted[y])/abs(ad.fitted[y]-ad.fitted[y-1]))*(ad.fitted[y]/l.fitted[y+1])
}

output[1,1]<-"adult survival"
output[1,2]<-mean(sens.ad.surv)
output[1,3]<-mean(elas.ad.surv)




sens.juv.surv<-numeric()
elas.juv.surv<-numeric()
for (y in 2:(nyear-1)){
sens.juv.surv[y-1]<-abs(l.fitted[y+1]-(l.fitted[y]))/abs(ju.fitted[y]-ju.fitted[y-1])
elas.juv.surv[y-1]<-abs(l.fitted[y+1]-(l.fitted[y]))/abs(ju.fitted[y]-ju.fitted[y-1])*(ju.fitted[y]/l.fitted[y+1])
}

output[2,1]<-"juvenile survival"
output[2,2]<-mean(sens.juv.surv)
output[2,3]<-mean(elas.juv.surv)






sens.fec<-numeric()
elas.fec<-numeric()
for (y in 2:(nyear-1)){
sens.fec[y-1]<-abs(l.fitted[y+1]-(l.fitted[y]))/abs(pr.fitted[y]-pr.fitted[y-1])
elas.fec[y-1]<-(abs(l.fitted[y+1]-l.fitted[y])/abs(pr.fitted[y]-pr.fitted[y-1]))*(pr.fitted[y]/l.fitted[y+1])
}

output[3,1]<-"productivity"
output[3,2]<-mean(sens.fec)
output[3,3]<-mean(elas.fec)




sens.rec<-numeric()
elas.rec<-numeric()
for (y in 3:(n.years)){
sens.rec[y-2]<-abs(l.fitted[y]-(l.fitted[y-1]))/abs(re.fitted[y]-re.fitted[y-1])
elas.rec[y-2]<-(abs(l.fitted[y]-l.fitted[y-1])/abs(re.fitted[y]-re.fitted[y-1]))*(re.fitted[y]/l.fitted[y])
}

output[4,1]<-"recruitment"
output[4,2]<-mean(sens.rec)
output[4,3]<-mean(elas.rec)





sens.em<-numeric()
elas.em<-numeric()
for (y in 3:(n.years)){
sens.em[y-2]<-abs(l.fitted[y]-(l.fitted[y-1]))/abs(em.fitted[y]-em.fitted[y-1])
elas.em[y-2]<-(abs(l.fitted[y]-l.fitted[y-1])/abs(em.fitted[y]-em.fitted[y-1]))*(em.fitted[y]/l.fitted[y])
}
sens.em[7]<-0
elas.em[7]<-0
output[5,1]<-"emigration"
output[5,2]<-mean(sens.em)
output[5,3]<-mean(elas.em)

names(output)<-c("parameter","sensitivity","elasticity")
write.table(output, "YNAL_demographics_sensitivity.csv", sep=",", row.names=F)













####################### do with all simulations ?? [abandoned] ####

sens.ad.surv<-array()
dim(ipm.YNAL$sims.list$S)
for (s in 1:32000)
for (y in 2:(n.years)){
sens.ad.surv[s,y]<-(ipm.YNAL$sims.list$S[s,2,y]-ipm.YNAL$sims.list$S[s-1,2,y])/((ipm.YNAL$sims.list$ntotmu[s,y]/ipm.YNAL$sims.list$ntotmu[s-1,y])-(ipm.YNAL$sims.list$ntotmu[s-1,y]/ipm.YNAL$sims.list$ntotmu[s-2,y]))
}
hist(sens.ad.surv)





