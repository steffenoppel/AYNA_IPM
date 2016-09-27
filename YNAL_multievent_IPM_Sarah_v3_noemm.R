##########################################################################
#
# YELLOW-NOSED ALBATROSS MULTI-STATE INTEGRATED POPULATION MODEL
#
##########################################################################
## written August 2016 by Sarah Converse 
## On top of ME model from October 2015 
## Variation on model by Cat Horswill

# Load necessary libraries
library(jagsUI)


#############################################################
#
# LOAD AND MANIPULATE DATA
#
##############################################################

#read in COUNT DATA 
setwd("C:\\Users\\sconverse\\Documents\\Albatross\\Peter Ryan-YNAL\\Analysis\\IPM")

#Number of breeding pairs and chicks fledged from census data
############################################################################# missing counts in 2015 
counts <- read.csv("YNAL_counts.csv")
nP<-counts[,2]  #breeding pairs
nF<-counts[,3]   #young females (assuming 1:1 sex ratio) 

#Stable age distribution
stable.rate <- c(0.075,0.070,0.064,0.060,0.045,0.034,0.026,0.019,0.015,0.011,0.009,0.007,0.005,0.004,0.003,0.150,0.107,0.296)
stable <- rmultinom(1,450,stable.rate)

#read in MARK EH
eh <- read.csv("YNdata_thru15_working.csv")
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



#############################################################
#
# SPECIFY THE MODEL
#
##############################################################

sink("YNAL.ipm.SJC.txt")
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
  int.p[g] ~ dunif(-15,15)
  
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
    
    S.t[3,t,i,1]<-0;  S.t[3,t,i,2]<-0;                  S.t[3,t,i,3]<-S[2,t]*B[t]*F[t];         S.t[3,t,i,4]<-S[2,t]*B[t]*(1-F[t]);        S.t[3,t,i,5]<-S[2,t]*(1-B[t]);  S.t[3,t,i,6]<-1-S[2,t];        
    
    S.t[4,t,i,1]<-0;  S.t[4,t,i,2]<-0;                  S.t[4,t,i,3]<-S[2,t]*B[t]*F[t];         S.t[4,t,i,4]<-S[2,t]*B[t]*(1-F[t]);        S.t[4,t,i,5]<-S[2,t]*(1-B[t]);  S.t[4,t,i,6]<-1-S[2,t];         
    
    S.t[5,t,i,1]<-0;  S.t[5,t,i,2]<-0;                  S.t[5,t,i,3]<-S[2,t]*B[t]*F[t];         S.t[5,t,i,4]<-S[2,t]*B[t]*(1-F[t]);        S.t[5,t,i,5]<-S[2,t]*(1-B[t]);  S.t[5,t,i,6]<-1-S[2,t];        
    
    S.t[6,t,i,1]<-0;  S.t[6,t,i,2]<-0;                  S.t[6,t,i,3]<-0;                        S.t[6,t,i,4]<-0;                           S.t[6,t,i,5]<-0;                S.t[6,t,i,6]<-1;        
    
  }
}
    
#PROCESS MODELS AND PRIORS
    
for(t in 1:(nyear-1)){
  for(m in 1:2){
    S[m,t] <- 1/(1+exp(-S.rand[m,t]))
    S.rand[m,t] ~ dnorm(int.S[m],tau.S[m])
  }
}
for(m in 1:2){
  int.S[m] ~ dunif(-15,15)
  tau.S[m] <- pow(sigma.S[m],-2)
  sigma.S[m] ~ dunif(0,10)
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
    
int.R ~ dunif(-15,15)
tau.R <- pow(sigma.R,-2)
sigma.R ~ dunif(0,10)
    
for(t in 1:(nyear-1)){
  B[t] <- 1/(1+exp(-B.rand[t]))
  B.rand[t] ~ dnorm(int.B,tau.B)
}
int.B ~ dunif(-15,15)
tau.B <- pow(sigma.B,-2)
sigma.B  ~ dunif(0,10)

for(t in 1:(nyear-1)){
  F[t] <- 1/(1+exp(-F.rand[t]))
  F.rand[t] ~ dnorm(int.F,tau.F)
}
int.F ~ dunif(-15,15)
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
# Likelihood for IPM
###################################################

for(r in 1:nyear){
  nF[r] ~ dnorm(n[1,r],tau.Fl) 
  nP[r] ~ dnorm(n.breeders[r],tau.Pr)
}
tau.Fl <- 100000
tau.Pr <- pow(sigma.Pr,-2) 
sigma.Pr ~ dunif(0,3)

for(i in 1:18){
 n[i,1] <- stable[i]
}
n.breeders[1] <- n[16,1] + n[17,1]

for(r in 1:(nyear-1)){

  ### all age and breeding groups 
  n[1,r+1]  <- n[16,r+1]                         #total fledglings                 
  n[2,r+1]  <- n[1,r]*S[1,r]*0.5                 #1yr NB - females only  
  n[3,r+1]  <- n[2,r]*S[1,r]                     #2yr NB
  n[4,r+1]  <- n[3,r]*S[1,r]                     #3yr NB
  n[5,r+1]  <- n[4,r]*S[1,r]*(1-R.age[4])        #4yr NB
  n[6,r+1]  <- n[5,r]*S[1,r]*(1-R.age[5])        #5yr NB
  n[7,r+1]  <- n[6,r]*S[1,r]*(1-R.age[6])        #6yr NB
  n[8,r+1]  <- n[7,r]*S[1,r]*(1-R.age[7])        #7yr NB
  n[9,r+1]  <- n[8,r]*S[1,r]*(1-R.age[8])        #8yr NB
  n[10,r+1] <- n[9,r]*S[1,r]*(1-R.age[9])        #9yr NB
  n[11,r+1] <- n[10,r]*S[1,r]*(1-R.age[10])      #10yr NB
  n[12,r+1] <- n[11,r]*S[1,r]*(1-R.age[11])      #11yr NB
  n[13,r+1] <- n[12,r]*S[1,r]*(1-R.age[12])      #12yr NB
  n[14,r+1] <- n[13,r]*S[1,r]*(1-R.age[13])      #13yr NB
  n[15,r+1] <- n[14,r]*S[1,r]*(1-R.age[14])      #14yr NB

  n.breeders[r+1] <- n[4,r]*S[1,r]*R.age[4] +    #Breeders 
               n[5,r]*S[1,r]*R.age[5] +
               n[6,r]*S[1,r]*R.age[6] +
               n[7,r]*S[1,r]*R.age[7] +
               n[8,r]*S[1,r]*R.age[8] +
               n[9,r]*S[1,r]*R.age[9] +
               n[10,r]*S[1,r]*R.age[10] +
               n[11,r]*S[1,r]*R.age[11] +
               n[12,r]*S[1,r]*R.age[12] +
               n[13,r]*S[1,r]*R.age[13] +
               n[14,r]*S[1,r]*R.age[14] +
               n[15,r]*S[1,r]*R.age[15] +
               n[16,r]*S[2,r]*B[1] +
               n[17,r]*S[2,r]*B[1] +
               n[18,r]*S[2,r]*B[2] 
              
  n[16,r+1] <- n.breeders[r+1]*F[r] 
  n[17,r+1] <- n.breeders[r+1]*(1-F[r])

  n[18,r+1] <- n[16,r]*S[2,r]*(1-B[1]) +                #Skipped breeders
                n[17,r]*S[2,r]*(1-B[1]) +
                n[18,r]*S[2,r]*(1-B[2])
}

}  # End Model
",fill=TRUE)
sink()


#####################################################################################################################################################
#
# SET UP SIMULATION AND RUN MODEL
# 
#####################################################################################################################################################

# DEFINE SPECIFICATIONS FOR RUNNING THE MODEL

# 1. MCMC specification
ni <- 20000
nt <- 1
nb <- 15000
nc <- 3


# 2. Combine all data into a list and specify parameters
dataset <- list (Y=Y,z=z.data,z.first=z.first,first=first,nind=nind,nyear=nyear,
                 age=age,nyear.in=nyear.in,yrs.in=yrs.in,nyear.out=nyear.out,
                 yrs.out=yrs.out, nF=nF, nP=nP, stable=stable)

parameters <- c("int.S","sigma.S","R.age","int.R","sigma.R","int.B","sigma.B","int.F",
                "sigma.F","n.breeders","sigma.Pr") 


# 3. Specify initialisation values
S.rand.st <- matrix(runif(2*(nyear-1),2,3),nrow=2,ncol=nyear-1)

inits <- function(){                                                                                                           
  list (z=z.start,S.rand=S.rand.st,int.S=runif(2,2,4),sigma.S=runif(2),int.B=runif(1,0,2),sigma.B=runif(1),B.rand=rep(0.8,nyear-1),F.rand=runif(nyear-1),int.F=runif(1,1,3),sigma.F=runif(1),mean.emm=0)
}                                                                             

############# RUN THE MODEL ################

beg.time <- Sys.time()

# Call JAGS from R 

ipm.YNAL  <- jags(data=dataset, inits=inits, parameters.to.save=parameters, "YNAL.ipm.SJC.txt", n.chains = nc, n.thin = nt,n.iter = ni, n.burnin = nb, parallel = TRUE)

end.time <- Sys.time()
end.time - beg.time
