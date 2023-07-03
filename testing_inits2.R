load("IPM_AEB_dat.Rdata")
const <- c(const, maxAge = const$n.occasions)
maxage <- const$maxAge
n.years.fec <- const$n.years.fec

#### INITIAL VALUES ####
iann.fec <- rep(0.4, n.years.fec)
imean.p.propensity <- 0.5
imean.p.det <- 2.5
imu.p.juv <- -4
iagebeta <- 2
isigma.p <- 0.1

imean.phi.juv <- 0.7
imean.phi.ad <- 0.9
isigma.phi <- 0.1

iN.ad.surv <- rep(0, n.years.fec)
iN.breed.ready <- rep(0, n.years.fec)
inestlings <- rep(0, n.years.fec)
iN.recruits <- rep(NA, n.years.fec)
iNtot.breed <- rep(NA, n.years.fec)
iJUV <- rep(NA, n.years.fec)
iN.atsea <- rep(NA, n.years.fec)
IMinits <- array(NA, dim = c(n.years.fec, maxage, 3))

IMinits[1,1,1] = max(c(rnorm(1, 263*0.5, 20)), 1) %>% round() 
IMinits[1,1,2] <- 0
IMinits[1,1,3] <- round(IMinits[1,1,1]) - IMinits[1,1,2]

IMinits[1,2,1] = max(c(rnorm(1, 275*0.5, 20)), 1)%>% round()
IMinits[1,2,2] <- rbinom(1, IMinits[1, 2, 1], plogis(imu.p.juv + iagebeta*sqrt(2)))
IMinits[1,2,3] <- round(IMinits[1,2,1]) - IMinits[1,2,2]

IMinits[1,3,1] = max(c(rnorm(1, 264*0.5, 20)), 1)%>% round()
IMinits[1,3,2] <- rbinom(1, IMinits[1, 3, 1], plogis(imu.p.juv + iagebeta*sqrt(3)))
IMinits[1,3,3] <- round(IMinits[1,3,1]) - IMinits[1,3,2]

IMinits[1,4,1] = max(c(rnorm(1, 177*0.5, 20)), 1)%>% round()
IMinits[1,4,2] <- rbinom(1, IMinits[1, 4, 1], plogis(imu.p.juv + iagebeta*sqrt(4)))
IMinits[1,4,3] <- round(IMinits[1,4,1]) - IMinits[1,4,2]

IMinits[1,5,1] = max(c(rnorm(1, 290*0.5, 20)), 1)%>% round()
IMinits[1,5,2] <- rbinom(1, IMinits[1, 5, 1], plogis(imu.p.juv + iagebeta*sqrt(5)))
IMinits[1,5,3] <- round(IMinits[1,5,1]) - IMinits[1,5,2]

IMinits[1,6,1] = max(c(rnorm(1, 90*0.5, 20)), 1)%>% round()
IMinits[1,6,2] <- rbinom(1, IMinits[1, 6, 1], plogis(imu.p.juv + iagebeta*sqrt(6)))
IMinits[1,6,3] <- round(IMinits[1,6,1]) - IMinits[1,6,2]

IMinits[1,7,1] = max(c(rnorm(1, 158*0.5, 20)), 1)%>% round()
IMinits[1,7,2] <- rbinom(1, IMinits[1, 7, 1], plogis(imu.p.juv + iagebeta*sqrt(7)))
IMinits[1,7,3] <- round(IMinits[1,7,1]) - IMinits[1,7,2]

for (age in 8:maxage) {
  #IMinits[1, age, 1] = rbinom(1, IMinits[1, age-1, 3], pow(0.9,(age-1)))
  # AEB change
  IMinits[1, age, 1] = rbinom(1, IMinits[1, age-1, 3], imean.phi.ad)
  IMinits[1, age, 2] <- rbinom(1, IMinits[1, age, 1], plogis(imu.p.juv + iagebeta*sqrt(age)))
  # print(paste(age, IMinits[1,age,1] - IMinits[1,age,2], sep = " "))
  IMinits[1, age, 3] <- IMinits[1,age,1] - IMinits[1,age,2]
}

IMinits[1, , ]
IMinits[,, 1]

iN.recruits[1] <- sum(IMinits[1,1:maxage,2]) 
iNtot.breed[1] <- max(rnorm(1, 640*0.5,sd = 20), 1) %>% round()#change here

iJUV[1] <- max(rnorm(1, 232*0.5, sd = 20), 1) %>% round()
iN.atsea[1] <- max(rnorm(1, 224*0.5,sd = 20), 1) %>% round() #change here

for (tt in 2:n.years.fec) {
  iN.ad.surv[tt] <- rbinom(1, iNtot.breed[tt-1]+iN.atsea[tt-1], imean.phi.ad)
  iN.breed.ready[tt] <- rbinom(1, iN.ad.surv[tt], imean.p.propensity)
  
  IMinits[tt,1,1] <- rbinom(1, iJUV[tt-1], imean.phi.juv)
  IMinits[tt, 1, 2] <- 0
  IMinits[tt, 1, 3] <- IMinits[tt,1,1] - IMinits[tt,1,2]
  
  for (age in 2:maxage) {
    IMinits[tt, age, 1] <- rbinom(1, IMinits[tt-1, age-1, 3], imean.phi.ad)
    #IMinits[tt, age, 2] <- min(IMinits[tt, age, 1], IMinits[tt, age-1, 3]) * p.juv.recruit.inits[age, tt] %>% round()
    IMinits[tt, age, 2] <- rbinom(1, IMinits[tt, age, 1], plogis(imu.p.juv + iagebeta*sqrt(age)))
    IMinits[tt, age, 3] <- IMinits[tt,age,1] - IMinits[tt,age,2]
  }
  
  iN.recruits[tt] <- sum(IMinits[tt,1:maxage,2])
  iNtot.breed[tt]<- iN.breed.ready[tt]+iN.recruits[tt]
  inestlings[tt] <- (iann.fec[tt] * 0.5 * iNtot.breed[tt]) %>% round()
  iJUV[tt] <- rpois(1, inestlings[tt])
  
  iN.atsea[tt] <- iN.ad.surv[tt]-iN.breed.ready[tt]
}

IMinits[,, 1]
IMinits[1, , ]

# checking if the inits are reasonable - this should be close to 0, on average
(iNtot.breed * 2) - rowSums(dat$y.count[1:13, ])
mean((iNtot.breed * 2 ) - rowSums(dat$y.count[1:13, ]))
(iNtot.breed * 2)
rowSums(dat$y.count[1:13, ])

######## state space version #######

# dat ####
load("AYNA_IPM_input.marray.Rdata")

CH_all <- bind_rows(AYNA_AD, AYNA_CHICK)
first <- apply(CH_all[, -1], 1, function(x){min(which(x != 0))})

age <- matrix(2, nrow = dim(CH_all)[1], ncol = dim(CH_all[,-1])[2])
for (i in (dim(AYNA_AD)[1] + 1):(dim(CH_all)[1])) {
  age[i, first[i]] <- 1
}

obs.mat.init.JUV <- matrix(
  c(
    0,1, # when banded as juvenile, assume observed and state = alive at sea
    1,0,
    1,0,
    1,0,
    1,0,
    1,0,
    1,0,
    1,0
  ), nrow = 8, ncol = 2, byrow = T
)
obs.mat.init.AD <- matrix(
  c(
    1,0,
    1,0,
    1,0,
    1,0,
    1,0,
    0,1, # when banded as an adult, assume observed and state = breeding in colony
    1,0,
    1,0
  ), nrow = 8, ncol = 2, byrow = T
)

Zdat <- matrix(NA, nrow = dim(CH_all)[1], ncol = dim(CH_all[,-1])[2])
init <- matrix(NA, nrow = dim(CH_all)[1], ncol = 8)
obs.mat.init <- array(NA, dim = c(dim(CH_all)[1], 8, 2))
for (i in 1:nrow(Zdat)) {
  if (age[i, first[i]] == 1) { # if age is 1 -alive at sea never bred
    Zdat[i, first[i]] <- 1 
    init[i, 1:8] <- c(1,0,0,0,0,0,0,0)
    obs.mat.init[i, , ] <- obs.mat.init.JUV
  } else { # if age is 2 - breeding in colony
    Zdat[i, first[i]] <- 6 
    init[i, 1:8] <- c(0,0,0,0,0,1,0,0)
    obs.mat.init[i, , ] <- obs.mat.init.AD
  }
}

Zinits <- matrix(NA, nrow = dim(CH_all)[1], ncol = dim(CH_all[,-1])[2])
for (i in 1:nrow(Zinits)) {
  if (age[i, first[i]] == 1) {
    CH_curr <- CH_all[i, -1] %>% as.numeric() 
    ever_observed <- any(CH_curr[(first[i]+1):(dim(CH_all[,-1])[2])] == 1)
    if (ever_observed) {
      all_inds <- (which(CH_curr[(1):(dim(CH_all[,-1])[2])] == 1))
      inds_2plus <- all_inds[2:length(all_inds)]
      Zinits[i, first[i]:(inds_2plus[1]-1)] <- 1
      Zinits[i, inds_2plus[1]-1] <- 5
      Zinits[i, min(inds_2plus[1]+1, dim(CH_all[,-1])[2]):(dim(CH_all[,-1])[2])] <- 7
      Zinits[i, inds_2plus] <- 6
    } else {
      for (t in (first[i]+1):(dim(CH_all[,-1])[2])) {
        curr_age <- t - first[i]
        if (curr_age < 10) {
          Zinits[i, t] <- 1
        } else if (curr_age == 10) {
          Zinits[i, t] <- 5
        } else if (curr_age > 10) {
          Zinits[i, t] <- 7
        }
      } 
    }
    Zinits[i, first[i]] <- NA
  } else { # if age is 2 - breeding in colony
    CH_curr <- CH_all[i, -1] %>% as.numeric() # grab encounter history
    Zinits[i, CH_curr == 0] <- 7
    Zinits[i, CH_curr == 1] <- 6
    Zinits[i, first[i]] <- NA
  }
}

sel = (first == dim(CH_all[,-1])[2])
sum(sel) # since this is zero don't have any birds banded in last season to remove for dHMM. 

dat_marginal <- list(
  Y = (CH_all[, -1] %>% as.matrix()) + 1, 
  init = init,
  obs.mat.init = obs.mat.init,
  y.count = rowSums(dat$y.count), 
  J = dat$J, 
  R = dat$R
) 

# const ####

const_marginal <- list(
  n.occasions = const$n.occasions, 
  offset = const$offset, 
  juv.poss = const$juv.poss,
  n.years.fec = const$n.years.fec,
  ageCat = age,
  n.inds = dim(CH_all)[1],
  first = first,
  maxAge = const$n.occasions
) 

save(dat_marginal, const_marginal, file = "IPM_AEB_dat_stateSpace_marginal_fixed.RData")

# inits #####

inits_marginal <- list(
  ann.fec = iann.fec,
  mean.p.fidelity = 0.8, 
  mean.p.propensity = imean.p.propensity, 
  mean.p.det = imean.p.det, 
  mu.p.juv = imu.p.juv, 
  agebeta = iagebeta,
  sigma.p = isigma.p, 
  mean.phi.juv = imean.phi.juv,
  mean.phi.ad = imean.phi.ad,
  sigma.phi = isigma.phi, 
  IM = IMinits,
  N.recruits = iN.recruits,
  N.ad.surv = iN.ad.surv,
  N.breed.ready = iN.breed.ready,
  Ntot.breed = iNtot.breed,
  N.atsea = iN.atsea,
  nestlings = inestlings,
  JUV = iJUV
) 

save(inits_marginal, file = "IPM_AEB_inits_stateSpace_marginal_fixed.RData")
