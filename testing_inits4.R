load("IPM_AEB_dat.Rdata")
const <- c(const, maxAge = const$n.occasions)
maxage <- const$maxAge
n.years.fec <- const$n.years.fec

#### INITIAL VALUES ####
iann.fec <- rep(0.5, n.years.fec)
#imean.p.fidelity <- c(0.9, 0.9)
imean.p.fidelity <- 0.9
imean.p.propensity <- 0.7
imean.p.atsea <- 0.05
imean.p.det <- 2.5
imu.p.juv <- -3
iagebeta <- 0.75
isigma.p <- 0.1

imean.phi.juv <- 0.4
imean.phi.ad <- 0.9
isigma.phi <- 0.1

iN.ad.surv <- rep(0, n.years.fec)
iN.breed.ready <- rep(0, n.years.fec)
inestlings <- rep(0, n.years.fec)

iN.recruits <- rep(NA, n.years.fec)
iNtot.breed <- rep(NA, n.years.fec)
iN.atsea <- rep(NA, n.years.fec)
iN.loaf <- rep(NA, n.years.fec)
iJUV <- rep(NA, n.years.fec)
IMinits <- array(NA, dim = c(n.years.fec, maxage, 3))

IMinits[1,1,1] = max(c(rnorm(1, 263*0.5, 20)), 1) %>% round() 
IMinits[1,1,2] <- 0
IMinits[1,1,3] <- round(IMinits[1,1,1]) - IMinits[1,1,2]

IMinits[1,2,1] = max(c(rnorm(1, 275*0.5, 20)), 1)%>% round()
IMinits[1,2,2] <- rbinom(1, IMinits[1, 2, 1], plogis(imu.p.juv + iagebeta*2))
IMinits[1,2,3] <- round(IMinits[1,2,1]) - IMinits[1,2,2]

IMinits[1,3,1] = max(c(rnorm(1, 264*0.5, 20)), 1)%>% round()
IMinits[1,3,2] <- rbinom(1, IMinits[1, 3, 1], plogis(imu.p.juv + iagebeta*3))
IMinits[1,3,3] <- round(IMinits[1,3,1]) - IMinits[1,3,2]

IMinits[1,4,1] = max(c(rnorm(1, 177*0.5, 20)), 1)%>% round()
IMinits[1,4,2] <- rbinom(1, IMinits[1, 4, 1], plogis(imu.p.juv + iagebeta*4))
IMinits[1,4,3] <- round(IMinits[1,4,1]) - IMinits[1,4,2]

IMinits[1,5,1] = max(c(rnorm(1, 290*0.5, 20)), 1)%>% round()
IMinits[1,5,2] <- rbinom(1, IMinits[1, 5, 1], plogis(imu.p.juv + iagebeta*5))
IMinits[1,5,3] <- round(IMinits[1,5,1]) - IMinits[1,5,2]

IMinits[1,6,1] = max(c(rnorm(1, 90*0.5, 20)), 1)%>% round()
IMinits[1,6,2] <- rbinom(1, IMinits[1, 6, 1], plogis(imu.p.juv + iagebeta*6))
IMinits[1,6,3] <- round(IMinits[1,6,1]) - IMinits[1,6,2]

IMinits[1,7,1] = max(c(rnorm(1, 158*0.5, 20)), 1)%>% round()
IMinits[1,7,2] <- rbinom(1, IMinits[1, 7, 1], plogis(imu.p.juv + iagebeta*7))
IMinits[1,7,3] <- round(IMinits[1,7,1]) - IMinits[1,7,2]

for (age in 8:maxage) {
  #IMinits[1, age, 1] = rbinom(1, IMinits[1, age-1, 3], pow(0.9,(age-1)))
  # AEB change
  IMinits[1, age, 1] = rbinom(1, IMinits[1, age-1, 3], imean.phi.ad)
  IMinits[1, age, 2] <- rbinom(1, IMinits[1, age, 1], plogis(imu.p.juv + iagebeta*age))
  # print(paste(age, IMinits[1,age,1] - IMinits[1,age,2], sep = " "))
  IMinits[1, age, 3] <- IMinits[1,age,1] - IMinits[1,age,2]
}

IMinits[1, , ]
IMinits[,, 1]

iN.recruits[1] <- sum(IMinits[1,1:maxage,2]) 
iNtot.breed[1] <- max(rnorm(1, 640*0.5,sd = 20), 1) %>% round()#change here
iN.atsea[1] <- max(rnorm(1, 32*0.5,sd = 20), 1) %>% round() #change here
iN.loaf[1] <- max(rnorm(1, 120*0.5,sd = 20), 1) %>% round() #change here
iJUV[1] <- max(rnorm(1, 232*0.5, sd = 20), 1) %>% round()

for (tt in 2:n.years.fec) {
  iN.ad.surv[tt] <- rbinom(1, iNtot.breed[tt-1]+iN.atsea[tt-1]+iN.loaf[tt-1], imean.phi.ad)
  iN.breed.ready[tt] <- rbinom(1, iN.ad.surv[tt], 1-imean.p.atsea)
  iN.atsea[tt] <- iN.ad.surv[tt] - iN.breed.ready[tt]
  iN.loaf[tt] <- rbinom(1, iN.breed.ready[tt], 1-imean.p.propensity)
  
  IMinits[tt,1,1] <- rbinom(1, iJUV[tt-1], imean.phi.juv)
  IMinits[tt, 1, 2] <- 0
  IMinits[tt, 1, 3] <- IMinits[tt,1,1] - IMinits[tt,1,2]
  
  for (age in 2:maxage) {
    IMinits[tt, age, 1] <- rbinom(1, IMinits[tt-1, age-1, 3], imean.phi.ad)
    #IMinits[tt, age, 2] <- min(IMinits[tt, age, 1], IMinits[tt, age-1, 3]) * p.juv.recruit.inits[age, tt] %>% round()
    IMinits[tt, age, 2] <- rbinom(1, IMinits[tt, age, 1], plogis(imu.p.juv + iagebeta*age))
    IMinits[tt, age, 3] <- IMinits[tt,age,1] - IMinits[tt,age,2]
  }
  
  iN.recruits[tt] <- sum(IMinits[tt,1:maxage,2])
  iNtot.breed[tt]<- iN.breed.ready[tt]-iN.loaf[tt]+iN.recruits[tt]
  inestlings[tt] <- (iann.fec[tt] * 0.5 * iNtot.breed[tt]) %>% round()
  iJUV[tt] <- rpois(1, inestlings[tt])
  
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

# REDUCED #######

AYNA_AD_reduced <- AYNA_AD %>% 
  group_by(across(-1)) %>% 
  summarise(mult = n())

AYNA_CHICK_reduced <- AYNA_CHICK %>% 
  group_by(across(-1)) %>% 
  summarise(mult = n())

CH_all_reduced <- bind_rows(AYNA_AD_reduced, AYNA_CHICK_reduced)
mult <- CH_all_reduced$mult
CH_all_reduced <- CH_all_reduced %>% select(-mult)
first <- apply(CH_all_reduced, 1, function(x){min(which(x != 0))})

age <- matrix(2, nrow = dim(CH_all_reduced)[1], ncol = dim(CH_all_reduced)[2])
for (i in (dim(AYNA_AD_reduced)[1] + 1):(dim(CH_all_reduced)[1])) {
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
    1,0,
    1,0,
    1,0
  ), nrow = 10, ncol = 2, byrow = T
)
obs.mat.init.AD <- matrix(
  c(
    1,0,
    1,0,
    1,0,
    1,0,
    1,0,
    1,0,
    1,0,
    0,1, # when banded as an adult, assume observed and state = breeding in colony
    1,0,
    1,0
  ), nrow = 10, ncol = 2, byrow = T
)

Zdat <- matrix(NA, nrow = dim(CH_all_reduced)[1], ncol = dim(CH_all_reduced)[2])
init <- matrix(NA, nrow = dim(CH_all_reduced)[1], ncol = 10)
obs.mat.init <- array(NA, dim = c(dim(CH_all_reduced)[1], 10, 2))
for (i in 1:nrow(Zdat)) {
  if (age[i, first[i]] == 1) { # if age is 1 -alive at sea never bred
    Zdat[i, first[i]] <- 1 
    init[i, 1:10] <- c(1,0,0,0,0,0,0,0,0,0)
    obs.mat.init[i, , ] <- obs.mat.init.JUV
  } else { # if age is 2 - breeding in colony
    Zdat[i, first[i]] <- 8 
    init[i, 1:10] <- c(0,0,0,0,0,0,0,1,0,0)
    obs.mat.init[i, , ] <- obs.mat.init.AD
  }
}

Zinits <- matrix(NA, nrow = dim(CH_all_reduced)[1], ncol = dim(CH_all_reduced)[2])
for (i in 1:nrow(Zinits)) {
  if (age[i, first[i]] == 1) {
    CH_curr <- CH_all_reduced[i,] %>% as.numeric() 
    ever_observed <- any(CH_curr[(first[i]+1):(dim(CH_all_reduced)[2])] == 1)
    if (ever_observed) {
      all_inds <- (which(CH_curr[(1):(dim(CH_all_reduced)[2])] == 1))
      inds_2plus <- all_inds[2:length(all_inds)]
      Zinits[i, first[i]:(inds_2plus[1]-1)] <- 1
      Zinits[i, inds_2plus[1]-1] <- 7
      Zinits[i, min(inds_2plus[1]+1, dim(CH_all_reduced)[2]):(dim(CH_all_reduced)[2])] <- 9
      Zinits[i, inds_2plus] <- 8
    } else {
      for (t in (first[i]+1):(dim(CH_all_reduced)[2])) {
        curr_age <- t - first[i]
        if (curr_age < 10) {
          Zinits[i, t] <- 1
        } else if (curr_age == 10) {
          Zinits[i, t] <- 7
        } else if (curr_age > 10) {
          Zinits[i, t] <- 9
        }
      } 
    }
    Zinits[i, first[i]] <- NA
  } else { # if age is 2 - breeding in colony
    CH_curr <- CH_all_reduced[i,] %>% as.numeric()  # grab encounter history
    Zinits[i, CH_curr == 0] <- 9
    Zinits[i, CH_curr == 1] <- 8
    Zinits[i, first[i]] <- NA
  }
}

sel = (first == dim(CH_all_reduced)[2])
sum(sel) # since this is zero don't have any birds banded in last season to remove for dHMM. 

longline <- longline %>% mutate(n_hooks = scale(n_hooks)) 
ave.since.2010 <- longline %>% filter(Year > 2010) %>% select(2) 
ave.since.2010 <- mean(ave.since.2010$n_hooks, na.rm = T)
longline <- longline %>% 
  mutate(n_hooks = if_else(Year %in% c(2020, 2021), ave.since.2010, n_hooks)) %>% 
  filter(Year >= 1985)
longline

dat_marginal <- list(
  Y = (CH_all_reduced %>% as.matrix()) + 1, 
  init = init,
  obs.mat.init = obs.mat.init,
  y.count = rowSums(dat$y.count), 
  J = dat$J, 
  R = dat$R,
  ICCAT.ll.e = longline$n_hooks,
  ICCAT.ll.mit = longline$mit.ICCAT,
  Nam.ll.mit = longline$mit.NAM,
  SA.ll.mit = longline$mit.RSA,
  Uru.ll.mit = longline$mit.URU
) 

# const ####

const_marginal <- list(
  n.occasions = const$n.occasions, 
  offset = const$offset, 
  juv.poss = const$juv.poss,
  n.years.fec = const$n.years.fec,
  ageCat = age,
  n.inds = dim(CH_all_reduced)[1],
  first = first,
  maxAge = const$n.occasions,
  mult = mult
) 

save(dat_marginal, const_marginal, file = "IPM_AEB_dat_stateSpace_marginal_loaf_reduced_COVARIATES.RData")

# inits #####

inits_marginal <- list(
  mean.fec = 0.5, 
  sigma.fec = 0.1,
  mean.p.fidelity = imean.p.fidelity, 
  mean.p.propensity = imean.p.propensity, 
  mean.p.atsea = imean.p.atsea,
  mean.p.det = imean.p.det, 
  mu.p.juv = imu.p.juv, 
  agebeta = iagebeta,
  sigma.p = isigma.p, 
  mean.phi.juv = imean.phi.juv,
  mean.phi.ad = imean.phi.ad,
  sigma.phi = isigma.phi, 
  IM = IMinits,
  #N.recruits = iN.recruits,
  N.ad.surv = iN.ad.surv,
  N.breed.ready = iN.breed.ready,
  Ntot.breed = iNtot.breed,
  N.atsea = iN.atsea,
  N.loaf = iN.loaf,
  #nestlings = inestlings,
  JUV = iJUV,
  sigma.beta1 = 0.1,
  sigma.beta2 = 0.1,
  sigma.beta3 = 0.1,
  sigma.beta4 = 0.1,
  sigma.beta5 = 0.1,
  beta.ICCAT.ll.e = c(0,0),
  beta.ICCAT.ll.mit = c(0,0),
  beta.Nam.ll.mit = c(0,0),
  beta.SA.ll.mit = c(0,0),
  beta.Uru.ll.mit = c(0,0)
) 

save(inits_marginal, file = "IPM_AEB_inits_stateSpace_marginal_loaf_reduced_COVARIATES.RData")

