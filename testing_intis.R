maxage <- const$maxAge
n.years.fec <- const$n.years.fec

#### INITIAL VALUES ####
iann.fec <- rep(0.4, n.years.fec)

isigma.obs <- rexp(1, 0.1)

imean.p.ad <- c(runif(1, 0.05, 0.5), runif(1, 0.2, 1))
imu.p.juv <- -4
iagebeta <- 2
isigma.p <- 0.1

imean.phi.juv <- 0.75
imean.phi.ad <- 0.95
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
IMinits[1,2,2] <- 0
IMinits[1,2,3] <- round(IMinits[1,2,1]) - IMinits[1,2,2]

IMinits[1,3,1] = max(c(rnorm(1, 264*0.5, 20)), 1)%>% round()
IMinits[1,3,2] <- 0
IMinits[1,3,3] <- round(IMinits[1,3,1]) - IMinits[1,3,2]

IMinits[1,4,1] = max(c(rnorm(1, 177*0.5, 20)), 1)%>% round()
IMinits[1,4,2] <- 0
IMinits[1,4,3] <- round(IMinits[1,4,1]) - IMinits[1,4,2]

IMinits[1,5,1] = max(c(rnorm(1, 290*0.5, 20)), 1)%>% round()
IMinits[1,5,2] <- 0
IMinits[1,5,3] <- round(IMinits[1,5,1]) - IMinits[1,5,2]

IMinits[1,6,1] = max(c(rnorm(1, 90*0.5, 20)), 1)%>% round()
IMinits[1,6,2] <- 0
IMinits[1,6,3] <- round(IMinits[1,6,1]) - IMinits[1,6,2]

IMinits[1,7,1] = max(c(rnorm(1, 158*0.5, 20)), 1)%>% round()
IMinits[1,7,2] <- 0
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
  iN.breed.ready[tt] <- rbinom(1, iN.ad.surv[tt], imean.p.ad)
  
  IMinits[tt,1,1] <- rbinom(1, iJUV[tt-1], imean.phi.juv)
  IMinits[tt, 1, 2] <- 0
  IMinits[tt, 1, 3] <- IMinits[tt,1,1] - IMinits[tt,1,2]
  
  for (age in 2:maxage) {
    IMinits[tt, age, 1] <- rbinom(1, IMinits[tt-1, age-1, 3], imean.phi.ad)
    # TODO wtf is this
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
(iNtot.breed * const$prop.sites[1:13, ] * 2 ) - dat$y.count[1:13, ]
mean((iNtot.breed * const$prop.site[1:13, ]  * 2 ) - dat$y.count[1:13, ])

# logic check
# there should not be negatives here
#   for (year in 2:n.years.fec) {
#     for(age in 2:maxage) {
#     print(paste(year, age, IMinits[year-1, age-1, 3] - IMinits[year, age, 1] , sep = ' '))
#   }
# }

inits <- list(
  
  # fecundity stuff
  ann.fec = iann.fec,
  
  # count stuff
  sigma.obs=matrix(rexp(const$n.sites.count*2, isigma.obs),ncol=2),
  
  # recapture stuff
  mean.p.ad = imean.p.ad, 
  mu.p.juv = rep(imu.p.juv, 2),
  agebeta = iagebeta,
  sigma.p = isigma.p, # fixed because of the IM 
  
  # survival stuff, TODO fix
  mean.phi.ad = imean.phi.ad , # fixed because of the IM
  mean.phi.juv = imean.phi.juv, # fixed because of the IM
  sigma.phi = isigma.phi, # fixed because of the IM
  
  # abundance stuff
  IM = IMinits, 
  N.ad.surv = iN.ad.surv,
  N.breed.ready = iN.breed.ready,
  nestlings = inestlings,
  N.recruits = iN.recruits,
  Ntot.breed = iNtot.breed,
  JUV = iJUV,
  N.atsea = iN.atsea
) 

