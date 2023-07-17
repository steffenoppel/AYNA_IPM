load("IPM_AEB_dat.Rdata")
const <- c(const, maxAge = const$n.occasions)
maxage <- const$maxAge
n.years.fec <- const$n.years.fec

#### INITIAL VALUES ####
iann.fec <- rep(0.5, n.years.fec)
#imean.p.fidelity <- c(0.9, 0.9)
imean.p.propensity <- 0.8
imean.p.atsea <- 0.35
imean.p.det <- 3
imu.p.juv <- -6
iagebeta <- 0.75
isigma.p <- 0.1

imean.phi.juv <- 0.8
iinflation.factor <- 0.01
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
  IMinits[1, age, 1] = rbinom(1, IMinits[1, age-1, 3], imean.phi.ad*(1+iinflation.factor))
  IMinits[1, age, 2] <- rbinom(1, IMinits[1, age, 1], plogis(imu.p.juv + iagebeta*age))
  # print(paste(age, IMinits[1,age,1] - IMinits[1,age,2], sep = " "))
  IMinits[1, age, 3] <- IMinits[1,age,1] - IMinits[1,age,2]
}

IMinits[1, , ]
IMinits[,, 1]

iN.recruits[1] <- sum(IMinits[1,1:maxage,2]) 
iNtot.breed[1] <- max(rnorm(1, 640*0.5,sd = 20), 1) %>% round()#change here
iN.atsea[1] <- max(rnorm(1, 20*0.5,sd = 20), 1) %>% round() #change here
iN.loaf[1] <- max(rnorm(1, 120*0.5,sd = 20), 1) %>% round() #change here
iJUV[1] <- max(rnorm(1, 232*0.5, sd = 20), 1) %>% round()

for (tt in 2:n.years.fec) {
  iN.ad.surv[tt] <- rbinom(1, iNtot.breed[tt-1]+iN.atsea[tt-1]+iN.loaf[tt-1], imean.phi.ad)
  iN.breed.ready[tt] <- rbinom(1, iN.ad.surv[tt], 1-imean.p.atsea)
  iN.atsea[tt] <- iN.ad.surv[tt] - iN.breed.ready[tt]
  iN.loaf[tt] <- rbinom(1, iN.breed.ready[tt], 1-imean.p.propensity)
  
  IMinits[tt,1,1] <- rbinom(1, iJUV[tt-1], imean.phi.juv*(1+iinflation.factor))
  IMinits[tt, 1, 2] <- 0
  IMinits[tt, 1, 3] <- IMinits[tt,1,1] - IMinits[tt,1,2]
  
  for (age in 2:maxage) {
    IMinits[tt, age, 1] <- rbinom(1, IMinits[tt-1, age-1, 3], imean.phi.ad*(1+iinflation.factor))
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

AYNA_contacts_6states <- read_csv("AYNA_contacts_6states.csv") 
table(AYNA_contacts_6states$BreedState, AYNA_contacts_6states$AGE) 

AYNA_contacts_6states <- AYNA_contacts_6states %>% 
  mutate(BreedState = if_else(BreedState != 1 & AGE == 0, 1, BreedState), # cannot be 5 or 6 if chick...
         BreedState = if_else(BreedState == 1 & ContAge > 0, 3, BreedState) # cannot be chick if old...
  )
table(AYNA_contacts_6states$BreedState, AYNA_contacts_6states$AGE) 

#7 = not seen
#6 = known breeder
#5 = non-breeder, who has been recorded as known breeder before
# OMIT #4 = non-breeder, who has NOT been recorded as known breeder before. Omitting because generally all previous records are unknown
#3 = seen alive, breeding status unknown
# OMIT #2 = dead recovery (only 5 cases – I can’t remember whether you included them in the model or whether we ignored that for simplicity)
#1 = chick

#The first year from which we have complete information is 2019, so not many years unfortunately (2019,2020,2021).
completeDat <- c("2019-20", "2020-21", "2021-22")

all.seasons <- paste(1985:2021, "-", 
                     c(86:99, 
                       "00", "01", "02", "03", "04", "05", "06", "07", "08", "09", 
                       10:22), 
                     sep = "")
all.seasons

AYNA_CHICK<- AYNA_contacts_6states %>% 
  filter(FIRST_AGE %in% c("Chick","Fledgling")) %>%
  filter(Contact_Type != "Recovery") %>% 
  group_by(BirdID,Contact_Season) %>%
  mutate(Contact_Season = factor(Contact_Season, levels = all.seasons, labels = all.seasons)) %>% 
  summarise(STATE=max(BreedState)) %>% 
  mutate(STATE = if_else(Contact_Season %in% completeDat, STATE, 3)) %>% # everything unknown prior to good data
  mutate(STATE = case_when(
    STATE == 1 ~ 3, # only at the first occasion, recode
    #STATE == 2 ~ 2, # does not exist in data
    STATE == 3 ~ 3,
    STATE == 4 ~ 3, #combining these two 
    STATE == 5 ~ 5,
    STATE == 6 ~ 6,
    STATE == 7 ~ 7
  )) %>% 
  ungroup() %>% 
  complete(BirdID, Contact_Season, fill = list(STATE = 7)) %>% 
  pivot_wider(names_from=Contact_Season, values_from=STATE, values_fill=7) %>% 
  arrange(BirdID) %>% 
  select("BirdID", sort(colnames(.)))

AYNA_AD<- AYNA_contacts_6states %>% 
  filter(FIRST_AGE %in% c("Adult")) %>%    ### filter after spread to ensure that years without any adult contacts (1984, 2003, 2005) are included in matrix
  filter(Contact_Type != "Recovery") %>% 
  group_by(BirdID,Contact_Season) %>%
  mutate(Contact_Season = factor(Contact_Season, levels = all.seasons, labels = all.seasons)) %>% 
  summarise(STATE=max(BreedState)) %>% 
  mutate(STATE = if_else(Contact_Season %in% completeDat, STATE, 3)) %>% # everything unknown prior to good data
  mutate(STATE = case_when(
    STATE == 1 ~ 3, # only at the first occasion, recode
    #STATE == 2 ~ 2, # does not exist in data
    STATE == 3 ~ 3,
    STATE == 4 ~ 3, #combining these two 
    STATE == 5 ~ 5,
    STATE == 6 ~ 6,
    STATE == 7 ~ 7
  )) %>% 
  ungroup() %>% 
  complete(BirdID, Contact_Season, fill = list(STATE = 7)) %>% 
  pivot_wider(names_from=Contact_Season, values_from=STATE, values_fill=7) %>% 
  arrange(BirdID) %>% 
  select("BirdID", sort(colnames(.)))

AYNA_CHICK_reduced <- AYNA_CHICK %>% 
  group_by(across(-1)) %>% 
  summarise(mult = n())

AYNA_AD_reduced <- AYNA_AD %>% 
  group_by(across(-1)) %>% 
  summarise(mult = n())

diffTimeObs <- cbind(apply(AYNA_CHICK_reduced[, -38], 1, function(x){min(which(x!=7))}), 
                     apply(AYNA_CHICK_reduced[, -38], 1, function(x){which(x!=7)[2]})) %>% 
  as.data.frame() %>% 
  mutate(age_at_first_obs = V2 - V1, 
  )
rowInds <- which(diffTimeObs$age_at_first_obs <= 15)

AYNA_CHICK_reduced <- AYNA_CHICK_reduced[rowInds, ]

CH_all_reduced <- bind_rows(AYNA_AD_reduced, AYNA_CHICK_reduced)
mult <- CH_all_reduced$mult
CH_all_reduced <- CH_all_reduced %>% select(-mult)
CH_all_reduced[CH_all_reduced == 3] <- 3 #3 unknown status -> 3
CH_all_reduced[CH_all_reduced == 5] <- 2 #5 known loafer -> 2
CH_all_reduced[CH_all_reduced == 6] <- 1 #6 known breeder -> 1
CH_all_reduced[CH_all_reduced == 7] <- 4 #7 undetected -> 4
first <- apply(CH_all_reduced, 1, function(x){min(which(x != 4))})

age <- matrix(2, nrow = dim(CH_all_reduced)[1], ncol = dim(CH_all_reduced)[2])
for (i in (dim(AYNA_AD_reduced)[1] + 1):(dim(CH_all_reduced)[1])) {
  age[i, first[i]] <- 1
}

obs.mat.init.JUV <- matrix(
  c(
    0,0,1,0, # when banded as juvenile, assume state = prebreeder, observed = state unknown
    0,0,1,0,
    0,0,1,0,
    0,0,1,0,
    0,0,1,0,
    0,0,1,0
  ), nrow = 6, ncol = 4, byrow = T
)

#Zdat <- matrix(NA, nrow = dim(CH_all_reduced)[1], ncol = dim(CH_all_reduced)[2])
init <- matrix(NA, nrow = dim(CH_all_reduced)[1], ncol = 6)
obs.mat.init <- array(NA, dim = c(dim(CH_all_reduced)[1], 6, 4))
for (i in 1:nrow(init)) {
  if (age[i, first[i]] == 1) { # if age is 1 - prebreeder
    #Zdat[i, first[i]] <- 1 # TODO
    init[i, 1:6] <- c(1,0,0,0,0,0) 
    obs.mat.init[i, , ] <- obs.mat.init.JUV
    CH_curr <- CH_all_reduced[i,] %>% as.numeric() 
    ever_observed <- any(CH_curr[(first[i]+1):(dim(CH_all_reduced)[2])] != 4)
    if (ever_observed) {
      ever_observed_loafing <- any(CH_curr[(first[i]+1):(dim(CH_all_reduced)[2])] == 2)
      ever_observed_breeding<- any(CH_curr[(first[i]+1):(dim(CH_all_reduced)[2])] == 1)
      if (ever_observed_loafing & ever_observed_breeding) {
    
        when_observed_loafing <- (which(CH_curr[(1):(dim(CH_all_reduced)[2])] == 2))
        when_observed_breeding <- (which(CH_curr[(1):(dim(CH_all_reduced)[2])] == 1))
        first_observed_breeding <- min(when_observed_breeding)
        if (min(when_observed_loafing) < first_observed_breeding) {
          to_fix <- when_observed_loafing[which(when_observed_loafing < first_observed_breeding)]
          CH_all_reduced[i, to_fix] <- 1 # 3 TODO fix here
        }
      }
      only_observed_loafing <- all((CH_curr[(first[i]+1):(dim(CH_all_reduced)[2])])[CH_curr[(first[i]+1):(dim(CH_all_reduced)[2])] != 4] == 2)
      if (only_observed_loafing) {
        to_fix <- (CH_all_reduced[i, ] == 2) %>% as.logical()
        CH_all_reduced[i, to_fix] <- 1 # 3 TODO fix here
      }
      only_observed_unk_after2019 <- all((CH_curr[(dim(CH_all_reduced)[2]-2):(dim(CH_all_reduced)[2])])[CH_curr[(dim(CH_all_reduced)[2]-2):(dim(CH_all_reduced)[2])] != 4] == 3)
      if (only_observed_unk) {
        to_fix <- (CH_all_reduced[i, ] == 3) %>% as.logical()
        CH_all_reduced[i, to_fix] <- 1 # 3 TODO fix here
      }
    }
  } else { # if age is 2 - breeding in colony
    #Zdat[i, first[i]] <- 8 
    if (CH_all_reduced[i, first[i]] == 1) { # seen breeding
      init[i, 1:6] <- c(0,0,0,0,1,0) 
      obs.mat.init.AD <- matrix(
        c(
          0,0,1,0,
          0,0,1,0,
          0,0,1,0,
          0,0,1,0,
          1,0,0,0, # when banded as an adult, assume state = prebreeder, observed = state unknown
          0,0,1,0
        ), nrow = 6, ncol = 4, byrow = T
      )
    } else if (CH_all_reduced[i, first[i]] == 2) { # seen loafing
      init[i, 1:6] <- c(0,0,1,0,0,0) 
      obs.mat.init.AD <- matrix(
        c(
          0,0,1,0,
          0,0,1,0,
          0,1,0,0,
          0,0,1,0,
          0,0,1,0, # when banded as an adult, assume state = prebreeder, observed = state unknown
          0,0,1,0
        ), nrow = 6, ncol = 4, byrow = T
      )
    } else if (CH_all_reduced[i, first[i]] == 3) { # seen unknown
      init[i, 1:6] <- c(0,0,0,0,1,0) 
      obs.mat.init.AD <- matrix(
        c(
          0,0,1,0,
          0,0,1,0,
          0,0,1,0,
          0,0,1,0,
          0,0,1,0, # when banded as an adult, assume state = prebreeder, observed = state unknown
          0,0,1,0
        ), nrow = 6, ncol = 4, byrow = T
      )
    } 
    obs.mat.init[i, , ] <- obs.mat.init.AD
  }
}

# Zinits <- matrix(NA, nrow = dim(CH_all_reduced)[1], ncol = dim(CH_all_reduced)[2])
# for (i in 1:nrow(Zinits)) {
#   if (age[i, first[i]] == 1) {
#     CH_curr <- CH_all_reduced[i,] %>% as.numeric() 
#     ever_observed <- any(CH_curr[(first[i]+1):(dim(CH_all_reduced)[2])] == 1)
#     if (ever_observed) {
#       all_inds <- (which(CH_curr[(1):(dim(CH_all_reduced)[2])] == 1))
#       inds_2plus <- all_inds[2:length(all_inds)]
#       Zinits[i, first[i]:(inds_2plus[1]-1)] <- 1
#       Zinits[i, inds_2plus[1]-1] <- 4
#       Zinits[i, min(inds_2plus[1]+1, dim(CH_all_reduced)[2]):(dim(CH_all_reduced)[2])] <- 5
#       Zinits[i, inds_2plus] <- 5
#     } else {
#       for (t in (first[i]+1):(dim(CH_all_reduced)[2])) {
#         curr_age <- t - first[i]
#         if (curr_age < 15) {
#           Zinits[i, t] <- 1
#         } else if (curr_age == 10) {
#           Zinits[i, t] <- 4
#         } else if (curr_age > 10) {
#           Zinits[i, t] <- 5
#         }
#       } 
#     }
#     Zinits[i, first[i]] <- NA
#   } else { # if age is 2 - breeding in colony
#     CH_curr <- CH_all_reduced[i,] %>% as.numeric()  # grab encounter history
#     Zinits[i, CH_curr == 0] <- 5
#     Zinits[i, CH_curr == 1] <- 5
#     Zinits[i, first[i]] <- NA
#   }
# }

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
  Y = (CH_all_reduced %>% as.matrix()), 
  init = init,
  obs.mat.init = obs.mat.init,
  y.count = rowSums(dat$y.count), 
  J = dat$J, 
  R = dat$R,
  ICCAT.ll.mit = longline$mit.ICCAT,
  Nam.ll.mit = longline$mit.NAM,
  SA.ll.mit = longline$mit.RSA,
  Uru.ll.mit = longline$mit.URU,
  phi_constraint_data = 1
) 

# const ####

const_marginal <- list(
  n.occasions = const$n.occasions, 
  offset = const$offset,
  n.years.fec = const$n.years.fec,
  ageCat = age,
  n.inds = dim(CH_all_reduced)[1],
  first = first,
  maxAge = const$n.occasions,
  mult = mult,
  gamma = c(rep(0, times = (length(all.seasons)-3)), rep(1, times = 3))
) 

save(dat_marginal, const_marginal, file = "IPM_AEB_dat_stateSpace_marginal_loaf_reduced_incolony_informed_COVARIATES.RData")

# inits #####

inits_marginal <- list(
  mean.fec = 0.5, 
  sigma.fec = 0.1,
  mean.p.propensity = imean.p.propensity, 
  mean.p.atsea = imean.p.atsea,
  mean.p.det = imean.p.det, 
  mu.p.juv = imu.p.juv, 
  agebeta = iagebeta,
  sigma.p = isigma.p, 
  mean.phi.juv = imean.phi.juv,
  mean.phi.ad = imean.phi.ad,
  inflation.factor = iinflation.factor,
  sigma.phi = isigma.phi, 
  IM = IMinits,
  N.ad.surv = iN.ad.surv,
  N.breed.ready = iN.breed.ready,
  Ntot.breed = iNtot.breed,
  N.atsea = iN.atsea,
  N.loaf = iN.loaf,
  JUV = iJUV,
  sigma.beta1 = 0.1,
  beta.mit = c(0,0),
  alpha = c(10, 10, 10, 10),
  w = c(0.25, 0.25, 0.25, 0.25)
) 

save(inits_marginal, file = "IPM_AEB_inits_stateSpace_marginal_loaf_reduced_incolony_informed_COVARIATES.RData")

