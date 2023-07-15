# LIBRARIES #####
library(nimble)
library(coda)
library(tidyverse)
library(tidybayes)
library(strex)
library(beepr)
library(postpack)
library(here)
library(ggdist)
library(beepr)
library(wesanderson)

# LOAD SAMPLES #####
load("~/Documents/AYNA_IPM/samples_statespace_marginal_loaf_reduced_incolony_new_COVARIATES_chain1.Rdata")
chain1 <- out1
load("~/Documents/AYNA_IPM/samples_statespace_marginal_loaf_reduced_incolony_new_COVARIATES_chain2.Rdata")
chain2 <- out1
load("~/Documents/AYNA_IPM/samples_statespace_marginal_loaf_reduced_incolony_new_COVARIATES_chain3.Rdata")
chain3 <- out1

out1 <- list(chain1 = chain1, chain2 = chain2, chain3 = chain3) %>% 
  as.mcmc.list()
nb <- dim(chain1)[1]/2
out1_wburnin <- lapply(out1, function(x) x[(nrow(x)-nb+1):nrow(x), ] %>% as.mcmc()) %>% 
  as.mcmc.list()
out1_wburnin_thinned <- post_thin(out1_wburnin, keep_iters = nb/10)

summ <- t(post_summ(out1_wburnin_thinned, get_params(out1_wburnin_thinned, type = "base_index"), 
                    neff = TRUE, Rhat = TRUE, probs = c(0.025, 0.5, 0.975))) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name")
beep(sound = 1)

# SET UP SCENARIOS ####

IM_samps <- post_subset(out1_wburnin_thinned, "IM", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains("[13,")) %>% 
  pivot_longer(everything()) %>% 
  mutate(age = str_nth_number(name, 2), 
         col = str_nth_number(name, 3)) %>% 
  group_by(age, col) %>% 
  mutate(index = row_number())

IM_samps <- array(data = IM_samps$value, 
      dim=c(length(unique(IM_samps$age)), 
            length(unique(IM_samps$col)), 
            length(unique(IM_samps$index))), 
      dimnames=list(unique(IM_samps$age), unique(IM_samps$col), unique(IM_samps$index))
)

N.recruits_samps <- post_subset(out1_wburnin_thinned, "N.recruits", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains("[13]"))
N.ad.surv_samps <- post_subset(out1_wburnin_thinned, "N.ad.surv", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains("[13]"))
N.breed.ready_samps <- post_subset(out1_wburnin_thinned, "N.breed.ready", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains("[13]"))
N.atsea_samps <- post_subset(out1_wburnin_thinned, "N.atsea", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains("[13]"))
N.loaf_samps <- post_subset(out1_wburnin_thinned, "N.loaf", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains("[13]"))
Ntot.breed_samps <- post_subset(out1_wburnin_thinned, "Ntot.breed", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains("[13]"))
nestlings_samps <- post_subset(out1_wburnin_thinned, "nestlings", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains("[13]"))
JUV_samps <- post_subset(out1_wburnin_thinned, "JUV", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains("[13]"))
Ntot_samps <- post_subset(out1_wburnin_thinned, "Ntot", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  select(contains("[13]"))

ages <- matrix(1:37, nrow = 15000, ncol = length(1:37), byrow = T)
p.juv.recruit <- post_subset(out1_wburnin_thinned, "mu.p.juv|agebeta", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  bind_cols(ages) %>% 
  mutate(across(3:39, ~ plogis(mu.p.juv + agebeta * .)))

# NOTE - have to take logit of mean.fec and then random deviate from sigma fec
ann.fec <- post_subset(out1_wburnin_thinned, "mean.fec|eps.fec|sigma.fec", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  rowwise() %>% 
  mutate(across(1:13, ~ plogis(log(mean.fec / (1-mean.fec)) + .))) %>% 
  mutate(rand.fec = rnorm(1, 0, sigma.fec)) %>% 
  mutate(rand.fec = plogis(log(mean.fec / (1-mean.fec)) + rand.fec))

# TODO - these will change

w <- post_subset(out1_wburnin_thinned, "w\\[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() 

# w <- matrix(NA, nrow = 15000, ncol = 4)
# for (i in 1:15000) {
#   w[i, 1:4] <- rdirch(1, c(0.4, 0.1, 0.1, 0.4))
#   #w[i, 1:4] <- c(0, 0.05, 0.05, 0.9)
# }

# TODO - beta mit
# for now let's assume beta.mit is positive, and less than 0.25
# since (log(0.9 / (1-0.9)) + 0.0*4) %>% plogis() vs log(0.9 / (1-0.9)) + 0.25*4) %>% plogis() is a 6% increase in survival...
beta.mit <- post_subset(out1_wburnin_thinned, "beta.mit", matrix = T, iters = F, chains = F) %>% 
  as.data.frame()

# beta.mit1 <- rbeta(15000, 10, 80)
# beta.mit2 <- rbeta(15000, 15, 80)

# NOTE - have to take logit of mean.phi and then random deviate from sigma phi
# mit.index[j] <- w1*ICCAT.ll.mit[j] + w2*Nam.ll.mit[j] + w3*SA.ll.mit[j] + w4*Uru.ll.mit[j]
# logit(phi.ad[j])  <- mu.ad  + eps.phi[j]  + beta.mit[2]*mit.index[j]
phi.ad <- post_subset(out1_wburnin_thinned, "mean.phi.ad|sigma.phi|eps.phi", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  bind_cols(w, beta.mit = beta.mit$`beta.mit[2]`)

phi.im <- post_subset(out1_wburnin_thinned, "mean.phi.im|sigma.phi|eps.phi", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  bind_cols(w, beta.mit = beta.mit$`beta.mit[1]`)

# NOTE - have to take logit of mean.phi and then random deviate from sigma phi
# mit.index[j] <- w1*ICCAT.ll.mit[j] + w2*Nam.ll.mit[j] + w3*SA.ll.mit[j] + w4*Uru.ll.mit[j]
# logit(phi.ad[j])  <- mu.ad  + eps.phi[j]  + beta.mit[2]*mit.index[j]
phi.juv <- post_subset(out1_wburnin_thinned, "mean.phi.juv|sigma.phi|eps.phi", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  bind_cols(w, beta.mit = beta.mit$`beta.mit[1]`)

mean.p.atsea_samps <- post_subset(out1_wburnin_thinned, "mean.p.atsea", matrix = T, iters = F, chains = F) %>% 
  as.data.frame()

mean.p.propensity_samps <- post_subset(out1_wburnin_thinned, "mean.p.propensity", matrix = T, iters = F, chains = F) %>% 
  as.data.frame()

# grab indices
inds <- sample(1:(dim(mean.p.propensity_samps)[1]), 10000, replace = F)

# scenarios
# 1. Status quo
# 3. Increase mitigation -> higher survival, juveniles stand to benefit more
# With some variation in effect size
# 5. Decrease mitigation/increase bycatch -> lower survival, juveniles stand to lose more
# With some variation in effect size

# IGNORE UNLESS ESTIMATES SHOW NO EFFECT
# # 4. Decrease mitigation/increase bycatch -> lower survival, uniform decrease across ages
# With some variation in effect size
# # 2. Increase mitigation -> higher survival, uniform increase across ages
# With some variation in effect size

# simulate future mitigation across fisheries

# in 2021 - NAM is at 1 - assume continues
#           RSA is at 0.9 - assume goes to 1
#           URU is at 1 - assume continues
#           ICCAT is at 0.1 - assume increases by 0.05 per year for 19 years (max 1)
n.future <- 50
mit.NAM <- rep(1, n.future)
mit.RSA <- rep(1, n.future)
mit.URU <- rep(1, n.future)
mit.ICCAT <- rep(1, n.future)
#mit.ICCAT <- c(seq(from = 0.1, to = 1, by = 0.05), rep(1, n.future-length(seq(from = 0.1, to = 1, by = 0.05))))
future.mit.INCREASE <- cbind(mit.NAM, mit.RSA, mit.URU, mit.ICCAT)

# in 2021 - NAM is at 1 
#           RSA is at 0.9 
#           URU is at 1 
#           ICCAT is at 0.1
future.mit.STATUSQUO <- matrix(c(1, 0.9, 1, 0.1), nrow = n.future, ncol = 4, byrow = T)
colnames(future.mit.STATUSQUO) <- c("mit.NAM", "mit.RSA", "mit.URU", "mit.ICCAT")

# in 2021 - NAM goes to 0 
#           RSA goes to 0
#           URU goes to 0
#           ICCAT goes to 0
future.mit.DECREASE <- matrix(c(0, 0, 0, 0), nrow = n.future, ncol = 4, byrow = T)
colnames(future.mit.DECREASE) <- c("mit.NAM", "mit.RSA", "mit.URU", "mit.ICCAT")

# RUN PVA ####
maxAge <- 37
IM <- array(NA, dim = c(3, n.future+1, maxAge, 3, 10000))
N.recruits <- array(NA, dim = c(3, n.future+1, 10000))
N.ad.surv <- array(NA, dim = c(3, n.future+1, 10000))
N.breed.ready <- array(NA, dim = c(3, n.future+1, 10000))
N.atsea <- array(NA, dim = c(3, n.future+1, 10000))
N.loaf <- array(NA, dim = c(3, n.future+1, 10000))
Ntot.breed <- array(NA, dim = c(3, n.future+1, 10000))
nestlings <- array(NA, dim = c(3, n.future+1, 10000))
JUV <- array(NA, dim = c(3, n.future+1, 10000))
Ntot <- array(NA, dim = c(3, n.future+1, 10000))

for (s in 1:3) { # TODO maybe adding more scenarios?
  if (s == 1) {
    future.mit <- future.mit.STATUSQUO
  } else if (s == 2) {
    future.mit <- future.mit.INCREASE
  } else if (s == 3) {
    future.mit <- future.mit.DECREASE
  }
  
  # set up first year
  IM[s, 1, , , ] <- IM_samps[, , inds]
  N.recruits[s, 1, ] <- (N.recruits_samps %>% unlist() %>% unname())[inds]
  N.ad.surv[s, 1, ] <- (N.ad.surv_samps %>% unlist() %>% unname())[inds]
  N.breed.ready[s, 1, ] <- (N.breed.ready_samps %>% unlist() %>% unname())[inds]
  N.atsea[s, 1, ] <- (N.atsea_samps %>% unlist() %>% unname())[inds]
  N.loaf[s, 1, ] <- (N.loaf_samps %>% unlist() %>% unname())[inds]
  Ntot.breed[s, 1, ] <- (Ntot.breed_samps %>% unlist() %>% unname())[inds]
  nestlings[s, 1, ] <- (nestlings_samps %>% unlist() %>% unname())[inds]
  JUV[s, 1, ] <- (JUV_samps %>% unlist() %>% unname())[inds]
  Ntot[s, 1, ] <- (Ntot_samps %>% unlist() %>% unname())[inds]
  
  # loop over future years
  for (t in 2:(n.future+1)) {
    print(paste("scenario", s, "year", t, sep = " "))
    
    tmp_mit.index <- w$`w[1]`[inds]*future.mit[t-1, 4] + w$`w[2]`[inds]*future.mit[t-1, 1] + w$`w[3]`[inds]*future.mit[t-1, 2] + w$`w[4]`[inds]*future.mit[t-1, 3]
    tmp_eps.phi <- rnorm(10000, 0, phi.ad$sigma.phi[inds])
    #tmp_eps.phi <- 0
    
    tmp_mu.ad <- log(phi.ad$mean.phi.ad[inds] / (1 - phi.ad$mean.phi.ad[inds]))
    tmp_beta.mit.ad <- phi.ad$beta.mit[inds]
    tmp_phi.ad <- ((tmp_mu.ad + tmp_eps.phi + tmp_beta.mit.ad*tmp_mit.index) %>% plogis())
    
    tmp_mu.juv <- log(phi.juv$mean.phi.juv[inds] / (1 - phi.juv$mean.phi.juv[inds]))
    tmp_beta.mit.juv <- phi.juv$beta.mit[inds]
    tmp_phi.juv <- ((tmp_mu.juv + tmp_eps.phi + tmp_beta.mit.juv*tmp_mit.index) %>% plogis()) 
    
    tmp_mu.im <- log(phi.im$mean.phi.im[inds] / (1 - phi.im$mean.phi.im[inds]))
    tmp_beta.mit.im <- phi.im$beta.mit[inds]
    tmp_phi.im <- ((tmp_mu.im + tmp_eps.phi + tmp_beta.mit.im*tmp_mit.index) %>% plogis())
    
    p.juv.recruit.f <- p.juv.recruit[inds, -c(1,2)]
    
    #mean.fec <- ann.fec$ann.fec[inds]
    #mean.fec <- ann.fec$mean.fec[inds]
    mean.fec <- ann.fec$rand.fec[inds] 
    
    mean.p.atsea <- mean.p.atsea_samps$mean.p.atsea[inds]
    mean.p.propensity <- mean.p.propensity_samps$mean.p.propensity[inds]
    
    IM[s,t,1,1, ] <- rbinom(10000, JUV[s, t-1, ], tmp_phi.juv)                                
    IM[s,t,1,2, ] <- 0 
    IM[s,t,1,3, ] <- IM[s,t,1,1, ] - IM[s,t,1,2, ]
    for(age in 2:maxAge) {
      IM[s,t,age,1,] <- rbinom(10000, IM[s,t-1,age-1,3, ], tmp_phi.im)  
      IM[s,t,age,2,] <- rbinom(10000, IM[s,t,age,1,], p.juv.recruit.f[, age])
      IM[s,t,age,3,] <- IM[s,t,age,1,] - IM[s,t,age,2,]
    }
    
    N.recruits[s,t,] <- colSums(IM[s,t,,2,])
    N.ad.surv[s,t,] <- rbinom(10000, Ntot.breed[s,t-1,] + N.atsea[s,t-1,] + N.loaf[s,t-1,], tmp_phi.ad)
    N.breed.ready[s,t,] <- rbinom(10000, N.ad.surv[s,t,], 1-mean.p.atsea)
    N.atsea[s,t,] <- N.ad.surv[s,t,] - N.breed.ready[s,t,]
    N.loaf[s,t,] <- rbinom(10000, N.breed.ready[s,t,], 1-mean.p.propensity)
    Ntot.breed[s,t,] <- N.breed.ready[s,t,] - N.loaf[s,t,] + N.recruits[s,t,]
    nestlings[s,t,] <- round(mean.fec * 0.5 * Ntot.breed[s,t,])
    JUV[s,t,] <- rpois(10000, nestlings[s,t,])
    Ntot[s,t,] <- colSums(IM[s,t,,3,]) + Ntot.breed[s,t,] + N.atsea[s,t,] + N.loaf[s,t,] 
  }
}

# SAVE RESULTS ####

save(future.mit.STATUSQUO, future.mit.INCREASE, future.mit.DECREASE,
     IM, N.recruits, N.ad.surv, N.breed.ready, N.atsea, N.loaf, 
     Ntot.breed, nestlings, JUV, Ntot, 
     file = "PVA_results.RData"
     )

# TEST PLOT #####

test_Ntot <- apply(Ntot, 2, rbind) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(scenario = case_when(
    as.numeric(rowname) %% 3 == 1 ~ 1,
    as.numeric(rowname) %% 3 == 2 ~ 2,
    as.numeric(rowname) %% 3 == 0 ~ 3
  )) %>% 
  select(-rowname) %>% 
  pivot_longer(-scenario) %>% 
  mutate(name = str_first_number(name) + 2020) %>% 
  group_by(scenario, name) %>% 
  mutate(index = row_number()) %>% 
  summarise(mean = mean(value), 
            lower = quantile(value, 0.025), 
            upper = quantile(value, 0.975))

ggplot(test_Ntot) +
  geom_ribbon(aes(x=name,
                  ymin=lower,
                  ymax=upper,
                  fill=as.character(scenario)),
              alpha=0.2) +
  geom_line(aes(x=name,
                y=mean,
                color=as.character(scenario))) +
  labs(x="",y="",title = "") 

