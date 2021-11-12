##########################################################################
#
# ATLANTIC YELLOW NOSED ALBATROSS SURVIVAL - PRIOR SENSIBILITY TEST
#
##########################################################################
# written by Steffen Oppel, December 2020
# motivated by paper Bannan et al 2020 about appropriate use of priors in ecology
# simulate survival and recapture probabilities from specified priors to ensure they cover plausible survival probabilities

### CONVERSION OF PRECISION SPECIFICATION
# in JAGS, precision in the normal distribution is specified by 1/variance
# in R, precision in normal distribution is specified by sqrt(variance)

precconv<-function(x){sqrt(1/x)}


### CREATE DATA FRAME WITH 10000 RANDOM VALUES DRAWN FROM PRIORS

mean.phi.ad <- runif(10000,0.7, 0.97)   # uninformative prior for all MONTHLY survival probabilities
lp.phi.ad <- log(mean.phi.ad/(1 - mean.phi.ad))    # logit transformed survival intercept
mean.phi.juv <- runif(10000,0.5, 0.9)   # uninformative prior for all MONTHLY survival probabilities
lp.phi.juv <- log(mean.phi.juv/(1 - mean.phi.juv))    # logit transformed survival intercept
sigma.surv <- runif(10000,0, 0.5)                     # Prior for standard deviation of survival
tau.surv <- sigma.surv^-2
surv.raneff <- rnorm(10000,0, precconv(tau.surv))

mean.p.ad <- runif(10000,0.2, 0.8)   # uninformative prior for all MONTHLY survival probabilities
lp.p.ad <- log(mean.p.ad/(1 - mean.p.ad))    # logit transformed survival intercept


#### SPECIFY EQUATION TO CALCULATE ADULT SURVIVAL PROBABILITY FROM PRIORS
FAKEDATA<-data.frame(lp.phi.ad,surv.raneff)
FAKEDATA %>% #full_join(INPUT, by='simul') %>%
  mutate(logit_phi=
           lp.phi.ad +      ### intercept for mean survival 
            surv.raneff) %>%
  mutate(phi=plogis(logit_phi)) %>%
  
  ggplot() + geom_histogram(aes(x=phi))


#### SPECIFY EQUATION TO CALCULATE JUVENILE SURVIVAL PROBABILITY FROM PRIORS
FAKEDATA<-data.frame(lp.phi.juv,surv.raneff)
FAKEDATA %>% #full_join(INPUT, by='simul') %>%
  mutate(logit_phi=
           lp.phi.juv +      ### intercept for mean survival 
           surv.raneff) %>%
  mutate(phi=plogis(logit_phi)) %>%
  
  ggplot() + geom_histogram(aes(x=phi))


#### SPECIFY EQUATION TO CALCULATE RECAPTURE PROBABILITY FROM PRIORS
FAKEDATA<-data.frame(lp.p.ad,surv.raneff)
FAKEDATA %>% #full_join(INPUT, by='simul') %>%
  mutate(logit_p=
           lp.p.ad +      ### intercept for mean survival 
           surv.raneff) %>%
  mutate(p=plogis(logit_p)) %>%
  
  ggplot() + geom_histogram(aes(x=p))



#### SPECIFY EQUATION TO CALCULATE AGE-SPECIFIC RETURN OF JUVENILES
agebeta<-seq(0.8,1.2, 0.01)
mu.p.juv<-seq(-5,-4,0.1)
ages<-1:15
FAKEDATA<-expand.grid(a=ages,b=agebeta,int=mu.p.juv)
FAKEDATA %>% #full_join(INPUT, by='simul') %>%
  mutate(logit_p=
           int+ a*b/2) %>%
  mutate(p=plogis(logit_p)) %>%
  
  ggplot() + geom_histogram(aes(x=p)) + facet_wrap(~a)


