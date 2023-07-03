
# ann.fec ####

(dat$J/dat$R) %>% mean()
(dat$J/dat$R) %>% sd()

alpha <- 10
beta <- 10
(mean <- alpha / (alpha + beta))
(sd <- sqrt((alpha*beta) / ((alpha + beta)^2 * (alpha + beta + 1))))

hist(rbeta(1000, 10, 10))
summary(rbeta(1000, 10, 10))

# mean.phi.juv ####

alpha <- 75.7/7
beta <- 24.3/7
hist(rbeta(10000, alpha, beta))
summary(rbeta(10000, alpha, beta))
mean(rbeta(10000, alpha, beta))
sd(rbeta(10000, alpha, beta))

# mean.phi.ad ####

alpha <- (91 + 6)/2
beta <- (9 + -5)/2
hist(rbeta(10000, alpha, beta))
summary(rbeta(10000, alpha, beta))
mean(rbeta(10000, alpha, beta))
sd(rbeta(10000, alpha, beta))

# sigma.phi ####

alpha <- (91 + 6)/2
beta <- (9 + -5)/2
mean.phi.ad <- rbeta(10000, alpha, beta)
mu.ad <- log(mean.phi.ad / (1-mean.phi.ad)) # Logit transformation
sigma.phi <- rexp(10000, 10)
eps.phi <- rnorm(10000, 0, sigma.phi) 
phi.ad <- plogis(mu.ad + eps.phi)
hist(phi.ad)
hist(mean.phi.ad, add = T, col = 2)

alpha <- 75.7/7
beta <- 24.3/7
mean.phi.juv <- rbeta(10000, alpha, beta)
mu.juv <- log(mean.phi.juv / (1-mean.phi.juv)) # Logit transformation
sigma.phi <- rexp(10000, 10)
eps.phi <- rnorm(10000, 0, sigma.phi) 
phi.juv <- plogis(mu.juv + eps.phi)
hist(phi.juv)
hist(mean.phi.juv, add = T, col = 2)

# sigma.obs #####

apply(dat$y.count, 2, min)
apply(dat$y.count, 2, max)
apply(dat$y.count, 2, mean)
apply(const$prop.sites, 2, mean)
hist(rexp(1000, 0.1))
summary(rexp(1000, 0.1))

# mean.p.ad ####
#mean.p.ad[1] ~ dunif(0.05, 0.5)	# low detection

alpha <- 25/3
beta <- 75/3
hist(rbeta(10000, alpha, beta))
summary(rbeta(10000, alpha, beta))
mean(rbeta(10000, alpha, beta))
sd(rbeta(10000, alpha, beta))

#mean.p.ad[2] ~ dunif(0.2, 1) # high detection   
alpha <- 100/10
beta <- 50/10
hist(rbeta(10000, alpha, beta))
summary(rbeta(10000, alpha, beta))
mean(rbeta(10000, alpha, beta))
sd(rbeta(10000, alpha, beta))

# mu.p.juv ####
mu.p.juv <- rnorm(10000, -3, sd = 0.25) 
hist(plogis(mu.p.juv))

agebeta <- rnorm(10000, 0.75, sd = 0.1) # Prior for shape of increase in juvenile recapture probability with age
hist(agebeta)

sigma.p <- rexp(10000, 10)
hist(sigma.p)
eps.p <- rnorm(10000, 0, sigma.p) 
hist(eps.p)

p.juv <- matrix(NA, nrow = 10000, ncol = 37)
for (age in 1:37) {
  p.juv[, age]  <- (mu.p.juv + agebeta*max(0,age-5) + eps.p) %>% plogis()
}
p.juv <- p.juv %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Age = str_first_number(name))

ggplot(p.juv) +
  stat_eye(aes(x = Age, y = value))
  

# mu.p.ad ####
mean.p.ad <- rbeta(10000, 8.3, 25)	# low detection
mu.p.ad <- log(mean.p.ad / (1-mean.p.ad)) # Logit transformation
hist(plogis(mu.p.ad + eps.p))
hist(mean.p.ad, add = T, col = 2)

mean.p.ad <- rbeta(10000, 10, 5) # high detection     
mu.p.ad <- log(mean.p.ad / (1-mean.p.ad)) # Logit transformation
hist(plogis(mu.p.ad + eps.p))
hist(mean.p.ad, add = T, col = 2)

####### old stuff #####

hist(plogis(rnorm(1000, -4, 1)), xlim = c(0, 1), col = "black")

hist(plogis(rnorm(1000, 0, .1)), xlim = c(0, 1), col = "red")

hist(plogis(rnorm(1000, 4, 1)), xlim = c(0, 1), col = "green")

test <- function(x){
  return(-3 + (x/2-1))
}
test(ages)

par(mfrow = c(5, 3))
for(age in 1:30) {
  hist(plogis(rnorm(1000, -4, 0.5) + rnorm(1000, 1, 0.01)*age/2), 
       xlim = c(0, 1), col = "black", main = age)
}
par(mfrow = c(1, 1))





