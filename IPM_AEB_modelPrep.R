# generating initial values

library(coda)
library(tidyverse)

#### LOAD AND COMBINE RESULTS ####

# for each chain, load the resulting MCMC output
assign("out-jags-inits", readRDS("out-jags-inits.RDS"))
assign("out-jags-noinits", readRDS("out-jags-noinits.RDS"))
assign("out-jags-inits-parallel", readRDS("out-jags-inits-parallel.RDS"))
assign("out-jags-noinits-parallel", readRDS("out-jags-noinits-parallel.RDS"))

`out-jags-inits` <- `out-jags-inits`$samples
`out-jags-noinits` <- `out-jags-noinits`$samples
`out-jags-inits-parallel` <- `out-jags-inits-parallel`$samples
`out-jags-noinits-parallel` <- `out-jags-noinits-parallel`$samples

`out-jags-inits` <- lapply(`out-jags-inits`, tail, 2000)
`out-jags-noinits` <- lapply(`out-jags-noinits`, tail, 2000)
`out-jags-inits-parallel` <- lapply(`out-jags-inits-parallel`, tail, 2000)
`out-jags-noinits-parallel` <- lapply(`out-jags-noinits-parallel`, tail, 2000)

assign("out-nimblle-inits", readRDS("out-nimblle-inits.RDS"))
assign("out-nimblle-noinits", readRDS("out-nimblle-noinits.RDS"))
assign("out3.parallel-1", readRDS("out3.parallel-1.RDS"))
assign("out3.parallel-2", readRDS("out3.parallel-2.RDS"))
assign("out3.parallel-3", readRDS("out3.parallel-3.RDS"))
assign("out4.parallel-1", readRDS("out4.parallel-1.RDS"))
assign("out4.parallel-2", readRDS("out4.parallel-2.RDS"))
assign("out4.parallel-3", readRDS("out4.parallel-3.RDS"))

# combine them and make an mcmc object
assign("out-nimblle-inits-parallel", 
       list(`out3.parallel-1`$chain1, `out3.parallel-2`$chain2, `out3.parallel-1`$chain3) %>% as.mcmc.list())
assign("out-nimblle-noinits-parallel", 
       list(`out4.parallel-1`$chain1, `out4.parallel-2`$chain2, `out4.parallel-1`$chain3) %>% as.mcmc.list())

rm(`out3.parallel-1`, `out3.parallel-2`, `out3.parallel-3`, 
   `out4.parallel-1`, `out4.parallel-2`, `out4.parallel-3`)

`out-nimblle-inits` <- lapply(`out-nimblle-inits`, tail, 2000)
`out-nimblle-noinits` <- lapply(`out-nimblle-noinits`, tail, 2000)
`out-nimblle-inits-parallel` <- lapply(`out-nimblle-inits-parallel`, tail, 2000)
`out-nimblle-noinits-parallel` <- lapply(`out-nimblle-noinits-parallel`, tail, 2000)

save.image(file='jags-nimble-comparison.RData')

# now you can treat it just like you ran the chains in succession
library(coda)

gd <- gelman.diag(`out-jags-inits`, multivariate = FALSE)$psrf
gd2 <- gd[is.na(gd[, 1]) | is.infinite(gd[, 1]), ]

View(gd) # not converging
View(gd2) # but looks good here
# ??????

gd <- gelman.diag(`out-jags-noinits`, multivariate = FALSE)$psrf
gd2 <- gd[is.na(gd[, 1]) | is.infinite(gd[, 1]), ]

View(gd) # not converging
View(gd2) # but looks good here

gd <- gelman.diag(`out-jags-noinits-parallel`, multivariate = FALSE)$psrf
gd2 <- gd[is.na(gd[, 1]) | is.infinite(gd[, 1]), ]

View(gd) # not converging
View(gd2) # but looks good here

gd <- gelman.diag(`out-jags-inits-parallel`, multivariate = FALSE)$psrf
gd2 <- gd[is.na(gd[, 1]) | is.infinite(gd[, 1]), ]

View(gd) # not converging
View(gd2) # but looks good here

gd <- gelman.diag(`out-nimblle-inits`, multivariate = FALSE)$psrf
gd2 <- gd[is.na(gd[, 1]) | is.infinite(gd[, 1]), ]

View(gd) # not converging
View(gd2) # but looks good here

gd <- gelman.diag(`out-nimblle-noinits`, multivariate = FALSE)$psrf
gd2 <- gd[is.na(gd[, 1]) | is.infinite(gd[, 1]), ]

View(gd) # not converging
View(gd2) # but looks good here

gd <- gelman.diag(`out-nimblle-noinits-parallel`, multivariate = FALSE)$psrf
gd2 <- gd[is.na(gd[, 1]) | is.infinite(gd[, 1]), ]

View(gd) # not converging
View(gd2) # but looks good here

gd <- gelman.diag(`out-nimblle-inits-parallel`, multivariate = FALSE)$psrf
gd2 <- gd[is.na(gd[, 1]) | is.infinite(gd[, 1]), ]

View(gd) # not converging
View(gd2) # but looks good here

`out-jags-noinits` <- lapply(`out-jags-noinits`, function(x) x[, order(colnames(x))]) %>% 
  lapply(function(x) x[, -1])
`out-nimblle-noinits` <- lapply(`out-nimblle-noinits`, function(x) x[, order(colnames(x))])

pdf("sadplots-jags.pdf")
#plot(out)
par(mfrow = c(3, 2))
traceplot(`out-jags-noinits`)
par(mfrow = c(1,1))
dev.off()

pdf("sadplots-nimble.pdf")
#plot(out)
par(mfrow = c(3, 2))
traceplot(`out-nimblle-noinits`)
par(mfrow = c(1,1))
dev.off()

nimble.gd <- gelman.diag(`nimble-maxAge15`, multivariate = FALSE)$psrf
nimble.gd2 <- gd[is.na(gd[, 1]) | is.infinite(gd[, 1]), ]

jags.gd <- gelman.diag(`jags-maxAge15`, multivariate = FALSE)$psrf
jags.gd2 <- gd[is.na(gd[, 1]) | is.infinite(gd[, 1]), ]

##########

nimble.gd <- gelman.diag(`out-nimble-randeffage`, multivariate = FALSE)$psrf
nimble.gd2 <- gd[is.na(gd[, 1]) | is.infinite(gd[, 1]), ]

jags.gd <- gelman.diag(`out-jags-randeffage`, multivariate = FALSE)$psrf
jags.gd2 <- gd[is.na(gd[, 1]) | is.infinite(gd[, 1]), ]

`out-jags-randeffage` <- lapply(`out-jags-randeffage`, function(x) x[, order(colnames(x))]) %>% 
  lapply(function(x) x[, -1])
`out-nimble-randeffage` <- lapply(`out-nimble-randeffage`, function(x) x[, order(colnames(x))])

pdf("sadplots-jags-re.pdf")
#plot(out)
par(mfrow = c(3, 2))
traceplot(`out-jags-randeffage`)
par(mfrow = c(1,1))
dev.off()

pdf("sadplots-nimble-re.pdf")
#plot(out)
par(mfrow = c(3, 2))
traceplot(`out-nimble-randeffage`)
par(mfrow = c(1,1))
dev.off()
