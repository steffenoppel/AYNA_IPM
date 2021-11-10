#### CODE FROM PAQUET ET AL. 2021 ####
## https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13667
#Appendix 5: R code to compute LTRE contributions of immigration for the three study cases (and the correlation between the number of immigrants and yearly population growth rate).
#The use of LTRE contribution is recommended given the limitations of the ad hoc correlation approach (see main text).

#### kestrel with poisson immigration parameterization (code is the same for the “NoConst” parameterization)

load("zj_kestrel_pois_rev.R")
zj.kestrel.pois<-zj
nyears<-dim(zj.kestrel.pois$N)[1]
samples<- dim(zj.kestrel.pois$N)[2]
chains<- dim(zj.kestrel.pois$N)[3]

###transient LTRE contribution of immigration rate to variance in population growth rate

## calculate immigration rate

Immrate.kestrel.pois<- array(0,dim=c(nyears-1,samples,chains))

for (t in 1:(nyears-1)){
  
  Immrate.kestrel.pois[t,,]<-zj.kestrel.pois$N.imm[t+1,,]/zj.kestrel.pois$N[t,,]
}#t

#proportion of "adults" i.e. 2+ years old birds, local + immigrants
prop_ad.kestrel.pois<-(zj.kestrel.pois$N.ad+zj.kestrel.pois$N.imm)/zj.kestrel.pois$N


#calculate the arithmetic mean of each time varying demographic parameter across years (ignore the first year to limit the influence of priors on initial population size and hence on Immigration rate at year 2)
mean_fec.kestrel.pois <- matrix(0,samples,chains)
mean_phi.kestrel.pois <- matrix(0,samples,chains)
mean_phij.kestrel.pois <- matrix(0,samples,chains)
mean_Imm.kestrel.pois <- matrix(0,samples,chains)
mean_prop_ad.kestrel.pois <- matrix(0,samples,chains)

for(c in 1:chains){
  for (i in 1:samples){
    mean_fec.kestrel.pois[i,c]<-mean(zj.kestrel.pois$fec[2:(nyears-1),i,c])
    mean_phi.kestrel.pois[i,c]<-mean(zj.kestrel.pois$phi[2,2:(nyears-1),i,c])
    mean_phij.kestrel.pois[i,c]<-mean(zj.kestrel.pois$phi[1,2:(nyears-1),i,c])
    mean_Imm.kestrel.pois[i,c]<-mean(Immrate.kestrel.pois[2:(nyears-1),i,c])
    mean_prop_ad.kestrel.pois[i,c]<-mean(prop_ad.kestrel.pois[2:(nyears-1),i,c])
  }}

###transient mean sensitivities

sens_fec.kestrel.pois <- matrix(0,samples,chains)
sens_phi.kestrel.pois <- matrix(0,samples,chains)
sens_phij.kestrel.pois <- matrix(0,samples,chains)
sens_Imm.kestrel.pois <- matrix(0,samples,chains)
sens_Nad.kestrel.pois <- matrix(0,samples,chains)


for(c in 1:chains){
  for (i in 1:samples){
    sens_fec.kestrel.pois[i,c]<-mean_phij.kestrel.pois[i,c]*mean_prop_ad.kestrel.pois[i,c]+mean_phij.kestrel.pois[i,c]*(1-mean_prop_ad.kestrel.pois[i,c])*zj.kestrel.pois$prob.breedyoung[1,i,c]
    sens_phi.kestrel.pois[i,c]<-1
    sens_phij.kestrel.pois[i,c]<-mean_fec.kestrel.pois[i,c]*mean_prop_ad.kestrel.pois[i,c]+mean_fec.kestrel.pois[i,c]*(1-mean_prop_ad.kestrel.pois[i,c])*zj.kestrel.pois$prob.breedyoung[1,i,c]
    sens_Imm.kestrel.pois[i,c]<-1
    sens_Nad.kestrel.pois[i,c]<-(mean_fec.kestrel.pois[i,c]*mean_phij.kestrel.pois[i,c]+mean_phi.kestrel.pois[i,c])-(mean_fec.kestrel.pois[i,c]*mean_phij.kestrel.pois[i,c]*mean_prop_ad.kestrel.pois[i,c]+mean_phi.kestrel.pois[i,c]+mean_fec.kestrel.pois[i,c]*mean_phij.kestrel.pois[i,c]*zj.kestrel.pois$prob.breedyoung[1,i,c]*(1-mean_prop_ad.kestrel.pois[i,c]))
  }}

###transient mean contributions

cont_fec.kestrel.pois <- matrix(0,samples,chains)
cont_phi.kestrel.pois <- matrix(0,samples,chains)
cont_phij.kestrel.pois <- matrix(0,samples,chains)
cont_Imm.kestrel.pois <- matrix(0,samples,chains)
cont_Nad.kestrel.pois<- matrix(0,samples,chains)
cont_tot.kestrel.pois <- matrix(0,samples,chains)

for(c in 1:chains){
  for (i in 1:samples){
    dp_stoch <- cbind(zj.kestrel.pois$fec[2:(nyears-1),i,c],zj.kestrel.pois$phi[2,2:(nyears-1),i,c],zj.kestrel.pois$phi[1,2:(nyears-1),i,c],Immrate.kestrel.pois[2:(nyears-1),i,c],prop_ad.kestrel.pois[2:(nyears-1),i,c])
    # Derive process variance and among demographic parameters using
    # 'shrinkage' estimates of vital rates and proportionate abundances:
    dp_varcov <- var(dp_stoch)
    sensvec <- c( sens_fec.kestrel.pois[i,c], sens_phi.kestrel.pois[i,c], sens_phij.kestrel.pois[i,c],sens_Imm.kestrel.pois[i,c],sens_Nad.kestrel.pois[i,c])
    # calculate demographic contributions
    contmatrix <- matrix(0,5,5)
    
    for (k in 1:5){
      for (m in 1:5){
        contmatrix[k,m] <- dp_varcov[k,m]*sensvec[k]*sensvec[m]
      }#m
    }#k
    contributions <- rowSums(contmatrix)
    
    cont_fec.kestrel.pois[i,c] <- contributions[1]
    cont_phi.kestrel.pois[i,c] <- contributions[2]
    cont_phij.kestrel.pois[i,c] <- contributions[3]
    cont_Imm.kestrel.pois[i,c] <- contributions[4]
    cont_Nad.kestrel.pois[i,c] <- contributions[5]
    cont_tot.kestrel.pois[i,c]<-sum(contributions[])
  }}

mean(cont_Imm.kestrel.pois)
quantile(cont_Imm.kestrel.pois,c(0.025,0.975))

#mean proportion of variation in growth rate explained by variation in immigration rate (can't be properly calculated on the posteriors)
mean(cont_Imm.kestrel.pois)/mean(cont_tot.kestrel.pois)

#Note that the LTRE contribution of "the other demographic rates, as shown in Fig. 4 can be obtained by subtracting the contribution of immigration from the total contribution
mean(cont_tot.kestrel.pois-cont_Imm.kestrel.pois)
quantile(cont_tot.kestrel.pois-cont_Imm.kestrel.pois,c(0.025,0.975))

####
###correlation between population growth rate and number of Immigrants

#yearly population growth rate
lam.kestrel.pois<- array(0,dim=c(nyears-1,samples,chains))
for (t in 1:(nyears-1)){
  
  lam.kestrel.pois[t,,]<-zj.kestrel.pois$N[t+1,,]/zj.kestrel.pois$N[t,,]
  
}

#calculate correlation coefficient between growth rate and number of immigrants

cor_Nimm.kestrel.pois <- matrix(0,samples,chains)
for(c in 1:chains){
  for (i in 1:samples){
    cor_Nimm.kestrel.pois[i,c] <-cor(zj.kestrel.pois$N.imm[3:nyears,i,c],lam.kestrel.pois[2:(nyears-1),i,c])
  }
}

mean(cor_Nimm.kestrel.pois,na.rm=T)
quantile(cor_Nimm.kestrel.pois,0.025,na.rm=T)
quantile(cor_Nimm.kestrel.pois,0.975,na.rm=T)
sum(cor_Nimm.kestrel.pois[!is.na(cor_Nimm.kestrel.pois)]>0)/(chains*samples)

#save the data as available in repository
Contribution_kestrel_poisson<-list(N.imm.kestrel.pois=zj.kestrel.pois$N.imm,lam.kestrel.pois=lam.kestrel.pois,cont_Imm.kestrel.pois=cont_Imm.kestrel.pois,cont_tot.kestrel.pois=cont_tot.kestrel.pois)
save(Contribution_kestrel_poisson,file="Contribution_kestrel_poisson.RData")
