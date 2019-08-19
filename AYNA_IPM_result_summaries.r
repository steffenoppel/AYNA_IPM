##########################################################################
#
# ATLANTIC YELLOW-NOSED ALBATROSS INTEGRATED POPULATION MODEL 2000-2018
#
##########################################################################
# based on output created in AYNA_futureIPM.r
# includes JAGS output from 4 scenarios of AYNA population trajectory


library(tidyverse)
library(jagsUI)
library(data.table)
#library(nimble)
filter<-dplyr::filter
select<-dplyr::select


#########################################################################
# LOAD MODEL OUTPUT FROM IPMs
#########################################################################

# Call JAGS from R
# AYNAscenario0 <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\AYNA_IPM_projection_scenario_0.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)
# AYNAscenarioM <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\AYNA_IPM_projection_scenarioM.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)
# AYNAscenarioB <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\AYNA_IPM_projection_scenarioB.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)
# AYNAscenarioMB <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\AYNA_IPM_projection_scenarioMB.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM")
load("AYNA_IPM_v5_converged_final.RData")





#########################################################################
# PRODUCE OUTPUT TABLES THAT COMBINE ALL 4 SCENARIOS
#########################################################################


## write output into file ##
export0<-as.data.frame(AYNAscenario0$summary) %>% select(c(1,5,2,3,7,8)) %>%
  setNames(c('Mean', 'Median','SD','lcl', 'ucl','Rhat')) %>%
  mutate(parameter=row.names(AYNAscenario0$summary)) %>%
  mutate(parameter=ifelse(grepl("ann.surv\\[1,",parameter,perl=T,ignore.case = T)==T,"juv.survival",parameter)) %>%
  mutate(parameter=ifelse(grepl("ann.surv\\[2,",parameter,perl=T,ignore.case = T)==T,"adult.survival",parameter)) %>%
  #mutate(parameter=ifelse(grepl("beta",parameter,perl=T,ignore.case = T)==T,"mean.survival",parameter)) %>%
  mutate(Year=c(seq(2000,2048,1), ## Ntot.breed
                seq(2005,2048,1), ## NOBS
                seq(2000,2018,1), ## ann.fec
                seq(2000,2018,1), ## skip.prob
                seq(2000,2018,1), ## imm.rec
                rep(seq(2000.5,2017.5,1),each=2), ##ann.surv adult and juvenile
                seq(2000,2047,1), ## lambda and fut.lambda
                rep(NA,8))) %>%  ## non-time varying parameters
  mutate(Scenario="no management")
tail(export0)

## check which parameters have not converged
hist(export0$Rhat)
export0 %>% filter(Rhat>1.1)


exportM<-as.data.frame(AYNAscenarioM$summary) %>% select(c(1,5,2,3,7,8)) %>%
  setNames(c('Mean', 'Median','SD','lcl', 'ucl','Rhat')) %>%
  mutate(parameter=row.names(AYNAscenario0$summary)) %>%
  mutate(parameter=ifelse(grepl("1,",parameter,perl=T,ignore.case = T)==T,"juv.survival",parameter)) %>%
  mutate(parameter=ifelse(grepl("2,",parameter,perl=T,ignore.case = T)==T,"adult.survival",parameter)) %>%
  mutate(parameter=ifelse(grepl("beta",parameter,perl=T,ignore.case = T)==T,"mean.survival",parameter)) %>%
  mutate(Year=c(seq(2000,2048,1), ## Ntot.breed
                seq(2005,2048,1), ## NOBS
                seq(2000,2018,1), ## ann.fec
                seq(2000,2018,1), ## skip.prob
                seq(2000,2018,1), ## imm.rec
                rep(seq(2000.5,2017.5,1),each=2), ##ann.surv adult and juvenile
                seq(2000,2047,1), ## lambda and fut.lambda
                rep(NA,8))) %>%  ## non-time varying parameters
  mutate(Scenario="mouse eradication")
tail(exportM)


## check which parameters have not converged
hist(exportM$Rhat)
exportM %>% filter(Rhat>1.1)



#write.table(export0,"AYNA_Gough_IPM_estimates_v2.csv", sep=",", row.names=F)



#########################################################################
# PRODUCE TABLE 1 THAT SUMMARISES DEMOGRAPHIC RATES
#########################################################################
setwd("C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\AYNA_IPM")
Table1<-export0 %>% mutate(parameter=ifelse(grepl("ann.fec",parameter,perl=T,ignore.case = T)==T,"fecundity",parameter))%>%
    mutate(parameter=ifelse(grepl("imm.rec",parameter,perl=T,ignore.case = T)==T,"recruitment",parameter))%>%
    mutate(parameter=ifelse(grepl("Ntot.breed",parameter,perl=T,ignore.case = T)==T,"pop.size",parameter)) %>%
  mutate(parameter=ifelse(grepl("fut.lambda",parameter,perl=T,ignore.case = T)==T,"future.growth.rate",parameter)) %>%
  mutate(parameter=ifelse(grepl("lambda",parameter,perl=T,ignore.case = T)==T,"population.growth.rate",parameter)) %>%
  group_by(parameter) %>%
  summarise(median=median(Median),lcl=median(lcl),ucl=median(ucl)) %>%
  filter(parameter %in% c("mean.survival","mean.fec","fecundity","mean.p","mean.skip","mean.rec","recruitment","population.growth.rate","future.growth.rate","pop.size"))
Table1
fwrite(Table1,"C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\AYNA_IPM\\Table1.csv")




#########################################################################
# PRODUCE OUTPUT GRAPH THAT SHOWS ESTIMATES FOR POPULATION TREND
#########################################################################
pdf("AYNA_IPM_pop_trend_Gough_2000_2048_v5.pdf", width=14, height=8)

## COLLATE DATA FOR ALL 4 SCENARIOS
rbind(export0,exportM[exportM$Year>2017,]) %>%
  #export0 %>%
  filter(grepl("Ntot.breed",parameter,perl=T,ignore.case = T)) %>%
  arrange(Scenario, Year) %>%


## CREATE PLOT FOR POP TREND AND SAVE AS PDF

  ggplot() + 
  geom_line(aes(y=Median, x=Year, color=Scenario), size=1)+   #
  geom_ribbon(aes(x=Year, ymin=lcl,ymax=ucl, fill=Scenario),alpha=0.5)+ #
  geom_point(data=AYNA.pop,aes(y=tot, x=Year),col="firebrick", size=2.5)+
  #geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  #facet_wrap(~Scenario, ncol=2) +
  ylab("AYNA pairs in Gough study areas") +
  scale_y_continuous(breaks=seq(0,3000,500), limits=c(0,3000))+
  scale_x_continuous(breaks=seq(2000,2040,5),limits=c(2000,2040))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())


dev.off()







#########################################################################
# PRODUCE OUTPUT GRAPH FOR SURVIVAL AND FECUNDITY
#########################################################################


## CREATE PLOT FOR SURVIVAL AND SAVE AS PDF
pdf("AYNA_IPM_survival_Gough_2000_2018.pdf", width=11, height=8)
export0 %>% filter(grepl("survival",parameter,perl=T,ignore.case = T)) %>%

  ggplot(aes(y=Median, x=Year, colour=parameter)) + geom_point(size=2.5)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  geom_hline(yintercept=export0$Median[export0$parameter=="mean.survival"][1], colour="#619CFF") +
  geom_hline(yintercept=export0$Median[export0$parameter=="mean.survival"][2], colour="#619CFF") +
  ylab("Apparent annual survival probability") +
  scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0,1))+
  scale_x_continuous(breaks=seq(2000,2020,2))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()





## CREATE PLOT FOR FECUNDITY AND SAVE AS PDF
pdf("AYNA_IPM_fecundity_Gough_2000_2018.pdf", width=11, height=8)
export0 %>% filter(grepl("ann.fec",parameter,perl=T,ignore.case = T)) %>%
  
  ggplot(aes(y=Median, x=Year)) + geom_point(size=2.5)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  geom_hline(yintercept=export0$Median[export0$parameter=="mean.fec"], colour="#619CFF") +
  ylab("Annual fecundity") +
  scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0,1))+
  scale_x_continuous(breaks=seq(2000,2020,2))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUANTIFY PROBABILITY OF EXTINCTION OVER TIME
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### no more extinction in v5b when model converged!

samplesout0<-as.data.frame(rbind(AYNAscenario0$samples[[1]],AYNAscenario0$samples[[2]],AYNAscenario0$samples[[3]],AYNAscenario0$samples[[4]])) %>%
  #select(-pop.growth.rate,-future.growth.rate,-mean.fec,-bycatch,-deviance,-mean.survival) %>%
  gather(key="parameter", value="value") %>%
  mutate(Scenario="no management")
samplesoutM<-as.data.frame(rbind(AYNAscenarioM$samples[[1]],AYNAscenarioM$samples[[2]],AYNAscenarioM$samples[[3]],AYNAscenarioM$samples[[4]])) %>%
  #select(-pop.growth.rate,-future.growth.rate,-mean.fec,-bycatch,-deviance,-`beta[1]`,-`beta[2]`) %>%
  gather(key="parameter", value="value") %>%
  mutate(Scenario="mouse eradication")


extprop <- rbind(samplesout0,samplesoutM) %>%
  mutate(Year=sub(".*\\[(.*)\\].*", "\\1", parameter, perl=TRUE)) %>%
  #mutate(Year=ifelse(grepl(",",Year),as.numeric(strsplit(Year,",", perl=TRUE)[2]),as.numeric(Year))) %>%
  filter(grepl("Ntot.breed",parameter)) %>%
  mutate(n=1, inc=ifelse(value<10,1,0)) %>%
  group_by(Scenario,Year) %>%
  summarise(ext.prob=sum(inc)/sum(n)) %>%
  mutate(Year=as.numeric(Year)+1999) 



### produce plot with multiple lines per year

#pdf("EV_extinction_probability_C3.pdf", width=10, height=7)
#postscript("Fig1_Balkan.eps", width=9, height=6)
#jpeg("Fig1_Balkan.jpg", width=9, height=6, units="in", res=600, quality=100)
#par(oma=c(0,0,0,0),mar=c(4.2,4.5,0,0.5), cex=1.2)
ggplot(data=extprop)+
  geom_line(aes(x=Year, y=ext.prob, color=Scenario), size=1)+

  ## format axis ticks
  scale_y_continuous(name="Probability of extinction (%)", limits=c(0,0.8),breaks=seq(0,0.8,0.2), labels=as.character(seq(0,80,20)))+
  scale_x_continuous(name="Year", breaks=seq(2000,2048,5), labels=as.character(seq(2000,2048,5)))+
  guides(color=guide_legend(title="Scenario"),fill=guide_legend(title="Scenario"))+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=12, color="black",angle=45, vjust = 1, hjust=1), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))

dev.off()






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CORRELATION OF POP GROWTH RATE WITH DEMOGRAPHIC PARAMETERS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



##################  ####################################
pdf("YNAL_lambda_correlations.pdf", width=9, height=9)
l.fitted<-l.lower<-l.upper<-ad.fitted<-ad.lower<-ad.upper<-ju.fitted<-ju.lower<-ju.upper<-pr.fitted<-pr.lower<-pr.upper<-bp.fitted<-bp.lower<-bp.upper<-rec.fitted<-rec.lower<-rec.upper<-numeric()
year<-c(2000:2018)
n.years<-length(year)
for (i in 1:(n.years-1)){
  l.fitted[i]<-quantile(AYNAscenario0$sims.list$lambda[,i], 0.5)
  l.lower[i]<-quantile(AYNAscenario0$sims.list$lambda[,i], 0.025)
  l.upper[i]<-quantile(AYNAscenario0$sims.list$lambda[,i], 0.975)
  
  pr.fitted[i]<-quantile(AYNAscenario0$sims.list$ann.fec[,i], 0.5)
  pr.lower[i]<-quantile(AYNAscenario0$sims.list$ann.fec[,i], 0.025)
  pr.upper[i]<-quantile(AYNAscenario0$sims.list$ann.fec[,i], 0.975)
  
  bp.fitted[i]<-quantile(AYNAscenario0$sims.list$skip.prob[,i], 0.5)
  bp.lower[i]<-quantile(AYNAscenario0$sims.list$skip.prob[,i], 0.025)
  bp.upper[i]<-quantile(AYNAscenario0$sims.list$skip.prob[,i], 0.975)
  
  rec.fitted[i]<-quantile(AYNAscenario0$sims.list$imm.rec[,i], 0.5)
  rec.lower[i]<-quantile(AYNAscenario0$sims.list$imm.rec[,i], 0.025)
  rec.upper[i]<-quantile(AYNAscenario0$sims.list$imm.rec[,i], 0.975)
  
  ad.fitted[i]<-quantile(AYNAscenario0$sims.list$ann.surv[,2,i], 0.5)
  ad.lower[i]<-quantile(AYNAscenario0$sims.list$ann.surv[,2,i], 0.025)
  ad.upper[i]<-quantile(AYNAscenario0$sims.list$ann.surv[,2,i], 0.975)
  
  ju.fitted[i]<-quantile(AYNAscenario0$sims.list$ann.surv[,1,i], 0.5)
  ju.lower[i]<-quantile(AYNAscenario0$sims.list$ann.surv[,1,i], 0.025)
  ju.upper[i]<-quantile(AYNAscenario0$sims.list$ann.surv[,1,i], 0.975)
  
  }


par(mfrow=c(2,2))

plot(l.fitted~ad.fitted, xlim=c(0.5,1), ylim=c(0.5,1.5), xlab="Adult survival probability",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(ad.lower,l.fitted,ad.upper,l.fitted ,col="gray", lty=1, lwd=0.5)
segments(ad.fitted,l.lower,ad.fitted,l.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted,ad.fitted,alternative = c("two.sided"),method = "spearman")
text(0.5,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

plot(l.fitted~ju.fitted, xlim=c(0.5,1), ylim=c(0.5,1.5), xlab="Juvenile survival probability",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(ju.lower,l.fitted,ju.upper,l.fitted ,col="gray", lty=1, lwd=0.5)
segments(ju.fitted,l.lower,ju.fitted,l.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted,ju.fitted,alternative = c("two.sided"),method = "spearman")
text(0.5,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

plot(l.fitted~pr.fitted, xlim=c(0,1), ylim=c(0.5,1.5), xlab="Annual productivity",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(pr.lower,l.fitted,pr.upper,l.fitted ,col="gray", lty=1, lwd=0.5)
segments(pr.fitted,l.lower,pr.fitted,l.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted,pr.fitted,alternative = c("two.sided"),method = "spearman")
text(0,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

plot(l.fitted[1:18]~longlineICCAT, xlim=c(-3,3), ylim=c(0.5,1.5), xlab="Fishing effort",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(longlineICCAT,l.lower[1:18],longlineICCAT,l.upper[1:18],col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted[1:18],longlineICCAT,alternative = c("two.sided"),method = "spearman")
text(0,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CORRELATION OF DEMOGRAPHIC PARAMETERS WITH FISHING EFFORT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf("YNAL_ICCAT_correlations.pdf", width=9, height=9)
par(mfrow=c(2,2))

plot(longlineICCAT~ad.fitted[1:18], xlim=c(0.5,1), ylim=c(-3,3), xlab="Adult survival probability",ylab="Fishing effort",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(ad.lower[1:18],longlineICCAT,ad.upper[1:18],longlineICCAT ,col="gray", lty=1, lwd=0.5)
test<-cor.test(longlineICCAT,ad.fitted[1:18],alternative = c("two.sided"),method = "spearman")
text(0.5,-2.5, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

plot(longlineICCAT~ju.fitted[1:18], xlim=c(0.5,1), ylim=c(-3,3), xlab="Juvenile survival probability",ylab="Fishing effort",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(ju.lower[1:18],longlineICCAT,ju.upper[1:18],longlineICCAT ,col="gray", lty=1, lwd=0.5)
test<-cor.test(longlineICCAT,ju.fitted[1:18],alternative = c("two.sided"),method = "spearman")
text(0.5,-2.5, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

plot(longlineICCAT~pr.fitted[1:18], xlim=c(0,1), ylim=c(-3,3), xlab="Annual productivity",ylab="Fishing effort",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(pr.lower[1:18],longlineICCAT,pr.upper[1:18],longlineICCAT ,col="gray", lty=1, lwd=0.5)
test<-cor.test(longlineICCAT,pr.fitted[1:18],alternative = c("two.sided"),method = "spearman")
text(0,-2.5, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

plot(longlineICCAT~l.fitted[1:18], ylim=c(-3,3), xlim=c(0.5,1.5), xlab="Population growth rate",ylab="Fishing effort",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(l.lower[1:18],longlineICCAT,l.upper[1:18],longlineICCAT,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted[1:18],longlineICCAT,alternative = c("two.sided"),method = "spearman")
text(0.5,-2.5, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

dev.off()
















