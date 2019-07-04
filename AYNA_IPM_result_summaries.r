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
load("AYNA_IPM_output_4scenarios.RData")





#########################################################################
# PRODUCE OUTPUT TABLES THAT COMBINE ALL 4 SCENARIOS
#########################################################################


## write output into file ##
export0<-as.data.frame(AYNAscenario0$summary) %>% select(c(1,5,2,3,7,8)) %>%
  setNames(c('Mean', 'Median','SD','lcl', 'ucl','Rhat')) %>%
  mutate(parameter=row.names(AYNAscenario0$summary)) %>%
  mutate(parameter=ifelse(grepl("ann.surv\\[1,",parameter,perl=T,ignore.case = T)==T,"juv.survival",parameter)) %>%
  mutate(parameter=ifelse(grepl("ann.surv\\[2,",parameter,perl=T,ignore.case = T)==T,"adult.survival",parameter)) %>%
  mutate(parameter=ifelse(grepl("beta",parameter,perl=T,ignore.case = T)==T,"mean.survival",parameter)) %>%
  mutate(Year=c(seq(2000,2028,1),rep(seq(2000,2018,1),3),rep(seq(2000.5,2017.5,1),each=2),rep(NA,7))) %>%
  mutate(Scenario="no management")
tail(export0)

exportM<-as.data.frame(AYNAscenarioM$summary) %>% select(c(1,5,2,3,7,8)) %>%
  setNames(c('Mean', 'Median','SD','lcl', 'ucl','Rhat')) %>%
  mutate(parameter=row.names(AYNAscenario0$summary)) %>%
  mutate(parameter=ifelse(grepl("1,",parameter,perl=T,ignore.case = T)==T,"juv.survival",parameter)) %>%
  mutate(parameter=ifelse(grepl("2,",parameter,perl=T,ignore.case = T)==T,"adult.survival",parameter)) %>%
  mutate(parameter=ifelse(grepl("beta",parameter,perl=T,ignore.case = T)==T,"mean.survival",parameter)) %>%
  mutate(Year=c(seq(2000,2028,1),rep(seq(2000,2018,1),3),rep(seq(2000.5,2017.5,1),each=2),rep(NA,7))) %>%
  mutate(Scenario="mouse eradication")
tail(exportM)

exportB<-as.data.frame(AYNAscenarioB$summary) %>% select(c(1,5,2,3,7,8)) %>%
  setNames(c('Mean', 'Median','SD','lcl', 'ucl','Rhat')) %>%
  mutate(parameter=row.names(AYNAscenario0$summary)) %>%
  mutate(parameter=ifelse(grepl("1,",parameter,perl=T,ignore.case = T)==T,"juv.survival",parameter)) %>%
  mutate(parameter=ifelse(grepl("2,",parameter,perl=T,ignore.case = T)==T,"adult.survival",parameter)) %>%
  mutate(parameter=ifelse(grepl("beta",parameter,perl=T,ignore.case = T)==T,"mean.survival",parameter)) %>%
  mutate(Year=c(seq(2000,2028,1),rep(seq(2000,2018,1),3),rep(seq(2000.5,2017.5,1),each=2),rep(NA,7))) %>%
  mutate(Scenario="bycatch mitigation")
tail(exportB)

exportMB<-as.data.frame(AYNAscenarioMB$summary) %>% select(c(1,5,2,3,7,8)) %>%
  setNames(c('Mean', 'Median','SD','lcl', 'ucl','Rhat')) %>%
  mutate(parameter=row.names(AYNAscenario0$summary)) %>%
  mutate(parameter=ifelse(grepl("1,",parameter,perl=T,ignore.case = T)==T,"juv.survival",parameter)) %>%
  mutate(parameter=ifelse(grepl("2,",parameter,perl=T,ignore.case = T)==T,"adult.survival",parameter)) %>%
  mutate(parameter=ifelse(grepl("beta",parameter,perl=T,ignore.case = T)==T,"mean.survival",parameter)) %>%
  mutate(Year=c(seq(2000,2028,1),rep(seq(2000,2018,1),3),rep(seq(2000.5,2017.5,1),each=2),rep(NA,7))) %>%
  mutate(Scenario="mouse eradication and bycatch mitigation")
tail(exportMB)




#write.table(export0,"AYNA_Gough_IPM_estimates_v2.csv", sep=",", row.names=F)



#########################################################################
# PRODUCE TABLE 1 THAT SUMMARISES DEMOGRAPHIC RATES
#########################################################################

export0 %>% mutate(parameter=ifelse(grepl("ann.fec",parameter,perl=T,ignore.case = T)==T,"fecundity",parameter))%>%
    mutate(parameter=ifelse(grepl("imm.rec",parameter,perl=T,ignore.case = T)==T,"recruitment",parameter))%>%
    mutate(parameter=ifelse(grepl("skip.",parameter,perl=T,ignore.case = T)==T,"breed.propensity",parameter))%>%
    mutate(parameter=ifelse(grepl("Ntot.breed",parameter,perl=T,ignore.case = T)==T,"pop.size",parameter)) %>%
  group_by(parameter) %>%
  summarise(median=mean(Median),lcl=mean(lcl),ucl=mean(ucl))








# samplesout0<-as.data.frame(rbind(AYNAscenario0$samples[[1]],AYNAscenario0$samples[[2]],AYNAscenario0$samples[[3]],AYNAscenario0$samples[[4]])) %>%
#   select(-pop.growth.rate,-future.growth.rate,-mean.fec,-bycatch,-deviance,-`beta[1]`,-`beta[2]`) %>%
#   gather(key="parameter", value="value") %>%
#   mutate(Scenario="no management")
# samplesoutM<-as.data.frame(rbind(AYNAscenarioM$samples[[1]],AYNAscenarioM$samples[[2]],AYNAscenarioM$samples[[3]],AYNAscenarioM$samples[[4]])) %>%
#   select(-pop.growth.rate,-future.growth.rate,-mean.fec,-bycatch,-deviance,-`beta[1]`,-`beta[2]`) %>%
#   gather(key="parameter", value="value") %>%
#   mutate(Scenario="mouse eradication")
# samplesoutMB<-as.data.frame(rbind(AYNAscenarioMB$samples[[1]],AYNAscenarioMB$samples[[2]],AYNAscenarioMB$samples[[3]],AYNAscenarioMB$samples[[4]])) %>%
#   select(-pop.growth.rate,-future.growth.rate,-mean.fec,-bycatch,-deviance,-`beta[1]`,-`beta[2]`) %>%
#   gather(key="parameter", value="value") %>%
#   mutate(Scenario="bycatch mitigation")
# samplesoutB<-as.data.frame(rbind(AYNAscenarioB$samples[[1]],AYNAscenarioB$samples[[2]],AYNAscenarioB$samples[[3]],AYNAscenarioB$samples[[4]])) %>%
#   select(-pop.growth.rate,-future.growth.rate,-mean.fec,-bycatch,-deviance,-`beta[1]`,-`beta[2]`) %>%
#   gather(key="parameter", value="value") %>%
#   mutate(Scenario="mouse eradication and bycatch mitigation")
# 
# demcorr<-rbind(samplesout0,samplesoutM,samplesoutMB,samplesoutB) %>%
#   mutate(Year=sub(".*\\[(.*)\\].*", "\\1", parameter, perl=TRUE)) %>%
#   mutate(Year=ifelse(grepl(",",Year),as.numeric(strsplit(Year,",", perl=TRUE)[2]),as.numeric(Year))) %>%
#   mutate(parameter=ifelse(grepl("ann.surv\\[1,",parameter,perl=T,ignore.case = T)==T,"juv.survival",parameter)) %>%
#   mutate(parameter=ifelse(grepl("ann.surv\\[2,",parameter,perl=T,ignore.case = T)==T,"adult.survival",parameter))%>%
#   mutate(parameter=ifelse(grepl("ann.fec",parameter,perl=T,ignore.case = T)==T,"fecundity",parameter))%>%
#   mutate(parameter=ifelse(grepl("imm.rec",parameter,perl=T,ignore.case = T)==T,"recruitment",parameter))%>%
#   mutate(parameter=ifelse(grepl("skip.",parameter,perl=T,ignore.case = T)==T,"breed.propensity",parameter))%>%
#   mutate(parameter=ifelse(grepl("Ntot.breed",parameter,perl=T,ignore.case = T)==T,"pop.size",parameter)) %>%
#   arrange(Year,parameter)
# 
# dim(demcorr)
# head(demcorr)
# unique(demcorr$parameter)
# unique(demcorr$Year)





#########################################################################
# PRODUCE OUTPUT GRAPH THAT SHOWS ESTIMATES FOR POPULATION TREND
#########################################################################
pdf("AYNA_IPM_pop_trend_Gough_2000_2029_management.pdf", width=14, height=8)

## COLLATE DATA FOR ALL 4 SCENARIOS
rbind(export0,exportM[exportM$Year>2017,],exportB[exportB$Year>2017,],exportMB[exportMB$Year>2017,]) %>% filter(grepl("Ntot.breed",parameter,perl=T,ignore.case = T)) %>%
  arrange(Scenario, Year) %>%


## CREATE PLOT FOR POP TREND AND SAVE AS PDF

  ggplot() + 
  geom_line(aes(y=Median, x=Year, color=Scenario), size=1)+
  geom_ribbon(aes(x=Year, ymin=lcl,ymax=ucl, fill=Scenario),alpha=0.2)+
  #geom_point(size=2.5)+ geom_line()+
  #geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  #facet_wrap(~Scenario, ncol=2) +
  ylab("Number of AYNA pairs in Gough study areas") +
  scale_y_continuous(breaks=seq(200,1500,100), limits=c(200,1500))+
  scale_x_continuous(breaks=seq(2000,2028,2))+
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
  geom_hline(yintercept=0.7809829, colour="#619CFF") +
  geom_hline(yintercept=0.9176491, colour="#619CFF") +
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
# CREATE HISTOGRAM OF FUTURE LAMBDA AND CONTRAST THE FOUR SCENARIOS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samplesout0<-as.data.frame(rbind(AYNAscenario0$samples[[1]],AYNAscenario0$samples[[2]],AYNAscenario0$samples[[3]],AYNAscenario0$samples[[4]])) %>% select(future.growth.rate) %>%
  mutate(Scenario="no management")
samplesoutM<-as.data.frame(rbind(AYNAscenarioM$samples[[1]],AYNAscenarioM$samples[[2]],AYNAscenarioM$samples[[3]],AYNAscenarioM$samples[[4]])) %>% select(future.growth.rate) %>%
  mutate(Scenario="mouse eradication")
samplesoutMB<-as.data.frame(rbind(AYNAscenarioMB$samples[[1]],AYNAscenarioMB$samples[[2]],AYNAscenarioMB$samples[[3]],AYNAscenarioMB$samples[[4]])) %>% select(future.growth.rate) %>%
  mutate(Scenario="bycatch mitigation")
samplesoutB<-as.data.frame(rbind(AYNAscenarioB$samples[[1]],AYNAscenarioB$samples[[2]],AYNAscenarioB$samples[[3]],AYNAscenarioB$samples[[4]])) %>% select(future.growth.rate) %>%
  mutate(Scenario="mouse eradication and bycatch mitigation")

medlam<-rbind(samplesout0,samplesoutM,samplesoutMB,samplesoutB) %>% group_by(Scenario) %>%
  summarise(lam=mean(future.growth.rate))

## CREATE HISTOGRAMS WITH MEDIAN LINES
pdf("AYNA_future_growth_rate_Gough_management.pdf", width=11, height=9)
rbind(samplesout0,samplesoutM,samplesoutMB,samplesoutB) %>%

ggplot()+
  facet_wrap(~Scenario, ncol=2) +
  geom_histogram(aes(x=future.growth.rate))+
  geom_vline(data=medlam,aes(xintercept=lam), color='firebrick', size=1)+
  xlab("Future population trend") +
  ylab("Frequency of simulations") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=16, color="black"),
        axis.text.x=element_text(size=16, color="black", vjust=0.5), 
        axis.title=element_text(size=16), 
        strip.text.x=element_text(size=16, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()






## CALCULATE QUANTILES AND LAMBDA as pop.size/lag(pop.size)















################## PLOTTING THE POP GROWTH RATE AGAINST THE DEMOGRAPHIC PARAMETERS ####################################
pdf("YNAL_lambda_correlations.pdf", width=9, height=9)
l.fitted<-l.lower<-l.upper<-ad.fitted<-ad.lower<-ad.upper<-ju.fitted<-ju.lower<-ju.upper<-pr.fitted<-pr.lower<-pr.upper<-bp.fitted<-bp.lower<-bp.upper<-rec.fitted<-rec.lower<-rec.upper<-numeric()
year<-c(2000:2018)
n.years<-length(year)
for (i in 1:(n.years-1)){
  l.fitted[i]<-quantile(AYNAscenario0$sims.list$Ntot.breed[,i+1]/AYNAscenario0$sims.list$Ntot.breed[,i], 0.5)
  l.lower[i]<-quantile(AYNAscenario0$sims.list$Ntot.breed[,i+1]/AYNAscenario0$sims.list$Ntot.breed[,i], 0.025)
  l.upper[i]<-quantile(AYNAscenario0$sims.list$Ntot.breed[,i+1]/AYNAscenario0$sims.list$Ntot.breed[,i], 0.975)
  
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


par(mfrow=c(3,2))

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

plot(l.fitted~rec.fitted, xlim=c(0,1), ylim=c(0.5,1.5), xlab="Immature recruitment",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(rec.lower,l.fitted,rec.upper,l.fitted ,col="gray", lty=1, lwd=0.5)
segments(rec.fitted,l.lower,rec.fitted,l.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted,rec.fitted,alternative = c("two.sided"),method = "spearman")
text(0,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)


plot(l.fitted~bp.fitted, xlim=c(0,1), ylim=c(0.5,1.5), xlab="Breeding propensity",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(bp.lower,l.fitted,bp.upper,l.fitted ,col="gray", lty=1, lwd=0.5)
segments(bp.fitted,l.lower,bp.fitted,l.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted,bp.fitted,alternative = c("two.sided"),method = "spearman")
text(0,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)


dev.off()























