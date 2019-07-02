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
  mutate(parameter=ifelse(grepl("1,",parameter,perl=T,ignore.case = T)==T,"juv.survival",parameter)) %>%
  mutate(parameter=ifelse(grepl("2,",parameter,perl=T,ignore.case = T)==T,"adult.survival",parameter)) %>%
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
# PRODUCE OUTPUT GRAPH THAT SHOWS ESTIMATES FOR POPULATION TREND
#########################################################################
pdf("AYNA_IPM_pop_trend_Gough_2000_2029_management.pdf", width=12, height=8)

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
# EXTRACT THE SAMPLES WHERE FUTURE LAMBDA IS POSITIVE 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samplesout<-as.data.frame(rbind(AYNAscenarioMB$samples[[1]],AYNAscenarioMB$samples[[2]],AYNAscenarioMB$samples[[3]],AYNAscenarioMB$samples[[4]]))
head(samplesout)

ggplot(samplesout)+
  geom_point(aes(x=mean.fec,y=future.growth.rate,width=2))+
  geom_hline(aes(yintercept=1), color='red', size=2)+   
  ylab("Future population trend") +
  xlab("Improvement in survival probability") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=16, color="black"),
        axis.text.x=element_text(size=16, color="black", vjust=0.5), 
        axis.title=element_text(size=16), 
        strip.text.x=element_text(size=16, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())


### find the cutoff for lambda>= 1
## with a confidence of 90%

samplesout %>% select(improve.surv,fut.lambda) %>%
  mutate(improve.surv=round(improve.surv,2)) %>%
  group_by(improve.surv) %>%
  summarise(lcl=quantile(fut.lambda,0.1)) %>%
  arrange(desc(lcl))






















