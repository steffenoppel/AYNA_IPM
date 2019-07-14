##########################################################################
#
# ATLANTIC YELLOW-NOSED ALBATROSS FISHERIES BYCATCH EFFORT
#
##########################################################################
# CONTRIBUTING DATA TO IPM
# full analysis is in AYNA_futureIPM.r
## MAJOR REVISION 11 July 2019: after chat with Cleo Small included AYNA distribution from tracking data to create fishing overlap index
## AYNA distribution provided by Ana Carneiro

library(raster)
library(tidyverse)
library(data.table)
library(sf)
filter<-dplyr::filter
select<-dplyr::select


#########################################################################
# LOAD FISHERY DATA PROVIDED BY JOEL RICE
#########################################################################
try(setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch"), silent=T)
nhooks<-fread("SH_longline_effort_30MAY_2019.csv")
head(nhooks)
nhooksSummary<-nhooks %>%
  filter(lat5>(-40.1)) %>% filter(lat5<(-20.1)) %>%
  filter(lon5<20.1) %>% filter(lon5>(-15.1)) %>%
  filter(yy>1999) %>%
  group_by(yy) %>%
  summarise(N=sum(hooks))

## scale
longlineJR<- (nhooksSummary$N-mean(nhooksSummary$N))/sd(nhooksSummary$N)


# #########################################################################
# # LOAD FISHERY DATA FROM ICCAT (n hooks 2000 - 2017)
# #########################################################################
# ### more detailed data from Namibia do not cover the same time horizon
# 
# ## this CSV file is based on a query that uses
#  QuadID==2 for SE Atlantic
#  QuadID==3 for SW Atlantic!
# ## see https://www.iccat.int/Data/t2ce-ENG.pdf
# 
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\BycatchData"), silent=T)
nhooks<-fread("N_Hooks2000_2017.csv")
head(nhooks)


#### MODIFY nhooks to match coordinates and time periods
ICCAT<-nhooks %>% mutate(Time=c(1,1,1,2,2,2,3,3,3,4,4,4)[TimePeriodID]) %>%
  mutate(Lat=-Lat, Lon=ifelse(QuadID==2,Lon,-Lon)) %>%
  rename(nhooks=SumOfEff1) %>%
  select(YearC,Time,Lat,Lon,nhooks)
head(ICCAT)
           

#### READ IN AYNA DISTRIBUTION
AYNAQ1<-raster("Atlantic Yellow-nosed Albatross_Gough_Q1.tif")
AYNAQ2<-raster("Atlantic Yellow-nosed Albatross_Gough_Q2.tif")
AYNAQ3<-raster("Atlantic Yellow-nosed Albatross_Gough_Q3.tif")
AYNAQ4<-raster("Atlantic Yellow-nosed Albatross_Gough_Q4.tif")


#### NEED TO SPATIALLY EXTRACT FISH EFFORT DATA FROM THESE RASTERS

sfICCAT<-st_as_sf(ICCAT, coords = c('Lon', 'Lat'), crs = 4326)
spICCAT<-as(sfICCAT,'Spatial')
ICCAT$AYNA1<-raster::extract(AYNAQ1,spICCAT)
ICCAT$AYNA2<-raster::extract(AYNAQ2,spICCAT)
ICCAT$AYNA3<-raster::extract(AYNAQ3,spICCAT)
ICCAT$AYNA4<-raster::extract(AYNAQ4,spICCAT)

# AYNAQ1<-rasterToPoints(AYNAQ1)
# AYNAQ1<-as.data.frame(AYNAQ1)
# names(AYNAQ1)<-c('Lon','Lat','AYNA')
# AYNAQ1$Time<-1
# AYNAQ2<-rasterToPoints(AYNAQ2)
# AYNAQ2<-as.data.frame(AYNAQ2)
# names(AYNAQ2)<-c('Lon','Lat','AYNA')
# AYNAQ2$Time<-2
# AYNAQ3<-rasterToPoints(AYNAQ3)
# AYNAQ3<-as.data.frame(AYNAQ3)
# names(AYNAQ3)<-c('Lon','Lat','AYNA')
# AYNAQ3$Time<-3
# AYNAQ4<-rasterToPoints(AYNAQ4)
# AYNAQ4<-as.data.frame(AYNAQ4)
# names(AYNAQ4)<-c('Lon','Lat','AYNA')
# AYNAQ4$Time<-4
# 
# AYNAdis<-rbind(AYNAQ1,AYNAQ2,AYNAQ3,AYNAQ4)
# head(AYNAdis)
# 


### MULTIPLY EFFORT AND DISTRIBUTION DATA AND CALCULATE FISHING OVERLAP INDEX ####
longline<-ICCAT %>% mutate(Eff=ifelse(Time==1,nhooks*AYNA1,
                            ifelse(Time==2,nhooks*AYNA2,
                                   ifelse(Time==3,nhooks*AYNA3,nhooks*AYNA4)))) %>%
  group_by(YearC) %>%
  summarise(EFF=sum(Eff, na.rm=T)) %>%
  mutate(mitigation=c(rep(1,13),0.8,0.6,0.4,0.2,0.1)) %>%   ### insert proportion of ships not using any mitigation measures
  mutate(MitEFF=EFF*mitigation)

fwrite(longline,"ICCAT_AYNA_overlay_nhooks_2000_2017.csv")


#### FIXED WITH QuadID==2
### format coordinates - these are unbelievably not specified as N or S but 'QuadID' indicates hemisphere
### extract data from 20-40 S and 20E to 10W

# nhooksSummary<-nhooks %>% #mutate(Lat=ifelse(grepl('n',DSetTypeID)==T,Lat,Lat*-1)) %>%
#   #mutate(Lon=ifelse(grepl('w',DSetTypeID)==T,Lon*-1,Lon)) %>%
#   filter(Lat<(40.1)) %>% filter(Lat>(19.9)) %>%
#   filter(Lon<20.1) %>%
#   group_by(YearC) %>%
#   summarise(N=sum(SumOfEff1))

## scale
longlineICCAT<- (longline$EFF-mean(longline$EFF))/sd(longline$EFF)


### load long-line fishing effort from Namibia
### this only goes back to 2009, but in some years is greater than what ICCAT report for SE Atlantic!!

NamLLeff<-fread("Longline_effort_Namibia.csv")
head(NamLLeff)

NamLLeffSummary<-NamLLeff %>%
  group_by(Year) %>%
  summarise(N=sum(`Number of HOOKS_SET`))

NamLLeffSummary



### load trawling effort from Namibia
### this combines both wet and frozen fish

NamTRwet<-fread("Trawl_effort_Namibia_wet.csv")
NamTRfrozen<-fread("Trawl_effort_Namibia_frozen.csv")
head(NamTRfrozen)
head(NamTRwet)
names(NamTRfrozen)<-names(NamTRwet)

NamTrawlSummary<-NamTRwet %>% bind_rows(NamTRfrozen) %>%
  group_by(YEAR) %>%
  summarise(N=sum(`DURATION(HOURS)`))
NamTrawlSummary







#########################################################################
# LOAD FISHERY DATA PROVIDED BY NAMIBIAN FISHERIES DEPARTMENT
#########################################################################
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\BycatchData"), silent=T)
nhooks<-fread("Namibia_DemersalLongLine.csv")
head(nhooks)
nhooksSummary<-nhooks %>%
  filter(Year>1999 & Year<2018) %>%
  group_by(Year) %>%
  summarise(N=sum(HOOKS_SET))

## scale
longlineNAM<- (nhooksSummary$N-mean(nhooksSummary$N))/sd(nhooksSummary$N)






### PLOT THE THREE DATASETS TOGETHER ON ONE GRAPH
plotdat<- data.frame(Year=rep(seq(2000,2017), 3), N=c(longlineJR, longlineICCAT, longlineNAM), Source=rep(c("Joel Rice","ICCAT","Namibia"), each=18))

ggplot(data=plotdat) +
  geom_line(aes(x=Year, y=N, col=Source), size=2) +
  scale_x_continuous(name="Year", limits=c(2000,2018), breaks=seq(2000,2018,2)) +
  scale_y_continuous(name="Longline fishing effort (standardized)", limits=c(-3,3), breaks=seq(-3,3,0.5)) +
  theme(panel.background=element_rect(fill="white", colour="black"),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black", angle=45, vjust=0.5),
        axis.title=element_text(size=20),
        strip.text.x=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())



### ESTIMATE CORRELATION WITH SURVIVAL ESTIMATES
## THERE IS A REASONABLE POSITIVE CORRELATION

try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM"), silent=T)
AYNA<-fread("AYNA_Gough_IPM_estimates.csv")
head(AYNA)
surv<-AYNA %>% filter(parameter=="adult.survival") %>%
  filter(Year!='mean') %>%
  select(Year,Mean,Median)

cor.test(surv$Median,log(nhooksSummary$N))
plot(surv$Median~log(nhooksSummary$N))
#cor.test(surv$Median,lag(nhooksSummary$N,2))
#cor.test(surv$Median,lead(nhooksSummary$N,2))
#cor.test(surv$Median[8:17],NamLLeffSummary$N)


# surv<-AYNA %>% filter(parameter=="adult.survival") %>%
#   filter(Year %in% c('2009','2010','2016','2017')) %>%
#   select(Year,Mean,Median)
#
# cor.test(surv$Mean,NamTrawlSummary$N)




# #########################################################################
# # INCORPORATE BYCATCH MITIGATION ADOPTION
# #########################################################################
#
# mitigation=c(rep(1,13),0.8,0.6,0.4,0.2,0.1)
# nhooksSummary$Neff<-nhooksSummary$N*mitigation
#
# ## scale
# longlineNAM<- (nhooksSummary$Neff-mean(nhooksSummary$Neff))/sd(nhooksSummary$Neff)

