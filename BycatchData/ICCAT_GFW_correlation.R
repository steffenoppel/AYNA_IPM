############################################################################
######## FISHING EFFORT DATA PREPARATION FOR INTEGRATED POPULATION MODEL ###
############################################################################

### written by Steffen Oppel in November 2021
### steffen.oppel@rspb.org.uk
library(tidyverse)
library(lubridate)
library(dplyr)
library(rgdal)
library(maptools)
library(raster)
library(sf)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select

# #########################################################################
# # LOAD FISHERY DATA FROM ICCAT (n hooks 2000 - 2017)
# #########################################################################
# ## this CSV file is based on a query that uses
#  QuadID==2 for SE Atlantic
#  QuadID==3 for SW Atlantic!
# ## see https://www.iccat.int/Data/t2ce-ENG.pdf
#### MODIFY nhooks to match coordinates and time periods
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\BycatchData"), silent=T)
fish <- read.csv("ICCAT_longline_effort_data2021.csv")
mitig <- read.csv("AYNA_bycatch_mitigation_props_perfleet.csv")
head(fish)
unique(fish$TimePeriodID) ### 1-12 signify months, 13-16 quarters, 17 is for whole year
unique(fish$QuadID)
unique(fish$Year)

## subset 2011-2019 data
fish <- subset(fish, Year > 2011)

#### SUMMARISE LONGLINE EFFORT BY GRID CELL AND YEAR QUARTER

LL_effort<- fish %>% 
  mutate(Lat=Lat*(-1), Lon=ifelse(QuadID==2,Lon,Lon*(-1))) %>% ## modify lats and longs for west and south
  mutate(quarter=if_else(TimePeriodID %in% c(1, 2, 3),"Q1",ifelse(
    TimePeriodID %in% c(4, 5, 6),"Q2",ifelse(TimePeriodID %in% c(7, 8, 9),"Q3",ifelse(TimePeriodID %in% c(10, 11, 12),"Q4","ALL"))))) %>%
  filter(Eff1Type=="NO.HOOKS") %>%
  group_by(Year, quarter,Lat,Lon) %>%
  summarise(TotEff=sum(Effort))

dim(LL_effort)

## split the annual data into quarterly data
LLeff_replace<-LL_effort %>% filter(quarter=="ALL") %>%
  mutate(quarter="Q1", TotEff=TotEff/4) %>%
  bind_rows(LL_effort %>% filter(quarter=="ALL") %>%  mutate(quarter="Q2", TotEff=TotEff/4)) %>%
  bind_rows(LL_effort %>% filter(quarter=="ALL") %>%  mutate(quarter="Q3", TotEff=TotEff/4)) %>%  
  bind_rows(LL_effort %>% filter(quarter=="ALL") %>%  mutate(quarter="Q4", TotEff=TotEff/4))

## combine the two
ICCAT<- LL_effort %>% filter(quarter!="ALL") %>% bind_rows(LLeff_replace)
dim(LL_effort)

## summary for comparison with other datasets
nhooksSummaryICCAT<-ICCAT %>%
  #filter(Year>1999 & Year<2018) %>%
  group_by(Year) %>%
  summarise(N=sum(TotEff)) %>%
  mutate(Data="ICCAT")


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


### MULTIPLY EFFORT AND DISTRIBUTION DATA AND CALCULATE FISHING OVERLAP INDEX ####
longline<-ICCAT %>% mutate(Eff=ifelse(quarter=="Q1",TotEff*AYNA1,
                                      ifelse(quarter=="Q2",TotEff*AYNA2,
                                             ifelse(quarter=="Q3",TotEff*AYNA3,TotEff*AYNA4)))) %>%
  group_by(Year) %>%
  summarise(n_hooks=sum(Eff, na.rm=T)) %>%
  full_join(mitig, by="Year")





# #########################################################################
# # LOAD FISHING EFFORT DATA FROM GLOBAL FISHING WATCH
# #########################################################################
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\BycatchData\\GFW"), silent=T)
GFWzips<-list.files(pattern="zip")
out<-data.frame()

for (fp in 1:length(GFWzips)){
  unzip(GFWzips[fp])
}

  #### read and combine data for a single year
  
  GFWcsvs<-list.files(pattern = "\\.csv$", full.names = TRUE)
  GFWcsvs<-GFWcsvs[1:365]
  
  for (f in 1:length(GFWcsvs)){
    x<-fread(GFWcsvs[f]) %>% filter(geartype %in% c("drifting_longlines","set_longlines")) %>%
      mutate(quarter=if_else(month(date) %in% c(1, 2, 3),"Q1",ifelse(
        month(date) %in% c(4, 5, 6),"Q2",ifelse(month(date) %in% c(7, 8, 9),"Q3","Q4"))))
    sfx<-st_as_sf(x, coords = c('cell_ll_lon', 'cell_ll_lat'), crs = 4326)
    if(unique(x$quarter)=="Q1"){x$AYNA<-raster::extract(AYNAQ1,sfx)}
    if(unique(x$quarter)=="Q2"){x$AYNA<-raster::extract(AYNAQ2,sfx)}
    if(unique(x$quarter)=="Q3"){x$AYNA<-raster::extract(AYNAQ3,sfx)}
    if(unique(x$quarter)=="Q4"){x$AYNA<-raster::extract(AYNAQ4,sfx)}
    xout<-x %>% filter(!is.na(AYNA))
    out<-rbind(out,as.data.frame(xout))   ### extract only fishing effort overlapping with 
    rm(x,xout)
  }
}



#### read and combine data for a single year

GFWcsvs<-list.files(pattern="csv")
GFWcsvs<-GFWcsvs[1:365]

for (f in 1:length(GFWcsvs)){
  x<-fread(GFWcsvs[f]) %>% filter(geartype %in% c("drifting_longlines","set_longlines")) %>%
    mutate(quarter=if_else(month(date) %in% c(1, 2, 3),"Q1",ifelse(
      month(date) %in% c(4, 5, 6),"Q2",ifelse(month(date) %in% c(7, 8, 9),"Q3","Q4"))))
  sfx<-st_as_sf(x, coords = c('cell_ll_lon', 'cell_ll_lat'), crs = 4326)
  if(unique(x$quarter)=="Q1"){x$AYNA<-raster::extract(AYNAQ1,sfx)}
  if(unique(x$quarter)=="Q2"){x$AYNA<-raster::extract(AYNAQ2,sfx)}
  if(unique(x$quarter)=="Q3"){x$AYNA<-raster::extract(AYNAQ3,sfx)}
  if(unique(x$quarter)=="Q4"){x$AYNA<-raster::extract(AYNAQ4,sfx)}
  xout<-x %>% filter(!is.na(AYNA))
  out<-rbind(out,as.data.frame(xout))   ### extract only fishing effort overlapping with 
  rm(x,xout)
}
head(out)












#########################################################################
### PLOT THE THREE DATASETS TOGETHER ON ONE GRAPH
#########################################################################
plotdat<-bind_rows(nhooksSummaryNAM,nhooksSummaryJR,nhooksSummaryICCAT) %>%
  group_by(Data) %>%
  mutate(scaledEff=scale(N)) %>%
  filter(Year>start)

ggplot(plotdat) +
  geom_line(aes(x=Year, y=scaledEff, col=Data), size=2) +
  #scale_x_continuous(name="Year", limits=c(1979,2019), breaks=seq(1979,2019,5)) +
  #scale_y_continuous(name="Longline fishing effort (standardized)", limits=c(-3,3), breaks=seq(-3,3,0.5)) +
  theme(panel.background=element_rect(fill="white", colour="black"),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black", angle=45, vjust=0.5),
        axis.title=element_text(size=20),
        strip.text.x=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())


