############################################################################
######## DATA PREPARATION FOR INTEGRATED POPULATION MODEL     ##############
############################################################################

### written by Steffen Oppel in December 2018
### steffen.oppel@rspb.org.uk
### uses existing CMR and breeding database to extract data

library(tidyverse)
library(lubridate)
library(data.table)



#############################################################################
##   1. SPECIFY THE SPECIES AND START YEAR FOR WHICH YOU WANT A SUMMARY ####
#############################################################################
SP<-"AYNA"
start<-2000




###################################################################################
##   2. READ IN DATA FROM DATABASES AND FILTER DATA FOR SPECIES OF INTEREST ####
###################################################################################

## run the RODBC import of nest and count data in a 32-bit version of R
system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database\\RODBC_count_import.r")), wait = TRUE, invisible = FALSE)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database"), silent=T)
load("GOUGH_seabird_data.RData")

## run the RODBC import of CMR data in a 32-bit version of R
system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\RODBC_CMR_import.R")), wait = TRUE, invisible = FALSE)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS"), silent=T)
load("GOUGH_seabird_CMR_data.RData")

## import summary file from Alex Bond (origin of data cannot be verified!!)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\AYNA_IPM"), silent=T)
backup<-fread("AYNA_summary_1982_2011.csv")

## filter data for the selected species
contacts<-contacts %>% filter(SpeciesCode==SP)
nests<-nests %>% filter(Species==SP)
counts<-counts %>% filter(Species==SP)


## look at data
head(contacts)  ## CMR data
head(nests)  ## nest monitoring data
head(counts)  ## seabird count data
head(backup)






#############################################################################
##   3. PREPARE THE BREEDING SUCCESS DATA FROM NEST RECORDS #################
#############################################################################

### remove only partially monitored nests
exclude <- nests %>% #filter(Year==2010) %>%
  filter(LastStage=="INCU") %>%
  filter(SUCCESS==1)


### summary of breeding success per year from nests
FECUND<-nests %>% filter(Species==SP) %>% mutate(count=1) %>%
  filter(!NestID %in% exclude$NestID) %>%
  filter(Year<2018) %>%
  group_by(Year) %>%
  summarise(n_nests=sum(count),BREED_SUCC=mean(SUCCESS, na.rm=T))
FECUND


### add missing years from backup data (unknown source)
FEC2<-backup %>% mutate(n_nests=NA) %>%
  mutate(BREED_SUCC=CHIC/INCU) %>%
  select(Year,n_nests,BREED_SUCC) %>%
  filter(!Year %in% FECUND$Year)       ### remove the years already covered by good data

FECUND<-rbind(FECUND,FEC2) %>% arrange(Year)
FECUND


### PLOT TO SPOT ANY OUTLIERS OF BREEDING SUCCESS
ggplot(FECUND, aes(x=Year,y=BREED_SUCC)) +geom_point(size=2, color='darkred')+geom_smooth(method='lm') 





#############################################################################
##   4. PREPARE THE POPULATION COUNT DATA FROM COUNT RECORDS ################
#############################################################################

### find years in which no 'INCU' were counted - we need 'AON' for those years
### fixed on 24 Dec 2018 in database


### summary of population counts of breeding pairs per year and colony
POPSIZE<-counts %>% filter(Species==SP) %>%
  filter(Colony %in% c("Area 1","Area 2","Area 3","Area 4","Area 5","Area 6","Area 7","Area 9","Area 10")) %>%
  mutate(Year=year(Date)) %>%
  filter(Breed_Stage=="INCU") %>%
  filter(Cohort %in% c("INCU","TERR","AON")) %>%
  #filter(Cohort %in% c("INCU","TERR")) %>%
  group_by(Year,Colony) %>%
  summarise(N=sum(Number, na.rm=T)) %>%
  spread(key=Colony, value=N)
POPSIZE[19:35,]


### add missing years from backup data (unknown source)
POP2<-backup %>% mutate(Colony="Area 1") %>%
  select(Year,Colony,INCU) %>%
  #filter(!Year %in% FECUND$Year) %>%       ### remove the years already covered by good data
  spread(key=Colony, value=INCU)







