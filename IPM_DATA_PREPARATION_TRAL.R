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
SP<-"TRAL"
start<-2000




###################################################################################
##   2. READ IN DATA FROM DATABASES AND FILTER DATA FOR SPECIES OF INTEREST ####
###################################################################################

## run the RODBC import of nest and count data in a 32-bit version of R
system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database\\RODBC_count_import.r")), wait = TRUE, invisible = FALSE)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database"), silent=T)
load("GOUGH_seabird_data.RData")

## run the RODBC import of CMR data in a 32-bit version of R
system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival\\RODBC_CMR_import.R")), wait = TRUE, invisible = FALSE)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival"), silent=T)
load("GOUGH_seabird_CMR_data.RData")

# ## import summary file from Alex Bond (origin of data cannot be verified!!)
# try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\AYNA_IPM"), silent=T)
# backup<-fread("AYNA_summary_1982_2011.csv")

## filter data for the selected species
contacts<-contacts %>% filter(SpeciesCode==SP)
nests<-nests %>% filter(Species==SP)
counts<-counts %>% filter(Species==SP)


## look at data
head(nests)  ## nest monitoring data
head(counts)  ## seabird count data
#head(backup)






#############################################################################
##   3. PREPARE THE BREEDING SUCCESS DATA FROM NEST RECORDS #################
#############################################################################
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdBreedingSuccess"), silent=T)

### REQUIRES NEST RECORDS FROM 2011 and 2013 to go into database!

### remove only partially monitored nests
exclude <- nests %>% #filter(Year==2010) %>%
  filter(LastStage=="INCU") %>%
  filter(SUCCESS==1)


### summary of breeding success per year from nests
FECUND<-nests %>% filter(Species==SP) %>% mutate(count=1) %>%
  filter(!NestID %in% exclude$NestID) %>%
  group_by(Year,Colony) %>%
  summarise(n_nests=sum(count),BREED_SUCC=mean(SUCCESS, na.rm=T))
FECUND



### PLOT TO SPOT ANY OUTLIERS OF BREEDING SUCCESS
ggplot(FECUND, aes(x=Year,y=BREED_SUCC)) +geom_point(size=2, color='darkred')+geom_smooth(method='lm') 
fwrite(FECUND,"TRAL_breed_success_2004_2017.csv")



#############################################################################
##   4. PREPARE THE POPULATION COUNT DATA FROM COUNT RECORDS ################
#############################################################################
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdBreedingSuccess"), silent=T)



### summary of population counts of breeding pairs per year and colony
POPSIZE<-counts %>% filter(Species==SP) %>%
  mutate(Colony= as.character(Colony)) %>%
  mutate(Colony= if_else(Colony=="Green Hill South","Green Hill",Colony)) %>%
  mutate(Colony= if_else(Colony=="Green Hill North","Green Hill",Colony)) %>%
  mutate(Year=year(Date)) %>%
  filter(Year>start) %>%
  filter(Breed_Stage=="INCU") %>%
  filter(Cohort %in% c("INCU","TERR","AON")) %>%
  group_by(Year,Colony) %>%
  summarise(N=sum(Number, na.rm=T)) %>%
  spread(key=Colony, value=N)
POPSIZE
fwrite(POPSIZE,"TRAL_INCU_counts_2001_2018.csv")





### summary of population counts of fledglings per year and colony
CHICKCOUNT<-counts %>% filter(Species==SP) %>%
  mutate(Colony= as.character(Colony)) %>%
  mutate(Colony= if_else(Colony=="Green Hill South","Green Hill",Colony)) %>%
  mutate(Colony= if_else(Colony=="Green Hill North","Green Hill",Colony)) %>%
  mutate(Year=year(Date)) %>%
  filter(Year>start) %>%
  filter(Breed_Stage %in% c("CHIC","FLED")) %>%
  filter(Cohort %in% c("CHIC","FLED")) %>%
  group_by(Year,Colony) %>%
  summarise(N=sum(Number, na.rm=T)) %>%
  spread(key=Colony, value=N)
CHICKCOUNT
fwrite(POPSIZE,"TRAL_CHIC_counts_2001_2018.csv")





#############################################################################
##   5. PREPARE THE MARK-RECAPTURE DATA FOR SURVIVAL ANALYSIS ###############
#############################################################################

### checked on 25 Dec 2018: CMR database does not have satisfactory level of age or breeding status assignment for vast majority of contacts
## including 'age' and 'breeding status' reduces number of contacts from ~49000 to ~6000
## including the side on which a federal band was applied reduces the AYNA contacts from 25201 to 7786
## the only AYNA ringed as 'Chick' are from the 2015-16 season
## given the gross inadequacies of past records it is safer to simply use 0/1 contacts and assume all birds are breeding
## created new query called 'metalside' which extracts the recorded side of the body for the metal ring - all birds ringed as chicks should have "L"


head(contacts)  ## CMR data
dim(contacts)

### REMOVE RECORDS FROM BEFORE THE SET START YEAR
contacts<-contacts %>%
  filter(year(Date_Time)>start) 
dim(contacts)


### FIND MISSING DATA FOR SEASON AND REPLACE BASED ON DATE
contacts<-contacts %>%
  mutate(Contact_Season=as.character(Contact_Season)) %>%
  mutate(MO=month(Date_Time)) %>%
  mutate(Season1=paste(Contact_Year,(as.numeric(substr(Contact_Year,3,4))+1),sep="-")) %>%
  mutate(Season2=paste(Contact_Year-1,substr(Contact_Year,3,4),sep="-")) %>%
  mutate(Contact_Season=if_else(is.na(Contact_Season), if_else(MO>6,Season1,Season2),Contact_Season)) %>%
  select(BirdID,Location,Contact_Season)
dim(contacts)



# ### DETERMINE STATE FOR MULTI-STATE MODEL - ABANDONED on 25 Dec 2018 because records were too incomplete
# ## State 1 = chick
# ## State 2 = non-breeder
# ## State 3 = breeder
# stateContacts<-contacts %>%
#   mutate(Age=as.character(Age)) %>%
#   mutate(Age=if_else(is.na(Age),"Adult",Age)) %>%       ### BASED ON ASSUMPTION FROM QUICK CHECK IN 2018 - only adults had age lacking
#   mutate(STATE=if_else(Age=="Chick",1,2)) %>%
#   mutate(STATE=if_else(STATE>1 & Breeding_Status %in% c("Unknown","Loafing","Holding","Building","Bird on empty nest"),2,3))
# head(stateContacts)


### CREATE SIMPLE ENCOUNTER HISTORY (0/1)
TRAL_EH<-contacts %>% select(BirdID,Contact_Season) %>%
  mutate(count=1) %>%
  group_by(BirdID,Contact_Season) %>%
  summarise(STATE=max(count)) %>%
  spread(key=Contact_Season, value=STATE, fill=0)
TRAL_EH


### TRY TO ASSIGN AGE AT FIRST MARK FROM SIDE OF BODY
head(metalside)
birdage<-metalside %>% mutate(minage=if_else(Side_Of_Body=="L",0,1)) %>%
  group_by(SpeciesCode,BirdID) %>%
  summarise(MinAge=min(minage)) %>%
  filter(SpeciesCode==SP)
dim(birdage)


### INSERT AGE INTO ENCOUNTER HISTORY
TRAL_EH$AGE<-birdage$MinAge[match(TRAL_EH$BirdID,birdage$BirdID)]

### EXPORT ENCOUNTER HISTORY
dim(TRAL_EH)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival"), silent=T)
fwrite(TRAL_EH[,c(1,21,2:20)],"TRAL_simple_encounter_history_2000_2018.csv")

