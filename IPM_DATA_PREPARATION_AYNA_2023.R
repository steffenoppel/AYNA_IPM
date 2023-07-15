############################################################################
######## DATA PREPARATION FOR INTEGRATED POPULATION MODEL     ##############
############################################################################

### written by Steffen Oppel in November 2021 - adapted from TRAL code https://github.com/steffenoppel/TRAL_IPM
### steffen.oppel@rspb.org.uk
### uses existing CMR and breeding database to extract data

## revised on 14 July 2023 after chat with Sarah and Abby
## include all breeding states

library(tidyverse)
library(lubridate)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select



#############################################################################
##   1. SPECIFY THE SPECIES AND START YEAR FOR WHICH YOU WANT A SUMMARY ####
#############################################################################
## SPECIFY THE SPECIES AND START YEAR FOR SURVIVAL MODEL
SP<-"AYNA"


#############################################################################
##   5. PREPARE THE MARK-RECAPTURE DATA FOR SURVIVAL ANALYSIS ###############
#############################################################################

## run the RODBC import of CMR data in a 32-bit version of R
#system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\RODBC_CMR_import_AYNA.R")), wait = TRUE, invisible = FALSE, intern = T)
system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\RODBC_CMR_import_AYNA.R")), wait = TRUE, invisible = FALSE, intern = T)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM"), silent=T)
try(setwd("C:/Users/sop/Documents/Steffen/RSPB/Gough"), silent=T)
load("GOUGH_seabird_CMR_data.RData")


### COPIED FROM C:\STEFFEN\RSPB\UKOT\Gough\ANALYSIS\SeabirdSurvival\TRAL_survival_marray.r

## filter data for the selected species
contacts<-contacts %>% filter(SpeciesCode==SP) %>% ## %>% filter(Location %in% c("Area 1","Not Specified")) - filter by location optional
  mutate(Contact_Year=ifelse(!is.na(Date_Time) & (Contact_Year != year(Date_Time)),year(Date_Time),Contact_Year)) %>% ## fix years that were initially used as 'season'
  mutate(Contact_Year=ifelse(((Age %in% c("Chick","Fledgling") & is.na(Date_Time) & (Contact_Year == as.integer(substr(Contact_Season,1,4))))),  #
                             as.integer(Contact_Year)+1,as.integer(Contact_Year))) ## fix years that were initially used as 'season'


ages<-ages %>% filter(SpeciesCode==SP)
bands<-bands %>% filter(SpeciesCode==SP)

head(contacts)  ## CMR data
dim(contacts)
unique(contacts$Contact_Year)
unique(contacts$Location)
unique(contacts$Breeding_StatusID)

contacts %>% filter(BirdID=="GO-16-18-556")


#############################################################################
##   6. AGE ASSIGNMENT OF BIRDS FOR SURVIVAL ANALYSIS ###############
#############################################################################

### EXTRACT AGE AT DEPLOYMENT FROM DATABASE
deploy_age<-contacts %>% arrange(BirdID, Date_Time,Contact_Year) %>%
  mutate(AGE=ifelse(Age %in% c("Chick","Fledgling"),0,1)) %>%
  arrange(BirdID,Contact_Year,Date_Time) %>%
  group_by(BirdID) %>%
  summarise(MIN_AGE=min(AGE,na.rm=T), MAX_AGE=max(AGE,na.rm=T), FIRST_AGE=first(Age,na.rm=T), FIRST_Date=first(Date_Time,na.rm=T), FIRST_YEAR=min(Contact_Year,na.rm=T)) %>% #filter(is.na(FIRST_YEAR))
  mutate(FIRST_AGE=ifelse(FIRST_AGE=="Unknown" & month(FIRST_Date)>6,"Adult", as.character(FIRST_AGE))) %>%  ### unknowns marked after June were not chicks
  mutate(FIRST_YEAR=ifelse(is.na(FIRST_YEAR),year(FIRST_Date),FIRST_YEAR))
head(deploy_age)
dim(deploy_age)
unique(deploy_age$FIRST_YEAR)


MISSAGE<-deploy_age %>% filter(is.na(FIRST_AGE)) %>%   left_join(bands, by="BirdID") %>%
  select(BirdID, Band_Number,MIN_AGE,FIRST_Date,FIRST_YEAR)
dim(MISSAGE)


### DEFINE A YEAR (=SEASON) FROM SEPT X to JUNE X+1 because birds breed during summer from Sept - April
## use Contact_Season as encounter occasion grouping variable

contacts %>% #filter(is.na(Date_Time)) %>%
  #filter(is.na(Contact_Year)) %>%
  filter(is.na(Contact_Season)) 


contacts<-contacts %>%
  mutate(Contact_Season=if_else(is.na(Contact_Season),if_else(Age=="Chick",paste(Contact_Year-1,"-",substr(Contact_Year,3,4), sep =""),
                                                              paste(Contact_Year,"-",as.integer(substr(Contact_Year,3,4))+1, sep = "")),
                                Contact_Season))
dim(contacts)
head(contacts)
sort(unique(contacts$Contact_Season))




### ASSIGN AGE TO BIRDS WHERE THIS IS NOT SPECIFIED
## include a column with continuous age 

contacts<-contacts %>%
  left_join(deploy_age, by="BirdID") %>%
  mutate(AGE=ifelse(Age=="Adult",1,ifelse(Age %in% c("Chick","Fledgling"),0,NA))) %>%    ### certain assignments based on provided age
  mutate(AGE=ifelse(is.na(AGE),ifelse(FIRST_AGE=="Adult",1,if_else(Contact_Year>FIRST_YEAR,1,0)),AGE)) %>%    ### certain assignments based on provided age
  mutate(AGE=ifelse(is.na(AGE), ifelse(Sex %in% c("Male","Female"),1,NA),AGE)) %>%       ### inferred assignment from sex info - only adults can be sexed
  mutate(ContAge=ifelse(FIRST_AGE %in% c("Chick","Fledgling"),Contact_Year-FIRST_YEAR,Contact_Year-FIRST_YEAR+5)) #%>%      ### continuous age since first deployment, at least 5 years for birds marked as 'adult'

contacts %>% filter(is.na(AGE))
contacts %>% filter(is.na(ContAge))
contacts %>% filter(ContAge==1)




################################################################################################
##   7. REMOVE BIRDS FROM OUTSIDE THE STUDY AREAS AND BEFORE 1985  ###############
##############################################################################################
unique(contacts$Location)
STUDY_AREAS<- c("Area 1","Not Specified","Between the base and seal beach","Area 8","Area 3","Area 3/10","Area 10", "Prion Cave","Tumbledown")


### START IN SEASON 1982-83 - first solid data in database - but all chicks vanished, hence shifted to 1985
start<-1985  ## for CMR data

##### INDIVIDUALS CAPTURED OUTSIDE STUDY AREAS ###
OUTMARKED<-contacts %>% filter(!(Location %in% STUDY_AREAS)) %>% filter(Contact_Type=="Capture") %>% group_by(BirdID,FIRST_AGE,Location) %>%
  summarise(n_birds=length(unique(BirdID)))

# checked individuals with metal ring 883884 - they were only recorded once and can be excluded as birds from outside study area
## GO103015 (858922)can be included because it is observed regularly in styudy area
OUTMARKED<-OUTMARKED %>% filter(BirdID!="GO103015")

##### NUMBER OF DEAD RECOVERIES
contacts %>% filter(Contact_Type=="Recovery")


########## CREATE A LOOP OVER EVERY BIRD TO CHECK WHETHER THEY WERE EVER RECORDED IN STUDY AREAS
allbirds<-unique(contacts$BirdID)
fixed_contacts<-contacts %>% filter(!(BirdID %in% OUTMARKED$BirdID))
length(unique(fixed_contacts$BirdID))
length(allbirds)


### CHECK WHAT BIRDS WERE RINGED AS CHICKS BEFORE 1985
oldchicks<-deploy_age %>% filter(FIRST_AGE=="Chick") %>% filter(FIRST_YEAR<start+1)
length(unique(oldchicks))

### REMOVE RECORDS FROM BEFORE THE SET START YEAR AND BIRDS FIRST MARKED IN LAST YEAR
unique(fixed_contacts$Contact_Year)
contacts<-fixed_contacts %>%
  #filter(Contact_Year>=start) %>%
  filter(Contact_Year>= ifelse(AGE==0,start+1,start)) %>% #necessary after 1982 when chicks are ringed in year after adults (same season)
  mutate(FIRST_AGE=if_else(BirdID %in% oldchicks$BirdID,"Adult",FIRST_AGE))
dim(contacts)
unique(contacts$FIRST_AGE)

sort(unique(contacts$Contact_Season))

all.seasons <- paste(1985:2021, "-", 
                     c(86:99, 
                       "00", "01", "02", "03", "04", "05", "06", "07", "08", "09", 
                       10:22), 
                     sep = "")
all.seasons
which(!(all.seasons %in% sort(unique(contacts$Contact_Season))))
no.contact.seasons <- all.seasons[which(!(all.seasons %in% sort(unique(contacts$Contact_Season))))]
no.contact.seasons


#############################################################################
##   9. CREATE MATRIX OF ENCOUNTERS AND AGES ###############
#############################################################################
head(contacts)
table(contacts$Breeding_StatusID)

#### CONVERT BREEDING STATUS INTO descending numeric states so that the max of those states can be used
AYNA_EXPORT<- contacts %>%
  
  mutate(BreedState=ifelse(Contact_Type == "Recovery",2,
                           ifelse(!is.na(Nest_Description),6,
                           ifelse(Breeding_StatusID %in% c(1,1899636611,1899636612,1899636613,1899636615,1899636617,1899636618),6,
                                  ifelse(Breeding_StatusID %in% c(2,1899636610,1899636616,1899636619,105568723),4,
                                         ifelse(Breeding_StatusID %in% c(-1447839356,-1525788936),3,
                                                ifelse(Age=="Chick",1,3)))))))

head(AYNA_EXPORT)
dim(AYNA_EXPORT)
  # mutate(BreedState=ifelse(is.na(Nest_Description),
  #                          ifelse(Breeding_StatusID %in% c(1,1899636611,1899636612,1899636613,1899636615,1899636617,1899636618),1,2),1))
  # 

contacts %>% filter(BreedState==3) %>% summarise(cutoff=max(Contact_Year))   ## check that max year is 2018 with uncertain states


#### recode all BreedState==4 into 4 and 5 depending on whether the bird had bred before

for (l in 1:dim(AYNA_EXPORT)[1]) {
  if(AYNA_EXPORT$BreedState[l]==4) {
    id=AYNA_EXPORT$BirdID[l]
    year=AYNA_EXPORT$Contact_Year[l]
    prev_rec<-AYNA_EXPORT %>% filter(BirdID==id) %>% filter(Contact_Year<year)
    AYNA_EXPORT$BreedState[l]=ifelse(max(prev_rec$BreedState>4),5,4)
  }
}

table(AYNA_EXPORT$BreedState)

fwrite(AYNA_EXPORT,"AYNA_contacts_6states.csv")



### check birds that are in pre-breeding state
AYNA_EXPORT %>% filter(BreedState==4)

AYNA_EXPORT %>% filter(BirdID=="GO2832-TS") %>% select(Contact_Year,Date_Time,Contact_Type,Age,Location,ContAge,BreedState) %>%
  arrange(Date_Time)
