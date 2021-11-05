############################################################################
######## IMPORT SEABIRD CMR DATA FROM DATABASE ##########
############################################################################

library(RODBC)

#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Gough\\DATA\\CMR_Database"), silent=T)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\CMR_Database"), silent=T)

## find name of most recent version ##
recversion<-list.files(pattern="Gough_CMR_BackEnd")

db <- odbcConnectAccess2007(recversion)
contacts<- sqlQuery(db, "SELECT * FROM EXPORT_AYNA_contacts_nest_info")
contactsbreed<- sqlQuery(db, "SELECT * FROM EXPORT_AYNA_contacts_breeding_status")
metalside<- sqlQuery(db, "SELECT * FROM qry_loc_of_metal_ring")
ages<- sqlQuery(db, "SELECT * FROM Deployment_Age_BirdID")
bands<-sqlQuery(db, "SELECT * FROM qry_federal_ring_per_BirdID")
odbcClose(db)

#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM"), silent=T)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM"), silent=T)
save.image("GOUGH_seabird_CMR_data.RData")
