library(dplyr)
library(tidyverse)
library(visdat)
library(ggplot2)
library(stringr)
library(data.table)
library(reshape)
library(ggmap)
library(stringdist)
library(fuzzyjoin)
library(zoo)


## 1) Start with the new HI Data 
##    in this case we use the new HI Data we received via Email from Jarda (JardaNew, includes new data from 2018-19)
##    and the MiceTable from Alice (JardaTableOld, includes up to 2017 data)




#### MICE TABLE ALICE **********************************************************

#### BASE TABLE ################################################################
## get the data with genotype info and HI (2014-19)
JardaTableNew   <-  read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/EmanuelData.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
JardaTableOld   <-  read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/MiceTable_fullEimeriaInfos_2014to2017.csv")

## filter out Column called "X" and Years 2010 & 2011
JardaTableNew <- JardaTableNew[!names(JardaTableNew) == "X"]
JardaTableNew <- JardaTableNew[!JardaTableNew$Year %in% c(2010, 2011),]
JardaTableOld    <- JardaTableOld[!names(JardaTableOld) == "X"]
JardaTableOld    <- JardaTableOld[!JardaTableOld$Year %in% c(2010, 2011),]

## Check if all rows are NA and delete these rows
which(!rowSums(!is.na(JardaTableNew)))
which(!rowSums(!is.na(JardaTableOld)))


## homogenize Column names, especially Mouse_ID should be universal!
setnames(JardaTableNew,
         old = c("PIN", "X_Longit", "Y_Latit"), 
         new = c("Mouse_ID", "Longitude", "Latitude"))


## homogenize Mouse_ID Row Names, old Mouse_IDs == SK_ and new Mouse_IDs == AA_
JardaTableNew$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_", x = JardaTableNew$Mouse_ID)
JardaTableNew <- JardaTableNew %>% distinct(Mouse_ID, .keep_all = T)

JardaTableOld$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_", x = JardaTableOld$Mouse_ID)
JardaTableOld$HI_NLoci <- gsub(pattern = "HI ", replacement = "", x = JardaTableOld$HI_NLoci)
JardaTableOld$HI_NLoci <- as.integer(JardaTableOld$HI_NLoci)


# remove Embryos (no interest for parasitic studies)
JardaTableNew   <- JardaTableNew[sapply(JardaTableNew$Mouse_ID, nchar) <= 7,]
JardaTableOld   <- JardaTableOld[sapply(JardaTableOld$Mouse_ID, nchar) <= 7,]

## check if any rows are NA only and delete those
which(!rowSums(!is.na(JardaTableNew)))
which(!rowSums(!is.na(JardaTableOld)))

## introduce a combined table (JardaTable) which includes Alice's old MiceTable (JardaTableOld) 
## and new HI data (JardaTableNew)
JardaTable <- left_join(JardaTableNew, JardaTableOld)
JardaTable <- JardaTable[sapply(JardaTable$Mouse_ID, nchar) <= 7,]
which(!rowSums(!is.na(JardaTable)))


## Number of Samples from HZ_BR, --> do we have the HI per year?
## all Samples remaining should be from HZ_BR
## all Samples remaining should be equipped with an HI value
table(JardaTable$Transect)
JardaTable %>% count(is.na(HI))
JardaTable %>% count(Year)

setnames(JardaTable,
         old = c("BW", "L", "LCd"),
         new = c("Body_weight", "Body_length", "Tail_length"))



## Check Data for Duplicates
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]

## Divide huge data frame into smaller selection of 
basics <- c("Mouse_ID", "Sex", "Longitude", "Latitude", "Year")

gen.loci <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
              "Es1", "Gpd1", "Idh1", "Mpi", "Np", "Sod1", "Es1C", "Gpd1C",
              "Idh1C", "MpiC", "NpC", "Sod1C", "Zfy2", "SRY1", "Y", "HI_NLoci",
              "HI")

## columns originate from Dis2014 and what is contained in JardaTableNew
dissection.cols <- c(#"Arrival", "Wean", "Death", "Dissection_date", "DaysInLab",
                      "Body_weight", "Body_length", "Tail_length", "Spleen", "Left_testis", "Right_testis",
                     "Testes", "Lepid", "SemVes", "Sperm", "Uterus", "Ovaria", "Grav",
                     "Litters", "NN", "MM", "FF", "Protocol")

## new cols consist of abbreviations for some parasites, must be specific 
## ASP == "Aspiculuris_tetraptera" or "Aspiculuris_Syphacia"?
## SYP == "Syphacia_obvelata"
## TM ==  "Trichuris_muris" (or "Taenia_martis"?)
parasite.cols <- c("Aspiculuris_tetraptera", "Syphacia_obvelata", "Trichuris_muris",
                   "Taenia_taeniformis", "Flea", "Mix_Syphacia_Aspiculuris",
                   "Heterakis_spumosa", "Mastophorus_muris", "Hymenolepis_microstoma",
                   "Catenotaenia_pusilla", "Cysticercus", "Ectoparasites",
                   "Worms_presence", "Hymenolepis_diminiuta", "Taenia_martis",
                   "Heligmosomoides_polygurus", "Taenia", "Aspiculuris_Syphacia",
                   "Trichuris", "Heterakis", "Mastophorus")


parasite_new.cols <- c("ASP", "SYP", "TM", "Taenia.taeniformis",
                       "Aspiculuris_tetraptera", "Syphacia_obvelata", "Trichuris_muris",
                       "Taenia_taeniformis", "Flea", "Mix_Syphacia_Aspiculuris",
                       "Heterakis_spumosa", "Mastophorus_muris", "Hymenolepis_microstoma",
                       "Catenotaenia_pusilla", "Cysticercus", "Ectoparasites",
                       "Worms_presence", "Hymenolepis_diminiuta", "Taenia_martis",
                       "Heligmosomoides_polygurus", "Taenia", "Aspiculuris_Syphacia",
                       "Trichuris", "Heterakis", "Mastophorus")


oocyst.cols <- c("counter", "Feces_g", "Date_count", "N_oocysts_sq1",
                 "N_oocysts_sq2", "N_oocysts_sq3",  "N_oocysts_sq4",
                 "N_oocysts_sq5", "N_oocysts_sq6", "N_oocysts_sq7",
                 "N_oocysts_sq8", "mean_neubauer", "PBS_dil_in_mL", 
                 "OPG", "Ncells")

EqPCR.cols <- c("delta_ct_ilwe_MminusE", "delta_ct_cewe_MminusE",
                ## from 2018 on AN IMPORTANT IMPROVEMENT!!!
                "MC.Eimeria")

EimGeno.cols <- c("n18S_Seq", "COI_Seq", "ORF470_Seq", "eimeriaSpecies")

## use basics columns only in JardaTable_selected, for later comparison with JardaTable
sel_JardaTable <- JardaTable[, colnames(JardaTable)%in%c(basics), ]


#### Complement data with previous tables ######################################
## Dissection Data from 2014
Dis2014 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ14_Dissections_237-523.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
Dis2014$Mouse_ID      <- paste0(Dis2014$ID, "_", Dis2014$PIN)
Dis2014$Ectoparasites <- as.integer(Dis2014$Ectoparasites)
### homogenize column names
setnames(Dis2014, 
         old = c("BW", "L", "LCd", 
                 "ASP", "SYP", 
                 "TM","Taenia.taeniformis"),
         new = c("Body_weight", "Body_length", "Tail_length", 
                 "Aspiculuris_tetraptera", "Syphacia_obvelata", 
                 "Trichuris_muris", "Taenia_taeniformis"))

## check missing Data
vis_miss(Dis2014)
#### fill in the missing Transect Data
Dis2014$Mouse_ID[is.na(Dis2014$Transect)]
Dis2014_CZ <- Dis2014 %>%
  filter(is.na(Transect)) %>%
  mutate(Transect = "CZ")
Dis2014_Transect_known <- Dis2014 %>%
  filter(!is.na(Transect))
Dis2014 <- bind_rows(Dis2014_Transect_known, Dis2014_CZ) %>%
  arrange(Mouse_ID)



### join Basics table with 2014 dissection Data
sel_JardaTable1 <- full_join(sel_JardaTable, Dis2014, by = c("Mouse_ID", "Sex", "Longitude", "Latitude", "Year"))
#vis_miss(sel_JardaTable1)



## Dissection Data from 2015
Dis2015 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ15_Mice_Parasite.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
Dis2015$Mouse_ID <- paste0(Dis2015$ID, "_", Dis2015$PIN)
## remove transported mice
Dis2015 <- Dis2015[!is.na(Dis2015$PIN),]
setnames(Dis2015, 
         old = c("BW", "L", "LCd", 
                 "ASP", "SYP", 
                 "TM","Taenia.taeniformis"),
         new = c("Body_weight", "Body_length", "Tail_length", 
                 "Aspiculuris_tetraptera", "Syphacia_obvelata", 
                 "Trichuris_muris", "Taenia_taeniformis"))

## check missing Data
#vis_miss(Dis2015)
#### fill in the missing Transect Data
Dis2015$Mouse_ID[is.na(Dis2015$Transect)]
Dis2015_PL <- Dis2015 %>%
  filter(is.na(Transect)) %>%
  mutate(Transect = "Allo_PL")
Dis2015_DE <- Dis2015 %>%
  filter(State == "Germany")
Dis2015 <- bind_rows(Dis2015_PL, Dis2015_DE) %>%
  arrange(Mouse_ID)
rm(Dis2015_DE)
rm(Dis2015_PL)
#### fill in the missing Code Data
Dis2015$Mouse_ID[is.na(Dis2015$Code)]
Dis2015_Code_known <- Dis2015 %>%
  filter(!is.na(Code))
Dis2015_Gola <- Dis2015 %>%
  filter(PIN %in% c(3469:3471)) %>%
  mutate(Code = "GOLA48")
Dis2015_Recl <- Dis2015 %>%
  filter(PIN %in% c(3472,3473)) %>%
  mutate(Code = "RECL8")
Dis2015_Trzes <- Dis2015 %>%
  filter(PIN %in% c(3474,3475)) %>%
  mutate(Code = "TRZES17")
Dis2015 <- bind_rows(Dis2015_Code_known, Dis2015_Gola, Dis2015_Recl, Dis2015_Trzes) %>%
  arrange(Mouse_ID)
rm(Dis2015_Code_known)
rm(Dis2015_Gola)
rm(Dis2015_Recl)
rm(Dis2015_Trzes)


sel_JardaTable2 <- full_join(sel_JardaTable1, Dis2015)
#vis_miss(sel_JardaTable2)




## Dissection Data from 2016
Dis2016 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ16_Dissection_1-211.csv", na.strings=c(""," ","NA"), stringsAsFactors = F)[-c(1:2),]
Dis2016$Head.taken.       <- as.integer(Dis2016$Head.taken.)
Dis2016$Ectoparasites     <- as.integer(Dis2016$Ectoparasites)
Dis2016$Tail_length       <- as.integer(Dis2016$Tail_length)
#remove transported mice
Dis2016 <- Dis2016[!is.na(Dis2016$Mouse_ID),]

sel_JardaTable3 <- full_join(sel_JardaTable2, Dis2016)
vis_miss(sel_JardaTable3)




## Dissection Data from 2017
Dis2017 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ17_Dissections_237-523.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
Dis2017$Ectoparasites     <- as.integer(Dis2017$Ectoparasites)
Dis2017$Tail_length       <- as.integer(Dis2017$Tail_length)
Dis2017$Spleen            <- as.integer(Dis2017$Spleen)
Dis2017$Left_epididymis   <- as.character(Dis2017$Left_epididymis)
Dis2017$Embryo_left       <- as.integer(Dis2017$Embryo_left)

## Dissection Data from 2018
Dis2018 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ18_Dissections.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)

## Dissection Data from 2019
Dis2019 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ19_Dissections.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)

## Genotype Data from 2010-2019
Gen_2010_2019 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ10-19_Genotypes.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)


sel_JardaTable2 <- sel_JardaTable2 %>%
  mutate(ID = ifelse(Year < 2016, "SK",
                     ifelse(Year >= 2016, "AA",
                            ifelse(NA))))











## combine dissection data from all years
HZ14_19_Dissections <- bind_rows(Dis2014, Dis2015, Dis2016, Dis2017, Dis2018, Dis2019)
which(!rowSums(!is.na(HZ14_19_Dissections)))
which(!colSums(!is.na(HZ14_19_Dissections)))
#vis_miss(HZ14_19_Dissections)

## eliminate duplicates
HZ14_19_Dissections$Mouse_ID[duplicated(HZ14_19_Dissections$Mouse_ID)]
HZ14_19_Dissections <- HZ14_19_Dissections %>% distinct(Mouse_ID, .keep_all = T)

## cut down the Column size to selected columns
HZ14_19_Dissections <- HZ14_19_Dissections[, colnames(HZ14_19_Dissections)%in%c(basics, dissection.cols, parasite_new.cols), ]

## Merged Table
JardaTable_merged_NewData <- left_join(JardaTable, HZ14_19_Dissections)
JardaTable_merged_NewData <- merge(JardaTable, HZ14_19_Dissections, by = "Mouse_ID", all = T)
#### any way to use/find "fillGapsAfterMerge()", a function Alice created?


  # Merge & complete
  mergedMiceTable <- merge(mergedMiceTable, diss2015.2,  
                           by = c("Mouse_ID"), all = T)
  mergedMiceTable <- fillGapsAfterMerge(mergedMiceTable)
  
  # some have their transect lost... TODO LATER
  mergedMiceTable$Mouse_ID[is.na(mergedMiceTable$Transect)]
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  ###################### 2016
  diss2016 <- read.csv(paste0(pathToMyData, "Field_data/HZ16_Mice_18-07-16_dissections.csv"), 
                       na.strings=c(""," ","NA"), stringsAsFactors = F)[-c(1:2),]
  names(diss2016)[names(diss2016) %in% "ID_mouse"] <- "Mouse_ID"
  diss2016$Mouse_ID <- as.character(diss2016$Mouse_ID)
  
  ## Add worms 
  worms16 <- read.csv(paste0(pathToMyData, "Field_data/HZ16_Worms.csv"), 
                      na.strings=c(""," ","NA"), stringsAsFactors = F)[-11]
  ## rename with homogeneity
  names(worms16) <- c("Mouse_ID", "Syphacia_obvelata", "Aspiculuris_tetraptera", "Mix_Syphacia_Aspiculuris",
                      "Heterakis_spumosa", "Mastophorus_muris", "Trichuris_muris", 
                      "Hymenolepis_microstoma", "Catenotaenia_pusilla", "Cysticercus")
  
  ## merge worms and dissection table 2016
  diss2016 <- merge(diss2016, worms16, all = TRUE)
  
  # rename for homogeneity
  names(diss2016)[names(diss2016) %in% "Ectoparasites"] <- "Flea"
  diss2016$Capture <- as.Date(diss2016$Capture, "%d.%m.%Y") 
  diss2016$Dissection <- as.Date(diss2016$Dissection, "%d.%m.%Y") 
  
  # merge
  mergedMiceTable <- merge(mergedMiceTable, diss2016, 
                           by = c("Mouse_ID"), all = T)
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  mergedMiceTable <- fillGapsAfterMerge(mergedMiceTable)
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  # add missing or wrong latitude/longitude
  loc2016 <- read.csv(paste0(pathToMyData, "Field_data/Cleaned_HMHZ_2016_All.csv" ),
                      stringsAsFactors=F)
  loc2016 <- loc2016[names(loc2016) %in% c("location", "GPS.coordinates.long", "GPS.coordinates.lat")]
  names(loc2016) <- c("Code", "longitude", "latitude")
  
  mergedMiceTable <- merge(mergedMiceTable, loc2016, by = "Code", all.x = TRUE)
  
  # 1. empty lon/lat
  mergedMiceTable$Longitude[is.na(mergedMiceTable$Longitude)] = 
    mergedMiceTable$longitude[is.na(mergedMiceTable$Longitude)]
  
  mergedMiceTable$Latitude[is.na(mergedMiceTable$Latitude)] = 
    mergedMiceTable$latitude[is.na(mergedMiceTable$Latitude)]
  
  # 2. wrong lon/lat
  mergedMiceTable$Latitude[!is.na(mergedMiceTable$Latitude) &
                             !is.na(mergedMiceTable$latitude) &
                             mergedMiceTable$Latitude != mergedMiceTable$latitude] =
    mergedMiceTable$latitude[!is.na(mergedMiceTable$Latitude) &
                               !is.na(mergedMiceTable$latitude) &
                               mergedMiceTable$Latitude != mergedMiceTable$latitude]
  
  
  mergedMiceTable$Longitude[!is.na(mergedMiceTable$Longitude) &
                              !is.na(mergedMiceTable$longitude) &
                              mergedMiceTable$Longitude != mergedMiceTable$longitude] =
    mergedMiceTable$longitude[!is.na(mergedMiceTable$Longitude) &
                                !is.na(mergedMiceTable$longitude) &
                                mergedMiceTable$Longitude != mergedMiceTable$longitude]
  # 3. delete duplicated rows
  mergedMiceTable <- unique(mergedMiceTable)
  
  ## **********************************************************
  ## 2017
  diss2017 <- read.csv(paste0(pathToMyData, "Field_data/HZ17_September_Mice_Dissection.csv"),
                       na.strings=c(""," ","NA"), stringsAsFactors = F)
  
  # correction excel bullshit
  diss2017$Feces_weight <- as.numeric(as.character(diss2017$Feces_weight))
  
  diss2017$Feces_weight[diss2017$Feces_weight > 100 & !is.na(diss2017$Feces_weight)] <-
    diss2017$Feces_weight[diss2017$Feces_weight > 100 & !is.na(diss2017$Feces_weight)] / 1000
  
  names(diss2017)[names(diss2017) == "Ectoparasites"] <- "Flea"      
  
  ## Add worms 
  worms17 <- read.csv2(paste0(pathToMyData, "Field_data/HZ17_September_Mice_Dissection_Jen_final.csv"),
                       stringsAsFactors = F)
  
  names(worms17)[names(worms17) %in% "Mesocestoides"] <- "Taenia_martis"
  
  ## merge worms and dissection table 2016
  diss2017 <- merge(diss2017, worms17, 
                    by = "Mouse_ID", all = TRUE)
  diss2017 <- fillGapsAfterMerge(diss2017)
  
  # merge
  mergedMiceTable <- merge(mergedMiceTable, diss2017, 
                           by = c("Mouse_ID"), all = T)
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  mergedMiceTable <- fillGapsAfterMerge(mergedMiceTable)
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  ## Uniformisation
  mergedMiceTable$Sex[grep("female*.", mergedMiceTable$Sex)] <- "F"
  mergedMiceTable$Sex[grep("male*.", mergedMiceTable$Sex)] <- "M"
  
  # Add old HI from previous jarda table
  oldJarda <- read.csv(paste0(pathToMyData, "/Field_data/HIforEH_May2017.csv"),
                       stringsAsFactors=F)
  
  # Uniformize IDs
  oldJarda$Mouse_ID <- oldJarda$PIN
  oldJarda$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_",x = oldJarda$PIN)
  
  # Add ommited samples
  toadd <- oldJarda[!oldJarda$Mouse_ID %in% mergedMiceTable$Mouse_ID,]
  
  # Add samples without HI in mergedmicetable but with one in oldJarda
  missing <- mergedMiceTable$Mouse_ID[is.na(mergedMiceTable$HI)]
  
  toadd <- rbind(toadd, oldJarda[oldJarda$Mouse_ID %in% missing,])
  
  # mergedMiceTable[mergedMiceTable$Mouse_ID %in% missing 
  #                 & is.na(mergedMiceTable$Longitude) &
  #                   mergedMiceTable$Mouse_ID %in% toadd$Mouse_ID,]
  
  # setnames and merge
  setnames(toadd,
           old = c("BW", "L", "LCd", 
                   "Aspiculuris", "Syphacia", 
                   "Trichuris","Taenia"),
           new = c("Body_weight", "Body_length", "Tail_length", 
                   "Aspiculuris_tetraptera", "Syphacia_obvelata", 
                   "Trichuris_muris", "Taenia_taeniformis"))
  
  toadd = toadd[names(toadd)%in%c("REGion", "mtBamH", "Zfy2", "SRY1", "Y", "X332", "X347",
                                  "X65", "Tsx", "Btk", "Syap1", "Es1", "Gpd1", "Idh1",
                                  "Mpi", "Np", "Sod1", "Es1C", "Gpd1C", "Idh1C", "MpiC", "NpC",
                                  "Sod1C", "HI_NLoci", "HI", "Mouse_ID")]
  
  mergedMiceTable = merge(mergedMiceTable, toadd, by = "Mouse_ID", all = T)
  
  mergedMiceTable = fillGapsAfterMerge(mergedMiceTable)
  
  # correct error tail length
  mergedMiceTable$Tail_length <- as.numeric(mergedMiceTable$Tail_length)
  
  # delete duplicated rows
  mergedMiceTable <- unique(mergedMiceTable)
  
  duplicated = mergedMiceTable$Mouse_ID[duplicated(mergedMiceTable$Mouse_ID)]
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  ############ Worms ############
  ## in WATWM dataset : Hymenolepis, Taenia, Rodentolepis, Mesocestoides,
  ## Calodium, Mastophorus, Trichuris, Heterakis, Aspiculuris+Syphacia
  
  # Hymenolepis
  mergedMiceTable$Hymenolepis <- rowSums(
    mergedMiceTable[c("Hymenolepis_microstoma", "Hymenolepis_diminiuta")], 
    na.rm = T)
  mergedMiceTable$Hymenolepis[with(mergedMiceTable, 
                                   is.na(mergedMiceTable["Hymenolepis_microstoma"]) &  
                                     is.na(mergedMiceTable["Hymenolepis_diminiuta"]))] <- NA
  
  # Taenia
  mergedMiceTable$Taenia <- rowSums(
    mergedMiceTable[c("Taenia_martis", "Taenia_taeniformis", 
                      "Catenotaenia_pusilla", "Cysticercus")], 
    na.rm = T)
  mergedMiceTable$Hymenolepis[with(mergedMiceTable, 
                                   is.na(mergedMiceTable["Taenia_martis"]) &  
                                     is.na(mergedMiceTable["Taenia_taeniformis"]) &
                                     is.na(mergedMiceTable["Catenotaenia_pusilla"]) &
                                     is.na(mergedMiceTable["Cysticercus"]))] <- NA
  
  # Aspiculuris_Syphacia
  mergedMiceTable$Aspiculuris_Syphacia <- rowSums(
    mergedMiceTable[c("Syphacia_obvelata", "Aspiculuris_tetraptera", "Mix_Syphacia_Aspiculuris")], 
    na.rm = T)
  mergedMiceTable$Aspiculuris_Syphacia[with(mergedMiceTable, 
                                            is.na(mergedMiceTable["Syphacia_obvelata"]) &  
                                              is.na(mergedMiceTable["Aspiculuris_tetraptera"]) &
                                              is.na(mergedMiceTable["Mix_Syphacia_Aspiculuris"]))] <- NA
  
  # Trichuris
  mergedMiceTable$Trichuris <- mergedMiceTable$Trichuris_muris
  
  # Heterakis
  mergedMiceTable$Heterakis <- mergedMiceTable$Heterakis_spumosa
  
  # Mastophorus
  mergedMiceTable$Mastophorus <- mergedMiceTable$Mastophorus_muris
  
  ## Dataframe to plot worms
  WormsDF <- mergedMiceTable[c("Mouse_ID", "Year", 
                               "Hymenolepis", "Taenia", "Aspiculuris_Syphacia", 
                               "Trichuris", "Heterakis", "Mastophorus")]
  WormsDF <- melt(WormsDF, id = c("Mouse_ID", "Year"))
  
  WormsDF$value <- as.numeric(as.character(WormsDF$value))
  
  ggplot(data=WormsDF, aes(x = variable, y=log10(value))) +
    geom_violin(aes(fill = variable))  +
    geom_jitter(size = 0.5, width = .2, alpha = .8) +
    theme_classic() +
    facet_wrap( ~ Year, nrow = 2) +
    theme(text = element_text(size = 15),
          axis.text = element_text(angle = 45, hjust = 1))+
    theme(legend.position="none")
  ## TODO 2014 and 2015 worms, not correct!!
  
  # Final cleaning, and save!
  mergedMiceTable$Longitude <- as.numeric(mergedMiceTable$Longitude)
  mergedMiceTable$Latitude <- as.numeric(mergedMiceTable$Latitude)
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  ## Remove useless mice:
  
  # wildpark Schorfheide (not needed, test)
  wsh <- c(paste0("AA_000", 1:9), paste0("AA_00", 10:46))
  # apodemus caught in 2016
  apd <- c("A_0001", "A_0002", "A_0003")
  # useless info
  useless <- c(wsh, apd)
  
  mergedMiceTable <- mergedMiceTable[!(mergedMiceTable$Mouse_ID %in% useless),]
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  # correct body length/weight
  mergedMiceTable$Body_length <- as.numeric(mergedMiceTable$Body_length)
  
  mergedMiceTable$Body_weight <- as.numeric(mergedMiceTable$Body_weight)
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  # Manual correction
  mergedMiceTable$Body_weight[!is.na(mergedMiceTable$Body_weight) &
                                mergedMiceTable$Body_weight > 100] <-
    mergedMiceTable$Body_weight[!is.na(mergedMiceTable$Body_weight) &
                                  mergedMiceTable$Body_weight > 100] / 1000
  
  mergedMiceTable$Body_length[which(mergedMiceTable$Body_length < 20)] <- 
    mergedMiceTable$Body_length[which(mergedMiceTable$Body_length < 20)] * 10
  
  # Body condition index as log body mass/log body length (Hayes et al. 2014)
  mergedMiceTable$BCI <- log(mergedMiceTable$Body_weight) / log(mergedMiceTable$Body_length)
  
  # Correct wrong HI (>1)
  mergedMiceTable$HI[mergedMiceTable$HI > 1 & !is.na(mergedMiceTable$HI)] <- 
    mergedMiceTable$HI[mergedMiceTable$HI > 1 & !is.na(mergedMiceTable$HI)]/1000
  
  # Correct wrong Long/Lat: 1. wrong multiplicative factor
  mergedMiceTable$Longitude[mergedMiceTable$Longitude > 100 &
                              !is.na(mergedMiceTable$Longitude)] <- 
    mergedMiceTable$Longitude[mergedMiceTable$Longitude > 100 &
                                !is.na(mergedMiceTable$Longitude)] / 1000
  
  mergedMiceTable$Latitude[mergedMiceTable$Latitude > 100 & 
                             !is.na(mergedMiceTable$Latitude)] <- 
    mergedMiceTable$Latitude[mergedMiceTable$Latitude > 100 & 
                               !is.na(mergedMiceTable$Latitude)] / 1000
  
  # Correct wrong Long/Lat: 2. wrong lat/lon inversion
  mergedMiceTable$Longitude.temp <- mergedMiceTable$Longitude
  
  mergedMiceTable$Longitude[!is.na(mergedMiceTable$Longitude) &
                              mergedMiceTable$Longitude >= 30] =
    mergedMiceTable$Latitude[!is.na(mergedMiceTable$Longitude) &
                               mergedMiceTable$Longitude >= 30]
  
  mergedMiceTable$Latitude[!is.na(mergedMiceTable$Latitude) &
                             mergedMiceTable$Latitude <= 20] =
    mergedMiceTable$Longitude.temp[!is.na(mergedMiceTable$Latitude) &
                                     mergedMiceTable$Latitude <= 20]
  
  mergedMiceTable = 
    mergedMiceTable[ , -which(names(mergedMiceTable) %in% c("Longitude.temp"))]
  
  # add farm (TODO better localisation)
  mergedMiceTable$farm <- paste0(mergedMiceTable$Longitude, mergedMiceTable$Latitude)
  
  # Cluster by localities: rounded to about 700 meters ???...
  # mergedMiceTable$Latitude <- round(mergedMiceTable$Latitude, 2)
  # mergedMiceTable$Longitude <- round(mergedMiceTable$Longitude, 2)
  
  ## remove empty rows
  mergedMiceTable <- mergedMiceTable[!is.na(mergedMiceTable$Mouse_ID),]
  
  ## remove duplicated rows
  mergedMiceTable <- unique(mergedMiceTable)
  
  ## Correct duplicated mice by hand
  duplicated = mergedMiceTable$Mouse_ID[duplicated(mergedMiceTable$Mouse_ID)]
  
  ## 26 June 2018, add Jarda new csv
  missingMice = mergedMiceTable$Mouse_ID[is.na(mergedMiceTable$HI)]
  
  newCsv = read.csv(paste0(pathToMyData, "Field_data/EmanuelData_26061018.csv"),
                    stringsAsFactors=F)
  newCsv$PIN = gsub("SK", "SK_", newCsv$PIN)
  
  ## Correct previous mistakes before merging for Latitude, Longitude, Transect, Sex, Code
  names(newCsv)[names(newCsv) %in% c("PIN", "X_Longit", "Y_Latit")] = 
    c("Mouse_ID", "Longitude", "Latitude")
  
  dataToAdd = newCsv[!newCsv$Mouse_ID %in% mergedMiceTable$Mouse_ID |
                       newCsv$Mouse_ID %in% missingMice,]
  
  mergedMiceTable =  merge(mergedMiceTable, dataToAdd, by = "Mouse_ID", all = T)
  
  mergedMiceTable = fillGapsAfterMerge(mergedMiceTable)
  
  # correct year manually
  mergedMiceTable$Year[
    mergedMiceTable$Mouse_ID %in% c("SK_2903", "SK_2904")] <- 2014
  mergedMiceTable$Year[
    mergedMiceTable$Mouse_ID %in% c("AA_0330", "AA_0450", "AA_0451", "AA_0452")] <- 2017
  
  # Add Eimeria information
  ## flotation
  flot <- read.csv(paste0(pathToMyData, "Eimeria_detection/FINALOocysts2015to2017.csv"))
  
  ## how many neubauer cells were counted 
  flot$Ncells <- apply(flot[paste0("N_oocysts_sq", 1:8)], 1, function(x) sum(!is.na(x)))
  
  flot$OPG <- rowSums(flot[,paste0("N_oocysts_sq", 1:8)], na.rm = T) / flot$Ncells * 10000 /
    (flot$PBS_dil_in_mL * flot$Feces_g)
  
  # flot$Mouse_ID[!flot$Mouse_ID %in% miceTable$Mouse_ID]
  # SK_3174 only missing :(
  mergedMiceTable <- merge(mergedMiceTable, flot, all = T)
  
  ## qPCr
  qpcr <- read.csv(paste0(pathToMyData, "/Eimeria_detection/FINALqpcrData_2016_2017_threshold5.csv"))
  #qpcr$Mouse_ID[!qpcr$Mouse_ID %in% miceTable$Mouse_ID]
  #all there :)
  mergedMiceTable <- merge(mergedMiceTable, 
                           qpcr[c("Mouse_ID", "delta_ct_ilwe_MminusE", "delta_ct_cewe_MminusE", "observer")],
                           all = T)
  
  ## species identification
  species <- read.csv(paste0(pathToMyData, "/Eimeria_detection/Eimeria_species_assignment_14_17.csv"))
  names(species)[names(species) %in% "Species"] <- "eimeriaSpecies"
  
  # space error damn!
  species$Mouse_ID <- gsub(" ", "", as.character(species$Mouse_ID))
  mergedMiceTable <- merge(mergedMiceTable, 
                           species[c("Mouse_ID", "n18S_Seq", "COI_Seq", "ORF470_Seq", "eimeriaSpecies")],
                           by = "Mouse_ID", all = T)
  
  # clean Year
  mergedMiceTable[grep("A_00", mergedMiceTable$Mouse_ID),"Year"] <- 2016
  
  ## TO BE CHANGED WITH TISSUE RESULTS
  mergedMiceTable$eimeriaSpecies <- as.character(mergedMiceTable$eimeriaSpecies)
  mergedMiceTable$eimeriaSpecies[mergedMiceTable$Mouse_ID %in% "AA_0111"] <- "Double_ferrisi_vermiformis"
  mergedMiceTable$eimeriaSpecies[mergedMiceTable$Mouse_ID %in% "AA_0244"] <- "Double_tbd"
  mergedMiceTable$eimeriaSpecies[mergedMiceTable$Mouse_ID %in% "AA_0245"] <- "Double_tbd"
  mergedMiceTable$eimeriaSpecies[mergedMiceTable$Mouse_ID %in% "AA_0436"] <- "Double_tbd"
  mergedMiceTable$eimeriaSpecies[mergedMiceTable$Mouse_ID %in% "AA_0497"] <- "Double_tbd"
  mergedMiceTable$eimeriaSpecies <- as.factor(mergedMiceTable$eimeriaSpecies)
  
  return(mergedMiceTable)
}

miceTable <- makeMiceTable()
# write.csv(miceTable, "data/MiceTable_fullEimeriaInfos_2014to2017.csv", row.names = F)

# Luke plays around from here :

# load new data

# 1. merge dissection table 2018 with miceTable

# 2. merge genotype table 2018 with new merged table (when Jarda send it)

## test and make sure it's correct (Luke + Alice talking about tests)

## Luke: Put working code (but not the tests) into the function

## Alice: Put tests into a testing script 




