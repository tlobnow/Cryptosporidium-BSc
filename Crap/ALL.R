library(dplyr)
library(tidyverse)
library(visdat)
library(ggplot2)
library(stringr)
library(data.table)
library(reshape)
library(ggmap)



## 1) Start with the new HI Data 
##    in this case we use the new HI Data we received via Email from Jarda (JardaNew, includes new data from 2018-19)
##    and the MiceTable from Alice (JardaTableOld, includes up to 2017 data)




#### MICE TABLE ALICE **********************************************************

#### BASE TABLE ################################################################
## get the data with genotype info and HI (2014-19)
JardaTableNew   <-  read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/EmanuelData.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
JardaTableOld   <-  read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/MiceTable_fullEimeriaInfos_2014to2017.csv")

## filter out Column called "X" and Years 2010 & 2011
JardaTableNew   <- JardaTableNew[!names(JardaTableNew) == "X"]
JardaTableNew   <- JardaTableNew[!JardaTableNew$Year %in% c(2010, 2011),]
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
vis_miss(JardaTable)
table(JardaTable$Transect)
JardaTable %>% count(is.na(HI))
JardaTable %>% count(Year)


# Homogenize Column Names
setnames(JardaTable,
         old = c("BW", "L", "LCd", "Ltestis", "Rtestis", "SemVes"),
         new = c("Body_weight", "Body_length", "Tail_length",
                 "Left_Testis", "Right_Testis", "Seminal_vesicles"),
         skip_absent = T)

         
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
                     "Litters", "NN", "MM", "FF", "Protocol", "Seminal_vesicles")

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
## every single info is available, so we will start with this and add data to it.
sel_JardaTable <- JardaTable[, colnames(JardaTable)%in%c(basics)]


# replace Year of AA_0330 from "2014" to "2017"
sel_JardaTable_AA_0330 <- sel_JardaTable %>% filter(Mouse_ID == "AA_0330") %>% mutate(Year = replace(Year, Year == 2014, 2017))
sel_JardaTable_not_AA_0330 <- sel_JardaTable %>% filter(Mouse_ID != "AA_0330")
sel_JardaTable <- full_join(sel_JardaTable_AA_0330, sel_JardaTable_not_AA_0330) %>% arrange(Mouse_ID)
rm(sel_JardaTable_AA_0330)
rm(sel_JardaTable_not_AA_0330)


#### Complement data with previous tables ######################################
## Dissection Data from 2014
Dis2014 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ14_Dissections_237-523.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
Dis2014$Mouse_ID      <- paste0(Dis2014$ID, "_", Dis2014$PIN)
Dis2014$Ectoparasites <- as.integer(Dis2014$Ectoparasites)

### homogenize column names
setnames(Dis2014, 
         old = c("BW", "L", "LCd", 
                 "ASP", "SYP", 
                 "TM","Taenia.taeniformis", "Locality", "SemVes"),
         new = c("Body_weight", "Body_length", "Tail_length", 
                 "Aspiculuris_tetraptera", "Syphacia_obvelata", 
                 "Trichuris_muris", "Taenia_taeniformis", "Address", "Seminal_vesicles"),
         skip_absent = T)

## check missing Data
    #vis_miss(Dis2014)
    # fill in the missing Transect Data
    Dis2014$Mouse_ID[is.na(Dis2014$Transect)]
    Dis2014_CZ <- Dis2014 %>%
      filter(is.na(Transect)) %>%
      mutate(Transect = "CZ")
    Dis2014_Transect_known <- Dis2014 %>%
      filter(!is.na(Transect))
    Dis2014 <- bind_rows(Dis2014_Transect_known, Dis2014_CZ) %>%
      arrange(Mouse_ID)
    rm(Dis2014_CZ)
    rm(Dis2014_Transect_known)

## join Basics table with 2014 dissection Data
sel_JardaTable1 <- full_join(sel_JardaTable, Dis2014) #vis_miss(sel_JardaTable1)

## Dissection Data from 2015
Dis2015 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ15_Mice_Parasite.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
Dis2015$Mouse_ID <- paste0(Dis2015$ID, "_", Dis2015$PIN)
    # remove transported mice
    Dis2015 <- Dis2015[!is.na(Dis2015$PIN),]
    # Homogenize Column Names
    setnames(Dis2015, 
         old = c("BW", "L", "LCd", 
                 "ASP", "SYP", 
                 "TM","Taenia.taeniformis", "SemVes"),
         new = c("Body_weight", "Body_length", "Tail_length", 
                 "Aspiculuris_tetraptera", "Syphacia_obvelata", 
                 "Trichuris_muris", "Taenia_taeniformis", "Seminal_vesicles"),
         skip_absent = T)
## separate Ectoparasites into 2 columns:
## Ectoparasites_logical = were Ectoparasites found
## Ectoparasites_num = Ectoparasites count
Dis2015$Ectoparasites_numeric <- Dis2015$Ectoparasites
Dis2015 <- Dis2015 %>%
  mutate(Ectoparasites_logical = ifelse(Ectoparasites_numeric == 0, FALSE,
                                        ifelse(Ectoparasites_numeric > 0, TRUE,
                                               ifelse(NA)))) %>%
  select(-c(Ectoparasites))

## check missing Data
    #vis_miss(Dis2015)
    # fill in the missing Transect Data
    Dis2015$Mouse_ID[is.na(Dis2015$Transect)]
    Dis2015_PL <- Dis2015 %>% filter(is.na(Transect)) %>% mutate(Transect = "Allo_PL")
    Dis2015_DE <- Dis2015 %>% filter(State == "Germany")
    Dis2015 <- bind_rows(Dis2015_PL, Dis2015_DE) %>% arrange(Mouse_ID)
    rm(Dis2015_DE)
    rm(Dis2015_PL)
    # fill in the missing Code Data
    Dis2015$Mouse_ID[is.na(Dis2015$Code)]
    Dis2015_Code_known <- Dis2015 %>% filter(!is.na(Code))
    Dis2015_Gola <- Dis2015 %>% filter(PIN %in% c(3469:3471)) %>% mutate(Code = "GOLA48")
    Dis2015_Recl <- Dis2015 %>% filter(PIN %in% c(3472,3473)) %>% mutate(Code = "RECL8")
    Dis2015_Trzes <- Dis2015 %>% filter(PIN %in% c(3474,3475)) %>% mutate(Code = "TRZES17")
    Dis2015 <- bind_rows(Dis2015_Code_known, Dis2015_Gola, Dis2015_Recl, Dis2015_Trzes) %>%
      arrange(Mouse_ID)
    rm(Dis2015_Code_known)
    rm(Dis2015_Gola)
    rm(Dis2015_Recl)
    rm(Dis2015_Trzes)

# join
sel_JardaTable2 <- full_join(sel_JardaTable1, Dis2015) #vis_miss(sel_JardaTable2)


## Dissection Data from 2016
Dis2016 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ16_Dissection_1-211.csv", na.strings=c(""," ","NA"), stringsAsFactors = F)[-c(1:2),]
Dis2016$Head.taken.       <- as.logical(Dis2016$Head.taken.)
Dis2016$Tail_length       <- as.double(Dis2016$Tail_length)

    # remove transported mice
    Dis2016 <- Dis2016[!is.na(Dis2016$Mouse_ID),]
    
    # separate values and remove old Testis column with wrong Testis input (contains Left_Testis / Right_Testis in 1 column)
    Dis2016_Testis <- Dis2016 %>% select(Mouse_ID, Testis) %>% filter(Testis != "NA") %>% filter(Testis != "0.181 both")
    Dis2016_Testis1 <- setDT(Dis2016_Testis)[, paste0("Testis", 1:2) := tstrsplit(Testis, "/")]
    # Homogenize Column Names
    setnames(Dis2016_Testis1,
             old = c("Testis1", "Testis2"),
             new = c("Left_Testis", "Right_Testis"))
    Dis2016 <- full_join(Dis2016, Dis2016_Testis1)
    # replace input of special case AA_0029: "0.181 both"
    Dis2016_AA_0029 <- Dis2016 %>%
      filter(Mouse_ID == "AA_0029") %>%
      mutate(Left_Testis = replace(Left_Testis, Left_Testis == "NA", 0.181),
             Right_Testis = replace(Right_Testis, Right_Testis == "NA", 0.181))
    Dis2016_not_AA_0029 <- Dis2016 %>% filter(Mouse_ID != "AA_0029")
    Dis2016 <- full_join(Dis2016_AA_0029, Dis2016_not_AA_0029)
    Dis2016 <- Dis2016 %>% select(-Testis)
    rm(Dis2016_Testis)
    rm(Dis2016_Testis1)
    rm(Dis2016_AA_0029)
    rm(Dis2016_not_AA_0029)
    Dis2016$Left_Testis <- as.double(Dis2016$Left_Testis)
    Dis2016$Right_Testis <- as.double(Dis2016$Right_Testis)
    # clean AA_0014 tail_length input
    Dis2016_AA_0014 <- Dis2016 %>% filter(Mouse_ID == "AA_0014") %>% mutate(Tail_length = replace(Tail_length, Tail_length == "NA", 55))
    Dis2016_not_AA_0014 <- Dis2016 %>% filter(Mouse_ID != "AA_0014")
    Dis2016 <- full_join(Dis2016_AA_0014, Dis2016_not_AA_0014) %>% arrange(Mouse_ID)
    rm(Dis2016_AA_0014)
    rm(Dis2016_not_AA_0014)
    
    # Homogenize Column Names
    setnames(Dis2016,
             old = c("Locality", "Notes", "Head.taken."),
             new = c("Address", "Note", "Head_taken"),
             skip_absent = T)
    
    
    Dis2016 <- Dis2016 %>%
      mutate(Ectoparasites_logical = ifelse(Ectoparasites == FALSE, FALSE,
                                            ifelse(Ectoparasites == "fleas", TRUE,
                                                   ifelse(NA))),
             Fleas = ifelse(Ectoparasites_logical == TRUE, TRUE,
                            ifelse(Ectoparasites_logical == FALSE, FALSE,
                                   ifelse(NA)))) %>%
      select(-Ectoparasites)
    
    
    
    # add worm Data
    Worm2016 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ16_Worms_47-211.csv", na.strings=c(""," ","NA"), stringsAsFactors = F)[-c(1:2),]
    setnames(Worm2016,
             old = c("SYP", "SyphaciaAspiculuris", "CP", "HM"),
             new = c("Syphacia_obvelata", "Aspiculuris_Syphacia", 
                     "Catenotaenia_pusilla", "Hymenolepis_microstoma"),
             skip_absent = T)
    
    Dis2016 <- full_join(Dis2016, Worm2016)

# join
sel_JardaTable3 <- full_join(sel_JardaTable2, Dis2016) #vis_miss(sel_JardaTable3)

## Dissection Data from 2017
Dis2017 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ17_Dissections_237-523.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
Dis2017$Spleen            <- as.double(Dis2017$Spleen)
Dis2017$Left_epididymis   <- as.double(Dis2017$Left_epididymis)
Dis2017$Embryo_left       <- as.double(Dis2017$Embryo_left)
Dis2017$Feces_weight      <- as.double(Dis2017$Feces_weight)
setnames(Dis2017,
         old = c("Left_testis_", "Right_testis", "Ectoparasites"),
         new = c("Left_Testis", "Right_Testis", "Ectoparasites_logical"))

Dis2017 <- Dis2017 %>% mutate(Ectoparasites_logical = ifelse(Ectoparasites_logical == FALSE, FALSE,
                                                             ifelse(Ectoparasites_logical == TRUE, TRUE,
                                                                    ifelse(Ectoparasites_logical == "TRUE (collected)", TRUE,
                                                                           ifelse(NA)))))




# join
sel_JardaTable4 <- full_join(sel_JardaTable3, Dis2017)


## Dissection Data from 2018
Dis2018 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ18_Dissections.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
Dis2018$Feces_weight <- as.double(Dis2018$Feces_weight)
setnames(Dis2018,
         old = c("Feces", "ASP", "SYP", 
                 "HET", "MART", "CP", 
                 "HD", "HM", 
                 "MM", "TM"),
         new = c("Feces_weight", "Aspiculuris_tetraptera", "Syphacia_obvelata",
                 "Heterakis_spumosa", "Taenia_martis", "Catenotaenia_pusilla",
                 "Hymenolepis_diminiuta", "Hymenolepis_microstoma", 
                 "Mastophorus_muris", "Trichuris_muris"),
         skip_absent = T)

# join
sel_JardaTable5 <- full_join(sel_JardaTable4, Dis2018)



## Dissection Data from 2019
Dis2019 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ19_Dissections.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
setnames(Dis2019,
         old = c("Feces", "ASP", "SYP", 
                 "HET", "MART", "CP", 
                 "HD", "HM", 
                 "MM", "TM"),
         new = c("Feces_weight", "Aspiculuris_tetraptera", "Syphacia_obvelata",
                 "Heterakis_spumosa", "Taenia_martis", "Catenotaenia_pusilla",
                 "Hymenolepis_diminiuta", "Hymenolepis_microstoma", 
                 "Mastophorus_muris", "Trichuris_muris"),
         skip_absent = T)

# join
sel_JardaTable6 <- full_join(sel_JardaTable5, Dis2019)
glimpse(sel_JardaTable6)
vis_miss(sel_JardaTable6)

## check for only NA rows and delete
which(!rowSums(!is.na(sel_JardaTable6)))




## Genotype Data from 2010-2019
Gen_2010_2019 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ10-19_Genotypes.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)


















## Add samples without HI in MiceTable but with one in oldJarda
missing <- MiceTable$Mouse_ID[is.na(MiceTable$HI)]

toadd <- JardaTableOld[!JardaTableOld$Mouse_ID %in% MiceTable$Mouse_ID,]


toadd <- rbind(toadd, JardaTableOld[JardaTableOld$Mouse_ID %in% missing,])

# MiceTable[MiceTable$Mouse_ID %in% missing 
#                 & is.na(MiceTable$Longitude) &
#             MiceTable$Mouse_ID %in% toadd$Mouse_ID,]

## setnames and merge
setnames(toadd,
         old = c("BW", "L", "LCd", 
                 "Aspiculuris", "Syphacia", 
                 "Trichuris","Taenia"),
         new = c("Body_weight", "Body_length", "Tail_length", 
                 "Aspiculuris_tetraptera", "Syphacia_obvelata", 
                 "Trichuris_muris", "Taenia_taeniformis"),
         skip_absent = T)

toadd = toadd[names(toadd)%in%c("REGion", "mtBamH", "Zfy2", "SRY1", "Y", "X332", "X347",
                                "X65", "Tsx", "Btk", "Syap1", "Es1", "Gpd1", "Idh1",
                                "Mpi", "Np", "Sod1", "Es1C", "Gpd1C", "Idh1C", "MpiC", "NpC",
                                "Sod1C", "HI_NLoci", "HI", "Mouse_ID")]

MiceTable <-  full_join(MiceTable, toadd)

MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>%
  arrange(Mouse_ID) %>%
  group_by(Mouse_ID) %>%
  fill(c(everything()), .direction = "downup") %>%
  ungroup()
MiceTable <- MiceTable %>% distinct(Mouse_ID, .keep_all = T) 
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
## select columns != all NA
MiceTable <- MiceTable %>% select(!which(!colSums(!is.na(MiceTable)))) #vis_miss(MiceTable, cluster = T, sort_miss = T)




#### WORMS  ###################################################################
## in WATWM dataset : Hymenolepis, Taenia, Rodentolepis, Mesocestoides,
## Calodium, Mastophorus, Trichuris, Heterakis, Aspiculuris+Syphacia


# Hymenolepis
#      MiceTable$Hymenolepis <- rowSums(MiceTable[c("Hymenolepis_microstoma", "Hymenolepis_diminiuta")], na.rm = T)
#      MiceTable$Hymenolepis[with(MiceTable,
#                                 is.na(MiceTable["Hymenolepis_microstoma"]) &  
#                                         is.na(MiceTable["Hymenolepis_diminiuta"]))] <- NA

# Taenia
#      MiceTable$Taenia <- rowSums(MiceTable[c("Taenia_martis", "Taenia_taeniformis", "Catenotaenia_pusilla", "Cysticercus")], na.rm = T)
#      MiceTable$Hymenolepis[with(MiceTable, 
#                                 is.na(MiceTable["Taenia_martis"]) &
#                                   is.na(MiceTable["Taenia_taeniformis"]) &
#                                   is.na(MiceTable["Catenotaenia_pusilla"]) &
#                                   is.na(MiceTable["Cysticercus"]))] <- NA

# Aspiculuris_Syphacia
#      MiceTable$Aspiculuris_Syphacia <- rowSums(MiceTable[c("Syphacia_obvelata", "Aspiculuris_tetraptera", "Mix_Syphacia_Aspiculuris")], na.rm = T)
#      MiceTable$Aspiculuris_Syphacia[with(MiceTable, 
#                                                is.na(MiceTable["Syphacia_obvelata"]) &  
#                                                  is.na(MiceTable["Aspiculuris_tetraptera"]) &
#                                                  is.na(MiceTable["Mix_Syphacia_Aspiculuris"]))] <- NA

# Trichuris
#      MiceTable$Trichuris <- MiceTable$Trichuris_muris

# Heterakis
#      MiceTable$Heterakis <- MiceTable$Heterakis_spumosa

# Mastophorus
#      MiceTable$Mastophorus <- MiceTable$Mastophorus_muris

## Dataframe to plot worms
#      WormsDF <- MiceTable[c("Mouse_ID", "Year", 
#                             "Hymenolepis", "Taenia", "Aspiculuris_Syphacia", 
#                             "Trichuris", "Heterakis", "Mastophorus")]
#      WormsDF <- melt(WormsDF, id = c("Mouse_ID", "Year"))
#      WormsDF$value <- as.numeric(as.character(WormsDF$value))
#      WormsDF %>%
#      ggplot(aes(variable, log10(value))) +
#        geom_violin(aes(fill = value))  +
#        geom_jitter(size = 0.5, width = .2, alpha = .8) +
#        theme_classic() +
#        facet_wrap( ~ Year, nrow = 2) +
#        theme(text = element_text(size = 15),
#              axis.text = element_text(angle = 45, hjust = 1))+
#        theme(legend.position="none")

## TODO 2014 and 2015 worms, not correct!! <- Alice
## still incorrect?? <- Tessa























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



library(dplyr)
library(tidyverse)
library(visdat)
library(ggplot2)
library(stringr)
library(data.table)
library(reshape)
library(ggmap)
library(lubridate)


#### Selected Columns we want to record for our samples ########################



#### 1) MAKE BASE TABLE ########################################################

## a) Load Data:  
JardaTableNew   <-  read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/EmanuelData.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
JardaTableNew   <- JardaTableNew[!names(JardaTableNew) == "X"]
JardaTableNew   <- JardaTableNew[!JardaTableNew$Year %in% c(2010, 2011),]

## b) Homogenize Column names: especially Mouse_ID should be universal!
setnames(JardaTableNew,
         old = c("PIN", "X_Longit", "Y_Latit"), 
         new = c("Mouse_ID", "Longitude", "Latitude"))

## homogenize Mouse_ID Row Names, old Mouse_IDs == SK_ and new Mouse_IDs == AA_
JardaTableNew$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_", x = JardaTableNew$Mouse_ID)

## remove Embryos (no interest for parasitic studies)
JardaTableNew   <- JardaTableNew[sapply(JardaTableNew$Mouse_ID, nchar) <= 7,]

## c) Check for missing data: --> with the old JardaTable from Alice
JardaTableOld   <-  read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/MiceTable_fullEimeriaInfos_2014to2017.csv")
JardaTableOld    <- JardaTableOld[!names(JardaTableOld) == "X"]
JardaTableOld    <- JardaTableOld[!JardaTableOld$Year %in% c(2010, 2011),]

## homogenize Mouse_ID Row Names, old Mouse_IDs == SK_ and new Mouse_IDs == AA_
JardaTableOld$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_", x = JardaTableOld$Mouse_ID)
JardaTableOld$HI_NLoci <- gsub(pattern = "HI ", replacement = "", x = JardaTableOld$HI_NLoci)
JardaTableOld$HI_NLoci <- as.integer(JardaTableOld$HI_NLoci)

## remove Embryos (no interest for parasitic studies)
JardaTableOld   <- JardaTableOld[sapply(JardaTableOld$Mouse_ID, nchar) <= 7,]



## d) Join: COMBINED TABLE (ALICE's TABLE + NEW DATA)
## old MiceTable (JardaTableOld)  +  new HI data (JardaTableNew)
JardaTable <- full_join(JardaTableNew, JardaTableOld)
rm(JardaTableNew)
rm(JardaTableOld)


## e) Fill in all missing data, delete Duplicates:
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]
JardaTable <- JardaTable %>%
  arrange(Mouse_ID) %>%
  group_by(Mouse_ID) %>%
  fill(c(everything()), .direction = "downup") %>%
  ungroup()
JardaTable <- JardaTable %>% distinct(Mouse_ID, .keep_all = T) 
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]
#vis_miss(JardaTable, sort_miss = T, cluster = T)

#### JARDATABLE ****

## a) Load Data
JardaTable <- JardaTable %>% select(!which(!colSums(!is.na(JardaTable))), -"Locality")

## b) Homogenize Column Names:
setnames(JardaTable,
         old = c("Body_weight", "Body_length", "Tail_length",
                 "BW", "L", "LCd", "Ltestis", "Rtestis", "SemVes", "Dissection"),
         new = c("Body_weight1", "Body_length1", "Tail_length1",
                 "Body_weight", "Body_length", "Tail_length",
                 "Left_Testis", "Right_Testis", "Seminal_vesicles", "Dissection_date"),
         skip_absent = T)
JardaTable$Embryo_left <- as.integer(JardaTable$Embryo_left)
#vis_miss(JardaTable, cluster = T, sort_miss = T)


#### Complement data with previous tables ######################################

#### 2014 DATA ****
## a) Load Data: Dissection Data from 2014 and columns that need type Conversions
Dis2014                     <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ14_Dissections_237-523.csv", na.strings=c(""," ","NA", "na"), stringsAsFactors = FALSE)
Dis2014$Mouse_ID            <- paste0(Dis2014$ID, "_", Dis2014$PIN)

## b) Homogenize column names:
setnames(Dis2014, 
         old = c("BW", "L", "LCd", 
                 "ASP", "SYP", "TM","Taenia.taeniformis", 
                 "SemVes", "Ectoparasites", "Left_testis", "Right_testis"),
         new = c("Body_weight", "Body_length", "Tail_length", 
                 "Aspiculuris_tetraptera", "Syphacia_obvelata", 
                 "Trichuris_muris", "Taenia_taeniformis", 
                 "Seminal_vesicles", "Ectoparasites_Count", 
                 "Left_Testis", "Right_Testis"),
         skip_absent = T)

## c) Check missing Data and fill in if possible:
# vis_miss(Dis2014)
## fill in the missing Transect Data
Dis2014$Mouse_ID[is.na(Dis2014$Transect)]
Dis2014_CZ <- Dis2014 %>% filter(is.na(Transect)) %>% mutate(Transect = "CZ")
Dis2014_Transect_known <- Dis2014 %>% filter(!is.na(Transect))
Dis2014 <- bind_rows(Dis2014_Transect_known, Dis2014_CZ) %>% arrange(Mouse_ID)
rm(Dis2014_CZ)
rm(Dis2014_Transect_known)


## d) Join Basics table with 2014 dissection Data
JardaTable <- full_join(JardaTable, Dis2014)


## e) Check Data for Duplicated Mouse_IDs, fill in the gaps, continue with non-NA columns
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]
JardaTable <- JardaTable %>%
  arrange(Mouse_ID) %>%
  group_by(Mouse_ID) %>%
  fill(c(everything()), .direction = "downup") %>%
  ungroup()
JardaTable <- JardaTable %>% distinct(Mouse_ID, .keep_all = T) 
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]    
## select columns != all NA
JardaTable <- JardaTable %>% select(!which(!colSums(!is.na(JardaTable))))
#vis_miss(JardaTable, sort_miss = T, cluster = T)



#### 2015 DATA ****
## a) Load Data: Dissection Data from 2015
Dis2015           <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ15_Mice_Parasite.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
Dis2015$Mouse_ID  <- paste0(Dis2015$ID, "_", Dis2015$PIN)
## remove transported mice
Dis2015 <- Dis2015[!is.na(Dis2015$PIN),]


## b) Homogenize Column Names
setnames(Dis2015, 
         old = c("BW", "L", "LCd", "ASP", "SYP", 
                 "TM","Taenia.taeniformis", "SemVes", "Left_testis", "Right_testis"),
         new = c("Body_weight", "Body_length", "Tail_length", 
                 "Aspiculuris_tetraptera", "Syphacia_obvelata", 
                 "Trichuris_muris", "Taenia_taeniformis", "Seminal_vesicles",
                 "Left_Testis", "Right_Testis"),
         skip_absent = T)

## Ectoparasites Column needs some adjustment (it was used as either logical or numeric column in the past)
## separate Ectoparasites into 2 columns:
## Ectoparasites_logical = were Ectoparasites found
## Ectoparasites_Count = how many Ectoparasites found
Dis2015$Ectoparasites_Count <- Dis2015$Ectoparasites
Dis2015 <- Dis2015 %>%
  mutate(Ectoparasites_Logical = ifelse(Ectoparasites_Count == 0, FALSE,
                                        ifelse(Ectoparasites_Count > 0, TRUE,
                                               ifelse(NA)))) %>% select(-c(Ectoparasites))


## d) Join:
JardaTable <- full_join(Dis2015, JardaTable)  

## e) Check Data for Duplicated Mouse_IDs, fill in the gaps, continue with non-NA columns
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]
JardaTable <- JardaTable %>%
  arrange(Mouse_ID) %>%
  group_by(Mouse_ID) %>%
  fill(c(everything()), .direction = "downup") %>%
  ungroup()
JardaTable <- JardaTable %>% distinct(Mouse_ID, .keep_all = T) 
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]
## select columns != all NA
JardaTable <- JardaTable %>% select(!which(!colSums(!is.na(JardaTable))))
#vis_miss(JardaTable, sort_miss = T, cluster = T)


#### 2016 DATA ****
## a) Load Data: Dissection Data from 2016
Dis2016                   <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ16_Dissection_1-211.csv", na.strings=c(""," ","NA"), stringsAsFactors = F)[-c(1:2),]
Dis2016$Head.taken.       <- as.logical(Dis2016$Head.taken.)
Dis2016$Tail_length       <- as.double(Dis2016$Tail_length)
Dis2016$Left_epididymis   <- as.double(Dis2016$Left_epididymis)
JardaTable$Embryo_left   <- as.integer(JardaTable$Embryo_left)
JardaTable$Tail_length   <- as.double(JardaTable$Tail_length)



Dis2016 <- Dis2016[!is.na(Dis2016$Mouse_ID),]
## add separate Parasite Data
Worm2016 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ16_Worms_47-211.csv", na.strings=c(""," ","NA"), stringsAsFactors = F)[-c(1:2),]



## b) Homogenize Column Names:
setnames(Dis2016,
         old = c("Notes", "Head.taken."),
         new = c("Note", "Head_taken"),
         skip_absent = T)
setnames(Worm2016,
         old = c("SYP", "SyphaciaAspiculuris", "CP", "HM", "Aspiculuris",
                 "Trichuris", "Heterakis", "Mastophorus"),
         new = c("Syphacia_obvelata", "Aspiculuris_Syphacia", 
                 "Catenotaenia_pusilla", "Hymenolepis_microstoma", "Aspiculuris_tetraptera",
                 "Trichuris_muris", "Heterakis_spumosa", "Mastophorus_muris"),
         skip_absent = T)
## Ectoparasites Column needs some adjustment
## Instead of logical, "FALSE" and "fleas" was supplied..
## Therefore continuation of "Ectoparasites_Logical" 
## and introduction of "Fleas" column
Dis2016 <- Dis2016 %>% 
  mutate(Ectoparasites_Logical = ifelse(Ectoparasites == FALSE, FALSE, ifelse(Ectoparasites == "fleas", TRUE, ifelse(NA))),
         Fleas = ifelse(Ectoparasites_Logical == TRUE, TRUE, ifelse(Ectoparasites_Logical == FALSE, FALSE, ifelse(NA)))) %>% select(-Ectoparasites)


## d) Join:
## Worm + Dissection data
Dis_Worm2016         <- full_join(Worm2016, Dis2016)
## Join Dis_Worm2016 with Jarda
JardaTable <- full_join(Dis_Worm2016, JardaTable)


## e) Check Data for Duplicated Mouse_IDs, fill in the gaps, continue with non-NA columns
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]
JardaTable <- JardaTable %>%
  arrange(Mouse_ID) %>%
  group_by(Mouse_ID) %>%
  fill(c(everything()), .direction = "downup") %>%
  ungroup()
JardaTable <- JardaTable %>% distinct(Mouse_ID, .keep_all = T) 
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]
## select columns != all NA
JardaTable <- JardaTable %>% select(!which(!colSums(!is.na(JardaTable))))
#vis_miss(JardaTable, sort_miss = T, cluster = T)


#### 2017 DATA ****
## a) Load Data: Dissection Data from 2017
Dis2017 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ17_Dissections_237-523.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
Dis2017$Spleen            <- as.double(Dis2017$Spleen)
Dis2017$Left_epididymis   <- as.double(Dis2017$Left_epididymis)
Dis2017$Embryo_left       <- as.double(Dis2017$Embryo_left)
Dis2017$Tail_length       <- as.double(Dis2017$Tail_length)
Dis2017$Feces_weight      <- as.numeric(as.character(Dis2017$Feces_weight))
Dis2017$Feces_weight[Dis2017$Feces_weight > 100 & !is.na(Dis2017$Feces_weight)] <- Dis2017$Feces_weight[Dis2017$Feces_weight > 100 & !is.na(Dis2017$Feces_weight)] / 1000

JardaTable$Feces_weight   <- as.numeric(as.character(JardaTable$Feces_weight))
JardaTable$Feces_weight[JardaTable$Feces_weight > 100 & !is.na(JardaTable$Feces_weight)] <- JardaTable$Feces_weight[JardaTable$Feces_weight > 100 & !is.na(JardaTable$Feces_weight)] / 1000

## b) Homogenize Column Names:
setnames(Dis2017,
         old = c("Left_testis_", "Right_testis", "Ectoparasites"),
         new = c("Left_Testis", "Right_Testis", "Ectoparasites_logical"))
Dis2017 <- Dis2017 %>% mutate(Ectoparasites_Logical = ifelse(Ectoparasites_logical == FALSE, FALSE,
                                                             ifelse(Ectoparasites_logical == TRUE, TRUE,
                                                                    ifelse(Ectoparasites_logical == "TRUE (collected)", TRUE,
                                                                           ifelse(NA))))) %>% select(-Ectoparasites_logical)



## d) Join:
JardaTable <- full_join(JardaTable, Dis2017)

## e) Check Data for Duplicated Mouse_IDs, fill in the gaps, continue with non-NA columns
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]
JardaTable <- JardaTable %>%
  arrange(Mouse_ID) %>%
  group_by(Mouse_ID) %>%
  fill(c(everything()), .direction = "downup") %>%
  ungroup()
JardaTable <- JardaTable %>% distinct(Mouse_ID, .keep_all = T) 
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]
## select columns != all NA
JardaTable <- JardaTable %>% select(!which(!colSums(!is.na(JardaTable))))
#vis_miss(JardaTable, sort_miss = T, cluster = T)

#### 2018 DATA ****
## a) Load Data: Dissection Data from 2018
Dis2018 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ18_Dissections.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
Dis2018$Feces <- as.double(Dis2018$Feces)

## b) Homogenize Column Names:
setnames(Dis2018,
         old = c("Feces", "ASP", "SYP", "HET", "MART", "CP", "HD", "HM", "MM", "TM",
                 "Ectoparasites"),
         new = c("Feces_weight", "Aspiculuris_tetraptera", "Syphacia_obvelata",
                 "Heterakis_spumosa", "Taenia_martis", "Catenotaenia_pusilla",
                 "Hymenolepis_diminiuta", "Hymenolepis_microstoma", 
                 "Mastophorus_muris", "Trichuris_muris", "Ectoparasites_Logical"),
         skip_absent = T)

## d) Join:
JardaTable <- full_join(Dis2018, JardaTable)

## e) Check Data for Duplicated Mouse_IDs, fill in the gaps, continue with non-NA columns
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]
JardaTable <- JardaTable %>%
  arrange(Mouse_ID) %>%
  group_by(Mouse_ID) %>%
  fill(c(everything()), .direction = "downup") %>%
  ungroup()
JardaTable <- JardaTable %>% distinct(Mouse_ID, .keep_all = T) 
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]
## select columns != all NA
JardaTable <- JardaTable %>% select(!which(!colSums(!is.na(JardaTable))))
#vis_miss(JardaTable, sort_miss = T, cluster = T)


#### 2019 DATA ****
## a) Load Data: Dissection Data from 2019
Dis2019 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ19_Dissections.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)

## b) Homogenize Column Names:
setnames(Dis2019,
         old = c("Feces", "ASP", "SYP", 
                 "HET", "MART", "CP", 
                 "HD", "HM", 
                 "MM", "TM", "Ectoparasites"),
         new = c("Feces_weight", "Aspiculuris_tetraptera", "Syphacia_obvelata",
                 "Heterakis_spumosa", "Taenia_martis", "Catenotaenia_pusilla",
                 "Hymenolepis_diminiuta", "Hymenolepis_microstoma", 
                 "Mastophorus_muris", "Trichuris_muris", "Ectoparasites_Logical"),
         skip_absent = T)

## d) Join:
JardaTable <- full_join(Dis2019, JardaTable)

## e) Check Data for Duplicated Mouse_IDs, fill in the gaps, continue with non-NA columns
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]
JardaTable <- JardaTable %>%
  arrange(Mouse_ID) %>%
  group_by(Mouse_ID) %>%
  fill(c(everything()), .direction = "downup") %>%
  ungroup()
JardaTable <- JardaTable %>% distinct(Mouse_ID, .keep_all = T) 
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]
## select columns != all NA
JardaTable <- JardaTable %>% select(!which(!colSums(!is.na(JardaTable))))
#vis_miss(JardaTable, sort_miss = T, cluster = T)


#### JOIN BASICS AND DISSECTION DATA WITH GENOTYPE INFO ########################

## a) Load Data: Genotype Data from 2010-2019
Gen_2010_2019   <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ10-19_Genotypes.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
Gen_2010_2019   <- Gen_2010_2019[!names(Gen_2010_2019) == "X"]
Gen_2010_2019   <- Gen_2010_2019 %>% select(- c(CH4.141I14, CH17.106))
Gen_2010_2019$EH_ID <- gsub(pattern = "SK", replacement = "SK_", x = Gen_2010_2019$EH_ID)


## b) Homogenize Column Names:
setnames(Gen_2010_2019,
         old = c("EH_ID"),
         new = c("Mouse_ID"),
         skip_absent = T)

## c) Check missing Data: 


## d) Join:
MiceTable <- full_join(Gen_2010_2019, JardaTable)

## e) Check Data for Duplicated Mouse_IDs, fill in the gaps, continue with non-NA columns
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>%
  arrange(Mouse_ID) %>%
  group_by(Mouse_ID) %>%
  fill(c(everything()), .direction = "downup") %>%
  ungroup()
MiceTable <- MiceTable %>% distinct(Mouse_ID, .keep_all = T) 
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
## select columns != all NA
MiceTable <- MiceTable %>% select(!which(!colSums(!is.na(MiceTable))))
#vis_miss(MiceTable, cluster = T)

## WRITE CSV write.csv(MiceTable, "MiceTable_Dirty.csv")


################################################################################
######################### MICE TABLE CLEANING STEPS ############################
## Date_count fix by Sep?
## improve Date uniformity of Trap_Date 


#### What Columns to keep?
basics <- c("Mouse_ID", "Transect", "Code", "Region", "Sex", "Longitude", "Latitude", "Year", "State")


gen.loci <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
              "Es1", "Gpd1", "Idh1", "Mpi", "Np", "Sod1", "Es1C", "Gpd1C",
              "Idh1C", "MpiC", "NpC", "Sod1C", "HI_NLoci",
              "HI", "Zfy2", "SRY1", "Y")

dissection.cols <- c("Body_Weight", "Body_Length", "Tail_Length", "Spleen", 
                     "Left_Testis", "Right_Testis", "Seminal_Vesicles_Weight", 
                     "Sperm", "Left_Epididymis", "Right_Epididymis", 
                     "Arrival", "Wean", "Death", "Dissection_date", "DaysInLab",
                     "Lepid", "Uterus", "Ovaria", "Grav", "Litters", "NN", 
                     "MM", "FF", "Protocol")

parasite.cols <- c("Aspiculuris_tetraptera", "Syphacia_obvelata", "Aspiculuris_Syphacia",
                   "Trichuris", "Taenia", "Flea", "Mix_Syphacia_Aspiculuris", "Heligmosomoides_polygurus",
                   "Heterakis","Heterakis", "Mastophorus","Hymenolepis",
                   "Ectoparasites", "Ectoparasites_Count", "Ectoparasites_Logical",
                   "Worms_presence")

oocyst.cols <- c("counter", "Feces_Weight", "Date_count", "N_oocysts_sq1",
                 "N_oocysts_sq2", "N_oocysts_sq3",  "N_oocysts_sq4",
                 "N_oocysts_sq5", "N_oocysts_sq6", "N_oocysts_sq7",
                 "N_oocysts_sq8", "mean_neubauer", "PBS_dil_in_mL", 
                 "OPG", "Ncells")

EqPCR.cols <- c("delta_ct_ilwe_MminusE", "delta_ct_cewe_MminusE",
                ## from 2018 on AN IMPORTANT IMPROVEMENT!!!
                "MC.Eimeria")

EimGeno.cols <- c("n18S_Seq", "COI_Seq", "ORF470_Seq", "eimeriaSpecies")


#### COLUMN SEPARATION #########################################################
MiceTable   <- read.csv("MiceTable_Dirty.csv")
MiceTable   <- MiceTable[!names(MiceTable) == "X.1"]

## Testis Separation
## Wrong data input for "Testis": instead of individual Left_Testis or 
## Right_Testis data, a combination of "Left_Testis/Right_Testis" was supplied
## Separation of that Column into 2 Columns, original "Testis"  col. discarded
MiceTable_Sep <- MiceTable %>% select(Mouse_ID, Testis) %>% filter(Testis != "NA")
MiceTable_Sep <- setDT(MiceTable_Sep)[, paste0("Testis", 1:2) := tstrsplit(Testis, "/")]
setnames(MiceTable_Sep, old = c("Testis1", "Testis2"), new = c("Left_Testis", "Right_Testis"), skip_absent = T)
MiceTable_Sep$Left_Testis <- as.double(MiceTable_Sep$Left_Testis)
MiceTable_Sep$Right_Testis <- as.double(MiceTable_Sep$Right_Testis)
## join
MiceTable <- full_join(MiceTable, MiceTable_Sep) %>% select(-Testis)
## remove Duplicates
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Testis_mass Separation
## Wrong data input for "Testis_mass": instead of individual Left_Testis or 
## Right_Testis data, a combination of "Left_Testis/Right_Testis" was supplied
## Separation of that Column into 2 Columns, original "Testis"  col. discarded
MiceTable_Sep <- MiceTable %>% select(Mouse_ID, Testis_mass) %>% filter(Testis_mass != "NA")
MiceTable_Sep <- setDT(MiceTable_Sep)[, paste0("Testis_mass", 1:2) := tstrsplit(Testis_mass, "/")]
setnames(MiceTable_Sep, old = c("Testis_mass1", "Testis_mass2"), new = c("Left_Testis1", "Right_Testis1"), skip_absent = T)
MiceTable_Sep$Left_Testis1 <- as.double(MiceTable_Sep$Left_Testis1)
MiceTable_Sep$Right_Testis1 <- as.double(MiceTable_Sep$Right_Testis1)
## join
MiceTable <- full_join(MiceTable, MiceTable_Sep) %>% select(-Testis_mass)
rm(MiceTable_Sep)
## remove Duplicates
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Sex Separation, contains Status info of some Samples
## aka young, 
MiceTable_Sex <- MiceTable %>% select(Mouse_ID, Sex) %>% filter(Sex != "NA")
MiceTable_Sex <- setDT(MiceTable_Sex)[, paste0("Sex", 1:2) := tstrsplit(Sex, " ")]
setnames(MiceTable_Sex, old = c("Sex1", "Sex2"), new = c("SEX", "Status2"), skip_absent = T)
## join
MiceTable <- full_join(MiceTable, MiceTable_Sex) %>% select(-Sex)
rm(MiceTable_Sex)
## remove Duplicates
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


#### COLUMN CORRECTION #########################################################

## similar Columns
## Aspiculuris_Syphacia == "Syphacia_obvelata", "Aspiculuris_tetraptera", "Mix_Syphacia_Aspiculuris"
MT_Aspiculuris <- MiceTable %>% select(Mouse_ID, Aspiculuris_Syphacia, "Syphacia_obvelata", "Aspiculuris_tetraptera", "Mix_Syphacia_Aspiculuris")

## Body_Weight
MT_Body_Weight <- MiceTable %>% select(Mouse_ID, Body_weight, Body_weight1)
MT_Body_Weight <- MT_Body_Weight %>% pivot_longer(names_to = "Temp", values_to = "Body_Weight", cols = c(Body_weight, Body_weight1)) %>% 
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Body_Weight) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Body_Weight$Mouse_ID[duplicated(MT_Body_Weight$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Body_Weight) %>% select(-c(Body_weight, Body_weight1))
rm(MT_Body_Weight)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 



## Body_Length == "Body_length", "Body_length1" -----> correct one super low value 8.9 or something.. --> 8.9
MT_Body_Length <- MiceTable %>% select(Mouse_ID, Body_length, Body_length1)
MT_Body_Length <- MT_Body_Length %>% pivot_longer(names_to = "Temp", values_to = "Body_Length", cols = c(Body_length, Body_length1)) %>% 
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Body_Length) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Body_Length$Mouse_ID[duplicated(MT_Body_Length$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Body_Length) %>% select(-c(Body_length, Body_length1))
rm(MT_Body_Length)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Ectoparasites == "Ectoparasites", "Ectoparasites_Logical"
MiceTable$Ectoparasites <- as.logical(MiceTable$Ectoparasites) 
MT_Ectoparasites <- MiceTable %>% select(Mouse_ID, Ectoparasites, Ectoparasites_Logical) %>% pivot_longer(names_to = "Temp", values_to = "Ectoparasites_Logical", cols = c(Ectoparasites, Ectoparasites_Logical)) %>%
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Ectoparasites_Logical) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Ectoparasites$Mouse_ID[duplicated(MT_Ectoparasites$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Ectoparasites) %>% select(-c(Ectoparasites))
rm(MT_Ectoparasites)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## Left_Epididymis
## Left_Epididymis == "Left_epididymis", "Left.epididymis.weight"
MiceTable$Left.epididymis.weight <- as.double(MiceTable$Left.epididymis.weight)

MT_Left_Epididymis <- MiceTable %>% select(Mouse_ID, Left_epididymis, Left.epididymis.weight) %>% pivot_longer(names_to = "Temp", values_to = "Left_Epididymis", cols = c(Left_epididymis, Left.epididymis.weight)) %>% 
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Left_Epididymis) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Left_Epididymis$Mouse_ID[duplicated(MT_Left_Epididymis$Mouse_ID)]  
## join
MiceTable <- full_join(MiceTable, MT_Left_Epididymis) %>% select(-c(Left_epididymis, Left.epididymis.weight))
rm(MT_Left_Epididymis)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Feces_weight == "Feces_weight", "Feces_g"
MT_Feces_Weight <- MiceTable %>% select(Mouse_ID, Feces_weight, Feces_g) %>% pivot_longer(names_to = "Temp",  values_to = "Feces_Weight", cols = c(Feces_weight, Feces_g)) %>%
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Feces_Weight) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Feces_Weight$Mouse_ID[duplicated(MT_Feces_Weight$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Feces_Weight) %>% select(-c(Feces_weight, Feces_g))
rm(MT_Feces_Weight)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## Fleas == "Flea", "Fleas"
MiceTable$Flea <- as.logical(MiceTable$Flea)
MT_Fleas <- MiceTable %>% select(Mouse_ID, Flea, Fleas) %>% pivot_longer(names_to = "Temp", values_to = "Fleas", cols = c(Flea, Fleas)) %>%
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Fleas)
MT_Fleas <- MT_Fleas %>% distinct(Mouse_ID, .keep_all = T) 
MT_Fleas$Mouse_ID[duplicated(MT_Fleas$Mouse_ID)]
## join
MiceTable <- full_join(MiceTable, MT_Fleas) %>% select(-c(Flea))
rm(MT_Fleas)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Head_taken == "Head_taken", "Head.taken."
MT_Head_Taken <- MiceTable %>% select(Mouse_ID, Head_taken, Head.taken.) %>%  pivot_longer(names_to = "Temp", values_to = "Head_Taken", cols = c(Head_taken, Head.taken.)) %>% 
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Head_Taken)
MT_Head_Taken <- MT_Head_Taken %>% distinct(Mouse_ID, .keep_all = T) 
MT_Head_Taken$Mouse_ID[duplicated(MT_Head_Taken$Mouse_ID)]   
## join
MiceTable <- full_join(MiceTable, MT_Head_Taken) %>% select(-c(Head.taken., Head_taken))
rm(MT_Head_Taken)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## Heterakis == "Heterakis", "Heterakis_spumosa"
MT_Heterakis <- MiceTable %>% select(Mouse_ID, Heterakis, Heterakis_spumosa) %>%pivot_longer(names_to = "Temp", values_to = "Heterakis", cols = c(Heterakis, Heterakis_spumosa)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Heterakis)
MT_Heterakis <- MT_Heterakis %>% distinct(Mouse_ID, .keep_all = T) 
MT_Heterakis$Mouse_ID[duplicated(MT_Heterakis$Mouse_ID)]  
## join
MiceTable <- full_join(MiceTable, MT_Heterakis) %>% select(-c(Heterakis_spumosa))
rm(MT_Heterakis)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## Hymenolepis == "Hymenolepis", "Hymenolepis_microstoma", "Hymenolepis_diminiuta"
MT_Hymenolepis <- MiceTable %>% select(Mouse_ID, Hymenolepis, Hymenolepis_diminiuta, Hymenolepis_microstoma)
MT_Hymenolepis <- MT_Hymenolepis %>% pivot_longer(names_to = "Temp",  values_to = "Hymenolepis", cols = c(Hymenolepis, Hymenolepis_diminiuta, Hymenolepis_microstoma)) %>%
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Hymenolepis)
MT_Hymenolepis <- MT_Hymenolepis %>% distinct(Mouse_ID, .keep_all = T) 
MT_Hymenolepis$Mouse_ID[duplicated(MT_Hymenolepis$Mouse_ID)]   
## join
MiceTable <- full_join(MiceTable, MT_Hymenolepis)
rm(MT_Hymenolepis)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## Latitude == "Latitude", "longitude" ***** MIX UP WITH LONG *****
MT_Latitude <- MiceTable %>% select(Mouse_ID, Latitude, longitude) %>% pivot_longer(names_to = "Temp",  values_to = "Latitude", cols = c(Latitude, longitude)) %>%
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%   fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Latitude)
MT_Latitude <- MT_Latitude %>% distinct(Mouse_ID, .keep_all = T) 
MT_Latitude$Mouse_ID[duplicated(MT_Latitude$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Latitude) %>% select(-c(longitude))
rm(MT_Latitude)

## Longitude == "Longitude", "latitude" ***** MIX UP WITH LAT *****
MT_Longitude <- MiceTable %>% select(Mouse_ID, Longitude, latitude)  %>% pivot_longer(names_to = "Temp", values_to = "Longitude", cols = c(Longitude, latitude)) %>% 
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Longitude)
MT_Longitude <- MT_Longitude %>% distinct(Mouse_ID, .keep_all = T) 
MT_Longitude$Mouse_ID[duplicated(MT_Longitude$Mouse_ID)]   
## join
MiceTable <- full_join(MiceTable, MT_Longitude) %>% select(-c(latitude))
rm(MT_Longitude)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## Liver == "Liver", "Liver_mass"
MT_Liver <- MiceTable %>% select(Mouse_ID, Liver_mass, Liver) %>% pivot_longer(names_to = "Temp",  values_to = "Liver", cols = c(Liver, Liver_mass)) %>%  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Liver)
MT_Liver <- MT_Liver %>% distinct(Mouse_ID, .keep_all = T) 
MT_Liver$Mouse_ID[duplicated(MT_Liver$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Liver) %>% select(-c(Liver_mass))
rm(MT_Liver)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## Mastophorus == "Mastophorus", "Mastophorus_muris", "Mastaphorus"
MT_Mastophorus <- MiceTable %>% select(Mouse_ID, Mastaphorus, Mastophorus, Mastophorus_muris) %>% pivot_longer(names_to = "Temp",  values_to = "Mastophorus", cols = c(Mastaphorus, Mastophorus, Mastophorus_muris)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>% select(Mouse_ID, Mastophorus)
MT_Mastophorus <- MT_Mastophorus %>% distinct(Mouse_ID, .keep_all = T) 
MT_Mastophorus$Mouse_ID[duplicated(MT_Mastophorus$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Mastophorus) %>% select(-c(Mastaphorus, Mastophorus_muris))
rm(MT_Mastophorus)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## Notes == "comments", "Note", "Notes", ("Embryo")
MT_Notes <- MiceTable %>% select(Mouse_ID, Note, Notes, comments)  %>% pivot_longer(names_to = "Temp",  values_to = "Notes", cols = c(Note, Notes, comments)) %>%  
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Notes)
MT_Notes <- MT_Notes %>% distinct(Mouse_ID, .keep_all = T) 
MT_Notes$Mouse_ID[duplicated(MT_Notes$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Notes) %>% select(-c(Note, comments))
rm(MT_Notes)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 
#vis_miss(MiceTable, cluster = T)

## Ovaria ==
## Right_Ovarium_Weight == "Right_ovarium", Right.ovarium.weight"
MT_Right_Ovarium_Weight <- MiceTable %>% select(Mouse_ID, Right_ovarium, Right.ovarium.weight) %>% pivot_longer(names_to = "Temp",  values_to = "Right_Ovarium_Weight", cols = c(Right_ovarium, Right.ovarium.weight)) %>%  
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Right_Ovarium_Weight) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Right_Ovarium_Weight$Mouse_ID[duplicated(MT_Right_Ovarium_Weight$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Right_Ovarium_Weight) %>% select(-c(Right_ovarium, Right.ovarium.weight))
rm(MT_Right_Ovarium_Weight)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## Left_Ovarium_Weight == "Left_ovarium", "Left.ovarium.weight"
MT_Left_Ovarium_Weight <- MiceTable %>% select(Mouse_ID, Left_ovarium, Left.ovarium.weight) %>% pivot_longer(names_to = "Temp",  values_to = "Left_Ovarium_Weight", cols = c(Left_ovarium, Left.ovarium.weight)) %>%  
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Left_Ovarium_Weight) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Left_Ovarium_Weight$Mouse_ID[duplicated(MT_Left_Ovarium_Weight$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Left_Ovarium_Weight) %>% select(-c(Left_ovarium, Left.ovarium.weight))
rm(MT_Left_Ovarium_Weight)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 



## Region == "Region", "REGion"
MT_Region <- MiceTable %>% select(Mouse_ID, Region, REGion) %>% pivot_longer(names_to = "Temp",  values_to = "Region", cols = c(Region, REGion)) %>%
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Region)
MT_Region <- MT_Region %>% distinct(Mouse_ID, .keep_all = T) 
MT_Region$Mouse_ID[duplicated(MT_Region$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Region) %>% select(-c(REGion))
rm(MT_Region)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## Sex
#        MT_Sex <- MiceTable %>% select(Mouse_ID, Sex, SEX)  %>% pivot_longer(names_to = "Temp",  values_to = "Sex", cols = c(Sex, SEX)) %>% 
#          arrange(Mouse_ID) %>%   group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Sex) %>% distinct(Mouse_ID, .keep_all = T) 
#        MT_Sex$Mouse_ID[duplicated(MT_Sex$Mouse_ID)]    
#        ## join
#        MiceTable <- full_join(MiceTable, MT_Sex) %>% select(-c(Sex2))
#        rm(MT_Sex)
#        ## eliminate duplicates in MiceTable
#        MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
#        MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## Status
MT_Status <- MiceTable %>% select(Mouse_ID, Status, Status2)  %>% pivot_longer(names_to = "Temp",  values_to = "Status", cols = c(Status, Status2)) %>% 
  arrange(Mouse_ID) %>%   group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Status) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Status$Mouse_ID[duplicated(MT_Status$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Status) %>% select(-c(Status2))
rm(MT_Status)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 



## Seminal_Vesicles == "SemVes", "Seminal_vesicles", "Seminal_vesicle_weight"
MT_Seminal_Vesicles_Weight <- MiceTable %>% select(Mouse_ID, Seminal_vesicles, Seminal.vesicle.weight)  %>% pivot_longer(names_to = "Temp",   values_to = "Seminal_Vesicles_Weight", cols = c(Seminal_vesicles, Seminal.vesicle.weight)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Seminal_Vesicles_Weight) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Seminal_Vesicles_Weight$Mouse_ID[duplicated(MT_Seminal_Vesicles_Weight$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Seminal_Vesicles_Weight) %>% select(-c(Seminal_vesicles, Seminal.vesicle.weight))
rm(MT_Seminal_Vesicles_Weight)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 



## Spleen == "Spleen", "Spleen_mass"
MiceTable$Spleen_mass <- as.double(MiceTable$Spleen_mass)
MT_Spleen <- MiceTable %>% select(Mouse_ID, Spleen, Spleen_mass)  %>% pivot_longer(names_to = "Temp",  values_to = "Spleen",cols = c(Spleen, Spleen_mass)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Spleen) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Spleen$Mouse_ID[duplicated(MT_Spleen$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Spleen) %>% select(-c(Spleen_mass))
rm(MT_Spleen)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Taenia == "Taenia_martis", "Taenia_taeniformis", "Catenotaenia_pusilla", "Cysticercus"
MT_Taenia <- MiceTable %>% select(Mouse_ID, Taenia_martis, Taenia_taeniformis, Catenotaenia_pusilla, Cysticercus, Taenia)
MT_Taenia <- MT_Taenia %>%
  pivot_longer(names_to = "Temp", 
               values_to = "Taenia",
               cols = c(Taenia_martis, Taenia_taeniformis, Catenotaenia_pusilla, Cysticercus, Taenia))
MT_Taenia <- MT_Taenia %>% 
  arrange(Mouse_ID) %>% 
  group_by(Mouse_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  ungroup() %>% 
  select(Mouse_ID, Taenia)
MT_Taenia <- MT_Taenia %>% distinct(Mouse_ID, .keep_all = T) 
MT_Taenia$Mouse_ID[duplicated(MT_Taenia$Mouse_ID)]    
MiceTable <- full_join(MiceTable, MT_Taenia) %>% select(-c(Taenia_martis, Taenia_taeniformis, Catenotaenia_pusilla, Cysticercus))
rm(MT_Taenia)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Testes == "Testes"
## Left_Testis  == "Left_Testis", "Left_Testis1", ("Left_testis"), "Left_Testis_mass"
MT_Left_Testis <- MiceTable %>% select(Mouse_ID, Left_Testis, Left_Testis1, Left_Testis_mass) %>% pivot_longer(names_to = "Temp",  values_to = "Left_Testis", cols = c(Left_Testis, Left_Testis1, Left_Testis_mass)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Left_Testis) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Left_Testis$Mouse_ID[duplicated(MT_Left_Testis$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Left_Testis) %>% select(-c(Left_Testis1, Left_Testis_mass))
rm(MT_Left_Testis)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Right_Testis == "Right_Testis", "Right_Testis1", ("Right_testis"), "Right_Testis_mass"
MT_Right_Testis <- MiceTable %>% select(Mouse_ID, Right_Testis, Right_Testis1, Right_Testis_mass) %>% pivot_longer(names_to = "Temp",  values_to = "Right_Testis", cols = c(Right_Testis, Right_Testis1, Right_Testis_mass)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Right_Testis) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Right_Testis$Mouse_ID[duplicated(MT_Right_Testis$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Right_Testis) %>% select(-c(Right_Testis1, Right_Testis_mass))
rm(MT_Right_Testis)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Tail_Length == "Tail_length", "Tail_length1"
MT_Tail_Length <- MiceTable %>% select(Mouse_ID, Tail_length, Tail_length1) %>% pivot_longer(names_to = "Temp",  values_to = "Tail_Length", cols = c(Tail_length, Tail_length1)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Tail_Length) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Tail_Length$Mouse_ID[duplicated(MT_Tail_Length$Mouse_ID)]   
## join
MiceTable <- full_join(MiceTable, MT_Tail_Length) %>% select(-c(Tail_length, Tail_length1))
rm(MT_Tail_Length)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Trap_Date == "Capture", "Trap_date"
MT_Trap_Date <- MiceTable %>% select(Mouse_ID, Capture, Trap_date) %>% pivot_longer(names_to = "Temp",  values_to = "Trap_Date", cols = c(Capture, Trap_date)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Trap_Date) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Trap_Date$Mouse_ID[duplicated(MT_Trap_Date$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Trap_Date) %>% select(-c(Capture, Trap_date))
rm(MT_Trap_Date)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Trichuris == "Trichuris" "Trichuris_muris"
MT_Trichuris <- MiceTable %>% select(Mouse_ID, Trichuris, Trichuris_muris) %>% pivot_longer(names_to = "Temp",  values_to = "Trichuris", cols = c(Trichuris, Trichuris_muris)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Trichuris) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Trichuris$Mouse_ID[duplicated(MT_Trichuris$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Trichuris) %>% select(-c(Trichuris_muris))
rm(MT_Trichuris)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## Year == "Year", "year"
MT_Year <- MiceTable %>% select(Mouse_ID, Year, year) %>% pivot_longer(names_to = "Temp",  values_to = "Year", cols = c(Year, year)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Year)  %>% distinct(Mouse_ID, .keep_all = T) 
MT_Year$Mouse_ID[duplicated(MT_Year$Mouse_ID)]    
## join
MiceTable <- full_join(MiceTable, MT_Year) %>% select(-c(year))
rm(MT_Year)
## eliminate duplicates in MiceTable
MiceTable$Mouse_ID[duplicated(MiceTable$Mouse_ID)]
MiceTable <- MiceTable %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


#### COLUMN INPUT CORRECTION ##################################################
## Sex
MiceTable <- MiceTable %>% mutate(SEX = replace(SEX, SEX == "female", "F"),
                                  SEX = replace(SEX, SEX == "male", "M"))
## Status
MiceTable <- MiceTable %>% mutate(Status = replace(Status, Status == "(pregnant)", "pregnant"),
                                  Status = replace(Status, Status == "(young)", "young"))

## State
MiceTable <- MiceTable %>% mutate(State = replace(State, State == "Germany", "DE"),
                                  State = replace(State, State == "D", "DE"),
                                  State = replace(State, State == "Poland", "PL"))

MiceTable <- MiceTable %>% mutate(Multiple_Mice_per_Box = ifelse(Mouse_ID %in% c("AA_0514", "AA_0515", "AA_0513","AA_0349", "AA_0454", "ZZ_0037", "ZZ_0038"), TRUE, FALSE))

## remove useless Samples
# Wildpark Schorfheide (not needed, test) + apodemus caught in 2016
wsh <- c(paste0("AA_000", 1:9), paste0("AA_00", 10:46))
apd <- c("A_0001", "A_0002", "A_0003")
# useless info
useless <- c(wsh, apd)
MiceTable <- MiceTable[!(MiceTable$Mouse_ID %in% useless),]


MiceTable <- MiceTable[, colnames(MiceTable) %in% c(basics, gen.loci, dissection.cols, parasite.cols, EimGeno.cols, oocyst.cols), ]      

write.csv(MiceTable, "MiceTable_fullColAdjust.csv")

## Body_Weight
MiceTable$Body_Weight[!is.na(MiceTable$Body_Weight) & MiceTable$Body_Weight > 100] <-
  MiceTable$Body_Weight[!is.na(MiceTable$Body_Weight) & MiceTable$Body_Weight > 100] / 1000

## Body_Length
MiceTable$Body_Length[which(MiceTable$Body_Length < 20)] <- 
  MiceTable$Body_Length[which(MiceTable$Body_Length < 20)] * 10

## BCI == Body Condition Index as log(body_weight)/ log(body_length), (Hayes et al. 2014)
MiceTable$BCI <- log(MiceTable$Body_Weight) / log(MiceTable$Body_Length)

## add "Farm" variable for better Localization
MiceTable$Farm <- paste0(MiceTable$Longitude, MiceTable$Latitude)

## remove empty rows
MiceTable <- MiceTable[!is.na(MiceTable$Mouse_ID),]
## remove duplicated rows
MiceTable <- unique(MiceTable)


## introduce "Missing Mice", Mouse_IDs for which we don't have HI data yet
missingMice = MiceTable$Mouse_ID[is.na(MiceTable$HI)]


#### MANUAL CORRECTION -- INDIVIDUAL SAMPLES AND ###############################


## 26 June 2018, add Jarda new csv


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





## replace Year of AA_0330 from "2014" to "2017"


## c) check missing Data
# vis_miss(Dis2015)
## fill in the missing Transect Data
## Polish Samples are missing Transect Info
MiceTable$Mouse_ID[is.na(MiceTable$Transect)]
MiceTable_PL      <- MiceTable %>% filter(is.na(Transect)) %>% mutate(Transect = "Allo_PL")
MiceTable_DE      <- MiceTable %>% filter(!is.na(Transect)) %>% mutate(Transect = "Allo_PL")
MiceTable         <- bind_rows(MiceTable_PL, MiceTable_DE) %>% arrange(Mouse_ID)
rm(MiceTable_DE)
rm(MiceTable_PL)
## some Adresses were not supplied with a Code
MiceTable$Mouse_ID[is.na(MiceTable$Code)]
MiceTable_Code_known <- MiceTable %>% filter(!is.na(Code))
MiceTable_Gola       <- MiceTable %>% filter(PIN %in% c(3469:3471)) %>% mutate(Code = "GOLA48")
MiceTable_Recl       <- MiceTable %>% filter(PIN %in% c(3472,3473)) %>% mutate(Code = "RECL8")
MiceTable_Trzes      <- MiceTable %>% filter(PIN %in% c(3474,3475)) %>% mutate(Code = "TRZES17")
MiceTable <- bind_rows(MiceTable_Code_known, MiceTable_Gola, MiceTable_Recl, MiceTable_Trzes) %>% arrange(Mouse_ID)
rm(MiceTable_Code_known)
rm(MiceTable_Gola)
rm(MiceTable_Recl)
rm(MiceTable_Trzes)

## replace input of special case AA_0029: "0.181 both"
MiceTable_AA_0029 <- MiceTable %>%
  filter(Mouse_ID == "AA_0029") %>%
  mutate(Left_Testis = replace(Left_Testis, Left_Testis == "NA", 0.181),
         Right_Testis = replace(Right_Testis, Right_Testis == "NA", 0.181))
MiceTable_not_AA_0029 <- MiceTable %>% filter(Mouse_ID != "AA_0029")

MiceTable <- full_join(MiceTable_AA_0029, MiceTable_not_AA_0029) %>% select(-Testis)
rm(MiceTable_Sep)
rm(MiceTable_AA_0029)
rm(MiceTable_not_AA_0029)
MiceTable$Left_Testis   <- as.double(MiceTable$Left_Testis)
MiceTable$Right_Testis  <- as.double(MiceTable$Right_Testis)

## clean AA_0014 tail_length input
MiceTable_AA_0014     <- MiceTable %>% filter(Mouse_ID == "AA_0014") %>% mutate(Tail_Length = replace(Tail_Length, Tail_Length == "NA", 55))
MiceTable_not_AA_0014 <- MiceTable %>% filter(Mouse_ID != "AA_0014")
MiceTable             <- full_join(MiceTable_AA_0014, MiceTable_not_AA_0014) %>% arrange(Mouse_ID)
MiceTable$Tail_Length <- as.double(MiceTable$Tail_Length)
rm(MiceTable_AA_0014)
rm(MiceTable_not_AA_0014)


## Spleen AA_0286 Input: "<0.01"
MiceTable_AA_0286 <- MiceTable %>% filter(Mouse_ID == "AA_0286") %>% mutate(Spleen = replace(Spleen, Spleen == "<0.01", 0.01))
MiceTable_not_AA_0286 <- MiceTable %>% filter(Mouse_ID != "AA_0286")
MiceTable <- full_join(MiceTable_AA_0286, MiceTable_not_AA_0286) %>% arrange(Mouse_ID)
rm(MiceTable_AA_0286)
rm(MiceTable_not_AA_0286)
## Feces Weight AA_0514: "With 513 and 515"
## Feces Weight AA_0515: "With 513 and 514"
## Feces Weight AA_0513: "0.849 BUT same cage as 514 and 515"
## Feces Weight AA_0349: "2 mice same cage ^!! 0.49"
## Feces Weight AA_0454: "0.26 but in the same cage as AA_0463!!!!"
## Feces Weight ZZ_0037: "0.82 but mixed with ZZ_0038"
## Feces Weight ZZ_0038: "0.82 but mixed with ZZ_0037"
MiceTable_wrong_Fec   <- MiceTable %>% filter(Mouse_ID %in% c("AA_0514", "AA_0515", "AA_0513","AA_0349", "AA_0454", "ZZ_0037", "ZZ_0038")) %>%
  mutate(Feces_Weight = replace(Feces_Weight, Feces_Weight == "NA", 0.849),
         Feces_Weight = replace(Feces_Weight, Feces_Weight == "NA", 0.849),
         Feces_Weight = replace(Feces_Weight, Feces_Weight == "NA", 0.849),
         Feces_Weight = replace(Feces_Weight, Feces_Weight == "NA", 0.49),
         Feces_Weight = replace(Feces_Weight, Feces_Weight == "NA", 0.26),
         Feces_Weight = replace(Feces_Weight, Feces_Weight == "NA", 0.82),
         Feces_Weight = replace(Feces_Weight, Feces_Weight == "NA", 0.82))
MiceTable_correct_Fec <- MiceTable %>% filter(Mouse_ID != "AA_0514", Mouse_ID != "AA_0515", Mouse_ID != "AA_0513", Mouse_ID != "AA_0349", Mouse_ID != "AA_0454", Mouse_ID != "ZZ_0037",  Mouse_ID != "ZZ_0038")
MiceTable <- full_join(MiceTable_wrong_Fec, MiceTable_correct_Fec) %>% arrange(Mouse_ID)
MiceTable$Feces_Weight <- as.double(MiceTable$Feces_Weight)
rm(MiceTable_wrong_Fec)
rm(MiceTable_correct_Fec)
## Tail_length in Samples AA_0245, AA_0279, AA_0363, AA_0375 --> fix or turn into NA?

## fix Ectoparasites_Logical, Sample AA_0246 Ectoparasites_Logical == "TRUE (collected)" --> TRUE

## fix Fleas, Sample AA_0246 Ectoparasites_Logical == "TRUE (collected)" --> TRUE




