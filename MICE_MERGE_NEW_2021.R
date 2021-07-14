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

#### 1) MAKE BASE TABLE ########################################################
    JardaTableNew   <-  read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/EmanuelData.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
    JardaTableNew   <- JardaTableNew[!names(JardaTableNew) == "X"]
    JardaTableNew   <- JardaTableNew[!JardaTableNew$Year %in% c(2010, 2011),]
    
## homogenize Column names, especially Mouse_ID should be universal!
    setnames(JardaTableNew,
             old = c("PIN", "X_Longit", "Y_Latit"), 
             new = c("Mouse_ID", "Longitude", "Latitude"))
    
## homogenize Mouse_ID Row Names, old Mouse_IDs == SK_ and new Mouse_IDs == AA_
    JardaTableNew$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_", x = JardaTableNew$Mouse_ID)
    JardaTableNew <- JardaTableNew %>% distinct(Mouse_ID, .keep_all = T)
    
# remove Embryos (no interest for parasitic studies)
    JardaTableNew   <- JardaTableNew[sapply(JardaTableNew$Mouse_ID, nchar) <= 7,]


    
## double-check with old data
    JardaTableOld   <-  read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/MiceTable_fullEimeriaInfos_2014to2017.csv")
    JardaTableOld    <- JardaTableOld[!names(JardaTableOld) == "X"]
    JardaTableOld    <- JardaTableOld[!JardaTableOld$Year %in% c(2010, 2011),]

## homogenize Mouse_ID Row Names, old Mouse_IDs == SK_ and new Mouse_IDs == AA_
    JardaTableOld$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_", x = JardaTableOld$Mouse_ID)
    JardaTableOld$HI_NLoci <- gsub(pattern = "HI ", replacement = "", x = JardaTableOld$HI_NLoci)
    JardaTableOld$HI_NLoci <- as.integer(JardaTableOld$HI_NLoci)

## remove Embryos (no interest for parasitic studies)
    JardaTableOld   <- JardaTableOld[sapply(JardaTableOld$Mouse_ID, nchar) <= 7,]



## COMBINED TABLE (ALICE's TABLE + NEW DATA)
    # old MiceTable (JardaTableOld)  +  new HI data (JardaTableNew)
    JardaTable <- left_join(JardaTableNew, JardaTableOld)
    JardaTable <- JardaTable[sapply(JardaTable$Mouse_ID, nchar) <= 7,]
    which(!rowSums(!is.na(JardaTable)))



## Selected Columns we want to record for our samples:
basics <- c("Mouse_ID", "Transect", "Code", "Region", "Sex", "Longitude", "Latitude", "Year")

gen.loci <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
              "Es1", "Gpd1", "Idh1", "Mpi", "Np", "Sod1", "Es1C", "Gpd1C",
              "Idh1C", "MpiC", "NpC", "Sod1C", "HI_NLoci",
              "HI")#"Zfy2", "SRY1", "Y")

dissection.cols <- c("Body_weight", "Body_length", "Tail_length", "Spleen", 
                     "Left_testis", "Right_testis", "SemVes", 
                     "Sperm")
                     #"Arrival", "Wean", "Death", "Dissection_date", "DaysInLab",
                     #"Testes","Lepid", "Seminal_vesicles", "Uterus", "Ovaria", "Grav", "Litters", "NN", 
                     #"MM", "FF", "Protocol")

parasite.cols <- c("Aspiculuris_tetraptera", "Syphacia_obvelata", "Trichuris_muris",
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


# Homogenize Column Names of JardaTable:
# in order to avoid duplicated column names, the empty columns are subtracted
JardaTable <- JardaTable %>% select(-c("Body_weight", "Body_length", "Tail_length"))

setnames(JardaTable,
         old = c("BW", "L", "LCd", "Ltestis", "Rtestis", "SemVes", "Dissection"),
         new = c("Body_weight", "Body_length", "Tail_length",
                 "Left_Testis", "Right_Testis", "Seminal_vesicles", "Dissection_date"),
         skip_absent = T)

## Check Data for Duplicated Mouse_IDs
JardaTable$Mouse_ID[duplicated(JardaTable$Mouse_ID)]

## There's a lot of empty columns at the moment.
## I prefer to empty the entire data frame and only leave columns with data
## So we are keeping the "basics", gen.loci", and "dissection.cols" as Sel_JardaTable (Selected cols)
Sel_JardaTable <- JardaTable[,colnames(JardaTable) %in% c(basics, gen.loci, dissection.cols), ]
#vis_miss(Sel_JardaTable)


# replace Year of AA_0330 from "2014" to "2017"
Sel_JardaTable_AA_0330 <- Sel_JardaTable %>% filter(Mouse_ID == "AA_0330") %>% mutate(Year = replace(Year, Year == 2014, 2017))
Sel_JardaTable_not_AA_0330 <- Sel_JardaTable %>% filter(Mouse_ID != "AA_0330")
Sel_JardaTable <- full_join(Sel_JardaTable_AA_0330, Sel_JardaTable_not_AA_0330) %>% arrange(Mouse_ID)
rm(Sel_JardaTable_AA_0330)
rm(Sel_JardaTable_not_AA_0330)


#### Complement data with previous tables ######################################

#### 2014 DATA ****
## a) Load Data: Dissection Data from 2014 and columns that need type Conversions
      Dis2014                     <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ14_Dissections_237-523.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
      Dis2014$Mouse_ID            <- paste0(Dis2014$ID, "_", Dis2014$PIN)
      Dis2014$Ectoparasites_Count <- as.integer(Dis2014$Ectoparasites_Count)

## b) Homogenize column names:
      setnames(Dis2014, 
         old = c("BW", "L", "LCd", 
                 "ASP", "SYP", "TM","Taenia.taeniformis", 
                 "Locality", "SemVes", "Ectoparasites", "Left_testis", "Right_testis"),
         new = c("Body_weight", "Body_length", "Tail_length", 
                 "Aspiculuris_tetraptera", "Syphacia_obvelata", 
                 "Trichuris_muris", "Taenia_taeniformis", 
                 "Address", "Seminal_vesicles", "Ectoparasites_Count", 
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
Sel_JardaTable1 <- full_join(Sel_JardaTable, Dis2014)
vis_miss(Sel_JardaTable1)


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

## c) check missing Data
      # vis_miss(Dis2015)
      ## fill in the missing Transect Data
      ## Polish Samples are missing Transect Info
      Dis2015$Mouse_ID[is.na(Dis2015$Transect)]
      Dis2015_PL      <- Dis2015 %>% filter(is.na(Transect)) %>% mutate(Transect = "Allo_PL")
      Dis2015_DE      <- Dis2015 %>% filter(State == "Germany")
      Dis2015         <- bind_rows(Dis2015_PL, Dis2015_DE) %>% arrange(Mouse_ID)
      rm(Dis2015_DE)
      rm(Dis2015_PL)
      # some Adresses were not supplied with a Code
      Dis2015$Mouse_ID[is.na(Dis2015$Code)]
      Dis2015_Code_known <- Dis2015 %>% filter(!is.na(Code))
      Dis2015_Gola       <- Dis2015 %>% filter(PIN %in% c(3469:3471)) %>% mutate(Code = "GOLA48")
      Dis2015_Recl       <- Dis2015 %>% filter(PIN %in% c(3472,3473)) %>% mutate(Code = "RECL8")
      Dis2015_Trzes      <- Dis2015 %>% filter(PIN %in% c(3474,3475)) %>% mutate(Code = "TRZES17")
      Dis2015 <- bind_rows(Dis2015_Code_known, Dis2015_Gola, Dis2015_Recl, Dis2015_Trzes) %>% arrange(Mouse_ID)
      rm(Dis2015_Code_known)
      rm(Dis2015_Gola)
      rm(Dis2015_Recl)
      rm(Dis2015_Trzes)

# join
Sel_JardaTable2 <- full_join(Sel_JardaTable1, Dis2015)  # vis_miss(Sel_JardaTable2)


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
Sel_JardaTable3 <- full_join(Sel_JardaTable2, Dis2016) #vis_miss(Sel_JardaTable3)

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
Sel_JardaTable4 <- full_join(Sel_JardaTable3, Dis2017)


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
Sel_JardaTable5 <- full_join(Sel_JardaTable4, Dis2018)



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
Sel_JardaTable6 <- full_join(Sel_JardaTable5, Dis2019)
glimpse(Sel_JardaTable6)
vis_miss(Sel_JardaTable6)

## check for only NA rows and delete
which(!rowSums(!is.na(Sel_JardaTable6)))




## Genotype Data from 2010-2019
Gen_2010_2019 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ10-19_Genotypes.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)






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




