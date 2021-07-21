library(dplyr)
library(tidyverse)
library(visdat)
library(ggplot2)
library(stringr)
library(data.table)
library(reshape)
library(ggmap)
library(lubridate)

## start with Data from Alice from previous Years
    Alice <- read.csv("/Users/FinnLo/Documents/Programming/R/HZ_SC_and_Raw_Data/Mouse_Eimeria_Field/data_products/MiceTableMusAliceArticle.csv")
    ## shorten the HI_NLoci Column to ensure it's integer for comparison with new data and visualization
    Alice$HI_NLoci <- gsub(pattern = "HI ", replacement = "", x = Alice$HI_NLoci)
    Alice$HI_NLoci <- as.integer(Alice$HI_NLoci)
    Alice$Mouse_ID <- gsub(pattern = "Sk3173", replacement = "SK_3173", x = Alice$Mouse_ID)
    
    ## make sure there are no Mouse_ID duplicates in the Data
    ## if so, then group by Mouse_ID, compare the (two) sets of the same Mouse_ID
    ## and fill missing columns, and delete the duplicates (keep distinct Mouse_IDs)
    Alice <- Alice %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

    ## Body_weight
    Alice$Body_weight[!is.na(Alice$Body_weight) & Alice$Body_weight > 100] <- Alice$Body_weight[!is.na(Alice$Body_weight) & Alice$Body_weight > 100] / 1000
    
    ## Body_length
    Alice$Body_length[which(Alice$Body_length < 20)] <- Alice$Body_length[which(Alice$Body_length < 20)] * 10
    
    ## BCI == Body Condition Index as log(body_weight)/ log(body_length), (Hayes et al. 2014)
    Alice$BCI <- log(Alice$Body_weight) / log(Alice$Body_length)
    
    # Correct wrong Long/Lat: 1. wrong multiplicative factor
    Alice$Longitude[Alice$Longitude > 100 & !is.na(Alice$Longitude)] <- Alice$Longitude[Alice$Longitude > 100 & !is.na(Alice$Longitude)] / 1000
    Alice$Latitude[Alice$Latitude > 100 &  !is.na(Alice$Latitude)] <- Alice$Latitude[Alice$Latitude > 100 &  !is.na(Alice$Latitude)] / 1000
    
    # Correct wrong Long/Lat: 2. wrong lat/lon inversion
    Alice$Longitude.temp <- Alice$Longitude
    Alice$Longitude[!is.na(Alice$Longitude) & Alice$Longitude >= 30] = Alice$Latitude[!is.na(Alice$Longitude) & Alice$Longitude >= 30]
    Alice$Latitude[!is.na(Alice$Latitude) & Alice$Latitude <= 20] = Alice$Longitude.temp[!is.na(Alice$Latitude) & Alice$Latitude <= 20]
    
    Alice = Alice[ , -which(names(Alice) %in% c("Longitude.temp"))]
    
    
    ## add "Farm" variable for better Localization
    Alice$Farm <- paste0(Alice$Longitude, Alice$Latitude, sep = " ")
    
    ## remove empty rows and remove duplicated rows
    Alice <- Alice[!is.na(Alice$Mouse_ID),]
    Alice <- unique(Alice)
    
    missingMice = Alice$Mouse_ID[is.na(Alice$HI)]
    
    
    
## add new Jarda Data
    newCSV = read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/EmanuelData.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)

    ## change Column names to standardizes names for merging
    setnames(newCSV, old = c("PIN", "X_Longit", "Y_Latit"), new = c("Mouse_ID", "Longitude", "Latitude"), skip_absent = T)
    newCSV$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_", x = newCSV$Mouse_ID)
    newCSV$Mouse_ID <- gsub(pattern = "Sk3173", replacement = "SK_3173", x = newCSV$Mouse_ID)
    
    
    
    ## data we want to include from Jarda Data
    new_Alice <-  full_join(Alice, newCSV)
    #new_Alice[,order(colnames(new_Alice))]
    
    ## fill, then remove duplicate rows
    new_Alice <- new_Alice %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 
    
    
    # correct Year
    new_Alice$Year[ new_Alice$Mouse_ID %in% c("SK_2903", "SK_2904")] <- 2014
    new_Alice$Year[ new_Alice$Mouse_ID %in% c("AA_0330", "AA_0450", "AA_0451", "AA_0452")] <- 2017
    
    
    # Add Eimeria information
    ## flotation --> do we have more recent flotation data???
    flot <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Eimeria_detection/FINALOocysts2015to2017.csv")
    ## how many neubauer cells were counted 
    flot$Ncells <- apply(flot[paste0("N_oocysts_sq", 1:8)], 1, function(x) sum(!is.na(x)))
    
    flot$OPG <- rowSums(flot[,paste0("N_oocysts_sq", 1:8)], na.rm = T) / flot$Ncells * 10000 / (flot$PBS_dil_in_mL * flot$Feces_g)
    
    new_Alice = full_join(new_Alice, flot)
    new_Alice <- new_Alice %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 
    
    
    
    ## Eimeria qPCr Data
    Eim_qpcr <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Eimeria_detection/FINALqpcrData_2016_2017_threshold5.csv")

    new_Alice <- full_join(new_Alice, Eim_qpcr[c("Mouse_ID", "delta_ct_ilwe_MminusE", "delta_ct_cewe_MminusE", "observer")])
      
    
    ## Species Identification
    Eim_species <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Eimeria_detection/Eimeria_species_assignment_14_17.csv")
    names(Eim_species)[names(Eim_species) %in% "Eim_Species"] <- "eimeriaSpecies"
    ## remove space error
    Eim_species$Mouse_ID <- gsub(" ", "", as.character(Eim_species$Mouse_ID))
    
    
    ## join
    new_Alice <- left_join(new_Alice, Eim_species[c("Mouse_ID", "n18S_Seq", "COI_Seq", "ORF470_Seq", "Species")])
     
    
    # clean Year
    new_Alice[grep("A_00", new_Alice$Mouse_ID),"Year"] <- 2016
    
    ## TO BE CHANGED WITH TISSUE RESULTS
    new_Alice$eimeriaSpecies <- as.character(new_Alice$eimeriaSpecies)
    new_Alice$eimeriaSpecies[new_Alice$Mouse_ID %in% "AA_0111"] <- "Double_ferrisi_vermiformis"
    new_Alice$eimeriaSpecies[new_Alice$Mouse_ID %in% "AA_0244"] <- "Double_tbd"
    new_Alice$eimeriaSpecies[new_Alice$Mouse_ID %in% "AA_0245"] <- "Double_tbd"
    new_Alice$eimeriaSpecies[new_Alice$Mouse_ID %in% "AA_0436"] <- "Double_tbd"
    new_Alice$eimeriaSpecies[new_Alice$Mouse_ID %in% "AA_0497"] <- "Double_tbd"
    new_Alice$eimeriaSpecies <- as.factor(new_Alice$eimeriaSpecies)
    
    ## find, fill, remove Duplicates
    new_Alice$Mouse_ID[duplicated(new_Alice$Mouse_ID)]
    new_Alice <- new_Alice %>% arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 
    
    
    
    
## add Crypto DNA Extraction Data
    ILWE_DNA_Extraction   <- read.csv("WD_07_21/DNA_Extraction_ILWE_2018_2019.csv")
    ILWE_DNA_Extraction   <- ILWE_DNA_Extraction[!names(ILWE_DNA_Extraction) == "X"]
    
## add Crypto_qPCR Data
    Crypto_qPCR  <- read.csv("WD_07_21/qPCR_Data_MouseID_forJoin.csv")
    Crypto_qPCR  <- Crypto_qPCR[!names(Crypto_qPCR) == "X"]
    Crypto_qPCR$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_", x = Crypto_qPCR$Mouse_ID)
    Crypto_qPCR$Mouse_ID <- gsub(pattern = "SK__", replacement = "SK_", x = Crypto_qPCR$Mouse_ID)
    Crypto_qPCR$Mouse_ID[duplicated(Crypto_qPCR$Mouse_ID)]
    Crypto_qPCR <- Crypto_qPCR %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 
    
    
## JOINS
    ## Alice and new Jarda Info
    ## contains Eimeria data already
    new_Alice
    
    ## introduce the Data that was generated for Crypto
    ## aka DNA Extraction values, and Ct values (from Eppendorf and ABI Machines) 
    ## first join DNA Extraction Data and qPCR Data
    Crypto <- full_join(Crypto_qPCR, ILWE_DNA_Extraction)
    Crypto$Mouse_ID[duplicated(Crypto$Mouse_ID)]
    Crypto <- Crypto %>% arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 
    
    ## add Crypto Data to the "new_Alice" Table
    new_Alice <- left_join(new_Alice, Crypto[c("Mouse_ID", "ILWE_DNA_Content_ng.microliter", "Ct_mean", "Oocyst_Predict")])
    new_Alice$Mouse_ID[duplicated(new_Alice$Mouse_ID)]
    new_Alice <- new_Alice %>% arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 
    
  

#### COLUMN CORRECTION #########################################################

## Address
    MT_Address <- new_Alice %>% select(Mouse_ID, Address, Locality)
    MT_Address <- MT_Address %>% pivot_longer(names_to = "Temp", values_to = "Address", cols = c(Address, Locality)) %>% 
      arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Address) %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Address$Mouse_ID[duplicated(MT_Address$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Address) %>% select(-c(Locality))
    rm(MT_Address)
    ## eliminate duplicates in new_Alice
    new_Alice$Mouse_ID[duplicated(new_Alice$Mouse_ID)]
    new_Alice <- new_Alice %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 
    

## Body_Weight
    MT_Body_Weight <- new_Alice %>% select(Mouse_ID, BW, Body_weight)
    MT_Body_Weight <- MT_Body_Weight %>% pivot_longer(names_to = "Temp", values_to = "Body_Weight", cols = c(Body_weight, BW)) %>% 
      arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Body_Weight) %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Body_Weight$Mouse_ID[duplicated(MT_Body_Weight$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Body_Weight) %>% select(-c(Body_weight, BW))
    rm(MT_Body_Weight)
   
    
## Body_Length == "Body_length", "Body_length1" -----> correct one super low value 8.9 or something.. --> 8.9
    MT_Body_Length <- new_Alice %>% select(Mouse_ID, L, Body_length)
    MT_Body_Length <- MT_Body_Length %>% pivot_longer(names_to = "Temp", values_to = "Body_Length", cols = c(L, Body_length)) %>% 
      arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Body_Length) %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Body_Length$Mouse_ID[duplicated(MT_Body_Length$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Body_Length) %>% select(-c(Body_length))
    rm(MT_Body_Length)

## Ectoparasites == "Ectoparasites", "Ectoparasites_Logical"
    new_Alice$Ectoparasites_Logical <- as.logical(new_Alice$Ectoparasites)
    MT_Ectoparasites_Logical <- new_Alice %>% select(Mouse_ID, Ectoparasites) %>% pivot_longer(names_to = "Temp", values_to = "Ectoparasites_Logical", cols = c(Ectoparasites)) %>% 
      arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Ectoparasites_Logical) %>% distinct(Mouse_ID, .keep_all = T) 
    new_Alice <- full_join(new_Alice, MT_Ectoparasites_Logical) %>% select(-c(Ectoparasites))
    rm(MT_Ectoparasites_Logical)
    
## Epididymis
    ## Left_Epididymis == "Left_epididymis", "Left.epididymis.weight"
    new_Alice$Left.epididymis.weight <- as.double(new_Alice$Left.epididymis.weight)
    MT_Left_Epididymis <- new_Alice %>% select(Mouse_ID, Left.epididymis.weight) %>% pivot_longer(names_to = "Temp", values_to = "Left_Epididymis", cols = c(Left.epididymis.weight)) %>% 
      arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Left_Epididymis) %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Left_Epididymis$Mouse_ID[duplicated(MT_Left_Epididymis$Mouse_ID)]  
    ## join
    new_Alice <- full_join(new_Alice, MT_Left_Epididymis) %>% select(-c(Left.epididymis.weight))
    rm(MT_Left_Epididymis)
   
## Feces_weight == "Feces_weight", "Feces_g", Feces
new_Alice$Feces_weight <- as.double(new_Alice$Feces_weight)
    MT_Feces_Weight <- new_Alice %>% select(Mouse_ID, Feces_weight, Feces_g) %>% pivot_longer(names_to = "Temp",  values_to = "Feces_Weight", cols = c(Feces_weight, Feces_g)) %>%
      arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Feces_Weight) %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Feces_Weight$Mouse_ID[duplicated(MT_Feces_Weight$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Feces_Weight) %>% select(-c(Feces_g, Feces_weight))
    rm(MT_Feces_Weight)
    
## Fleas == "Flea", "Fleas"
    new_Alice <- new_Alice %>% mutate(Fleas_Count = ifelse(Flea %in% c("0", "1", "2", "3", "4", "5", "6", "9", "11", "12"), as.numeric(Flea),
                                             ifelse(NA)),
                                      Fleas_Logical = ifelse(Flea == "fleas", TRUE,
                                                             ifelse(Flea == "TRUE", TRUE,
                                                                    ifelse(Flea == "TRUE (collected)", TRUE,
                                                                           ifelse(Fleas_Count > 0, TRUE,
                                                                                  ifelse(Fleas_Count == 0, FALSE,
                                                                                         ifelse(Flea == "FALSE", FALSE,
                                                                                                ifelse(NA))))))))
    new_Alice <- new_Alice %>% select(-Flea)

    
## Head_Taken == "Head_taken", "Head.taken."

    MT_Head_Taken <- new_Alice %>% select(Mouse_ID, Head.taken.) %>%  pivot_longer(names_to = "Temp", values_to = "Head_Taken", cols = c(Head.taken.)) %>% 
      arrange(Mouse_ID) %>% group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Head_Taken)
    MT_Head_Taken <- MT_Head_Taken %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Head_Taken$Mouse_ID[duplicated(MT_Head_Taken$Mouse_ID)]   
    ## join
    new_Alice <- full_join(new_Alice, MT_Head_Taken) %>% select(-c(Head.taken.))
    rm(MT_Head_Taken)
   
    
## Heterakis == "Heterakis", "Heterakis_spumosa"
    MT_Heterakis <- new_Alice %>% select(Mouse_ID, Heterakis, Heterakis_spumosa) %>%pivot_longer(names_to = "Temp", values_to = "Heterakis", cols = c(Heterakis, Heterakis_spumosa)) %>% 
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Heterakis)
    MT_Heterakis <- MT_Heterakis %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Heterakis$Mouse_ID[duplicated(MT_Heterakis$Mouse_ID)]  
    ## join
    new_Alice <- full_join(new_Alice, MT_Heterakis) %>% select(-c(Heterakis_spumosa))
    rm(MT_Heterakis)
    
 
## Hymenolepis == "Hymenolepis", "Hymenolepis_microstoma", "Hymenolepis_diminiuta"
    MT_Hymenolepis <- new_Alice %>% select(Mouse_ID, Hymenolepis, Hymenolepis_diminiuta, Hymenolepis_microstoma)
    MT_Hymenolepis <- MT_Hymenolepis %>% pivot_longer(names_to = "Temp",  values_to = "Hymenolepis", cols = c(Hymenolepis, Hymenolepis_diminiuta, Hymenolepis_microstoma)) %>%
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Hymenolepis)
    MT_Hymenolepis <- MT_Hymenolepis %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Hymenolepis$Mouse_ID[duplicated(MT_Hymenolepis$Mouse_ID)]   
    ## join
    new_Alice <- full_join(new_Alice, MT_Hymenolepis)
    rm(MT_Hymenolepis)
  
    
## Latitude == "Latitude", "longitude" ***** MIX UP WITH LONG *****
    MT_Latitude <- new_Alice %>% select(Mouse_ID, Latitude, longitude) %>% pivot_longer(names_to = "Temp",  values_to = "Latitude", cols = c(Latitude, longitude)) %>%
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%   fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Latitude)
    MT_Latitude <- MT_Latitude %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Latitude$Mouse_ID[duplicated(MT_Latitude$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Latitude) %>% select(-c(longitude))
    rm(MT_Latitude)
    

## Longitude == "Longitude", "latitude" ***** MIX UP WITH LAT *****
    MT_Longitude <- new_Alice %>% select(Mouse_ID, Longitude, latitude)  %>% pivot_longer(names_to = "Temp", values_to = "Longitude", cols = c(Longitude, latitude)) %>% 
      arrange(Mouse_ID) %>% group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Longitude)
    MT_Longitude <- MT_Longitude %>% distinct(Mouse_ID, .keep_all = T) 
    ## join
    new_Alice <- full_join(new_Alice, MT_Longitude) %>% select(-c(latitude))
    rm(MT_Longitude)
 
    
## Liver == "Liver", "Liver_mass"
    MT_Liver <- new_Alice %>% select(Mouse_ID, Liver_mass) %>% pivot_longer(names_to = "Temp",  values_to = "Liver", cols = c(Liver_mass)) %>%  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Liver)
    MT_Liver <- MT_Liver %>% distinct(Mouse_ID, .keep_all = T) 
    ## join
    new_Alice <- full_join(new_Alice, MT_Liver) %>% select(-c(Liver_mass))
    rm(MT_Liver)


## Mastophorus == "Mastophorus", "Mastophorus_muris", "Mastaphorus"
    MT_Mastophorus <- new_Alice %>% select(Mouse_ID, Mastophorus, Mastophorus_muris) %>% pivot_longer(names_to = "Temp",  values_to = "Mastophorus", cols = c(Mastophorus, Mastophorus_muris)) %>% 
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>% select(Mouse_ID, Mastophorus)
    MT_Mastophorus <- MT_Mastophorus %>% distinct(Mouse_ID, .keep_all = T) 
    ## join
    new_Alice <- full_join(new_Alice, MT_Mastophorus) %>% select(-c(Mastophorus_muris))
    rm(MT_Mastophorus)
 
    
## Notes == "comments", "Note", "Notes", ("Embryo")
    MT_Notes <- new_Alice %>% select(Mouse_ID, Note, Notes, comments)  %>% pivot_longer(names_to = "Temp",  values_to = "Notes", cols = c(Note, Notes, comments)) %>%  
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Notes)
    MT_Notes <- MT_Notes %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Notes$Mouse_ID[duplicated(MT_Notes$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Notes) %>% select(-c(Note, comments))
    rm(MT_Notes)
    ## eliminate duplicates in new_Alice
    new_Alice$Mouse_ID[duplicated(new_Alice$Mouse_ID)]
    new_Alice <- new_Alice %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Ovaria ==
## Right_Ovarium_Weight == "Right_ovarium", Right.ovarium.weight"
    MT_Right_Ovarium_Weight <- new_Alice %>% select(Mouse_ID, Right.ovarium.weight) %>% pivot_longer(names_to = "Temp",  values_to = "Right_Ovarium_Weight", cols = c(Right.ovarium.weight)) %>%  
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Right_Ovarium_Weight) %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Right_Ovarium_Weight$Mouse_ID[duplicated(MT_Right_Ovarium_Weight$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Right_Ovarium_Weight) %>% select(-c(Right.ovarium.weight))
    rm(MT_Right_Ovarium_Weight)
    

## Left_Ovarium_Weight == "Left_ovarium", "Left.ovarium.weight"
    MT_Left_Ovarium_Weight <- new_Alice %>% select(Mouse_ID, Left.ovarium.weight) %>% pivot_longer(names_to = "Temp",  values_to = "Left_Ovarium_Weight", cols = c(Left.ovarium.weight)) %>%  
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Left_Ovarium_Weight) %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Left_Ovarium_Weight$Mouse_ID[duplicated(MT_Left_Ovarium_Weight$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Left_Ovarium_Weight) %>% select(-c(Left.ovarium.weight))
    rm(MT_Left_Ovarium_Weight)
  

## Region == "Region", "REGion"
    MT_Region <- new_Alice %>% select(Mouse_ID, Region, REGion) %>% pivot_longer(names_to = "Temp",  values_to = "Region", cols = c(Region, REGion)) %>%
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Region)
    MT_Region <- MT_Region %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Region$Mouse_ID[duplicated(MT_Region$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Region) %>% select(-c(REGion))
    rm(MT_Region)
    ## eliminate duplicates in new_Alice
    new_Alice$Mouse_ID[duplicated(new_Alice$Mouse_ID)]
    new_Alice <- new_Alice %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

    
## Seminal_Vesicles_Weight == SemVes, Seminal.vesicle.weight, Seminal_Vesicles_Weight
    MT_Seminal_Vesicles_Weight <- new_Alice %>% select(Mouse_ID, SemVes, Seminal.vesicle.weight) %>% pivot_longer(names_to = "Temp",  values_to = "Seminal_Vesicles_Weight", cols = c(SemVes, Seminal.vesicle.weight)) %>%
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Seminal_Vesicles_Weight)
    MT_Seminal_Vesicles_Weight <- MT_Seminal_Vesicles_Weight %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Seminal_Vesicles_Weight$Mouse_ID[duplicated(MT_Seminal_Vesicles_Weight$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Seminal_Vesicles_Weight) %>% select(-c(SemVes, Seminal.vesicle.weight))
    rm(MT_Seminal_Vesicles_Weight)


## Spleen == "Spleen", "Spleen_mass"
    new_Alice$Spleen_mass <- as.double(new_Alice$Spleen_mass)
    MT_Spleen <- new_Alice %>% select(Mouse_ID, Spleen, Spleen_mass)  %>% pivot_longer(names_to = "Temp",  values_to = "Spleen",cols = c(Spleen, Spleen_mass)) %>% 
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Spleen) %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Spleen$Mouse_ID[duplicated(MT_Spleen$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Spleen) %>% select(-c(Spleen_mass))
    rm(MT_Spleen)
    ## eliminate duplicates in new_Alice
    new_Alice$Mouse_ID[duplicated(new_Alice$Mouse_ID)]
    new_Alice <- new_Alice %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

    
## Taenia == "Taenia_martis", "Taenia_taeniformis", "Catenotaenia_pusilla", "Cysticercus"
    MT_Taenia <- new_Alice %>% select(Mouse_ID, Taenia_martis, Taenia_taeniformis, Catenotaenia_pusilla, Cysticercus, Taenia)
    MT_Taenia <- MT_Taenia %>%
      pivot_longer(names_to = "Temp", 
                   values_to = "Taenia",
                   cols = c(Taenia_martis, Taenia_taeniformis, Catenotaenia_pusilla, Cysticercus, Taenia))
    MT_Taenia <- MT_Taenia %>%  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Taenia)
    MT_Taenia <- MT_Taenia %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Taenia$Mouse_ID[duplicated(MT_Taenia$Mouse_ID)]    
    new_Alice <- full_join(new_Alice, MT_Taenia) %>% select(-c(Taenia_martis, Taenia_taeniformis, Catenotaenia_pusilla, Cysticercus))
    rm(MT_Taenia)
    ## eliminate duplicates in new_Alice
    new_Alice$Mouse_ID[duplicated(new_Alice$Mouse_ID)]
    new_Alice <- new_Alice %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 
    

## Testis_mass Separation
    ## Wrong data input for "Testis_mass": instead of individual Left_Testis or 
    ## Right_Testis data, a combination of "Left_Testis/Right_Testis" was supplied
    ## Separation of that Column into 2 Columns, original "Testis"  col. discarded
    new_Alice_Sep <- new_Alice %>% select(Mouse_ID, Testis_mass) %>% filter(Testis_mass != "NA")
    new_Alice_Sep <- setDT(new_Alice_Sep)[, paste0("Testis_mass", 1:2) := tstrsplit(Testis_mass, "/")]
    setnames(new_Alice_Sep, old = c("Testis_mass1", "Testis_mass2"), new = c("Left_Testis1", "Right_Testis1"), skip_absent = T)
    new_Alice_Sep$Left_Testis1 <- as.double(new_Alice_Sep$Left_Testis1)
    new_Alice_Sep$Right_Testis1 <- as.double(new_Alice_Sep$Right_Testis1)
    ## join
    new_Alice <- full_join(new_Alice, new_Alice_Sep) %>% select(-Testis_mass)
    rm(new_Alice_Sep)
    
## Testes == "Testes"
    ## Left_Testis  == "Left_Testis1", "Left_Testis_mass"
    MT_Left_Testis <- new_Alice %>% select(Mouse_ID, Left_Testis1, Left_Testis_mass) %>% pivot_longer(names_to = "Temp",  values_to = "Left_Testis", cols = c(Left_Testis1, Left_Testis_mass)) %>% 
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Left_Testis) %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Left_Testis$Mouse_ID[duplicated(MT_Left_Testis$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Left_Testis) %>% select(-c(Left_Testis_mass, Left_Testis1))
    rm(MT_Left_Testis)
   
    ## Right_Testis == Right_Testis1", "Right_Testis_mass"
    MT_Right_Testis <- new_Alice %>% select(Mouse_ID, Right_Testis1, Right_Testis_mass) %>% pivot_longer(names_to = "Temp",  values_to = "Right_Testis", cols = c(Right_Testis1, Right_Testis_mass)) %>% 
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Right_Testis) %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Right_Testis$Mouse_ID[duplicated(MT_Right_Testis$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Right_Testis) %>% select(-c(Right_Testis_mass, Right_Testis1))
    rm(MT_Right_Testis)
 

## Tail_Length == "Tail_length", "LCd"
    MT_Tail_Length <- new_Alice %>% select(Mouse_ID, Tail_length, LCd) %>% pivot_longer(names_to = "Temp",  values_to = "Tail_Length", cols = c(Tail_length, LCd)) %>% 
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Tail_Length) %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Tail_Length$Mouse_ID[duplicated(MT_Tail_Length$Mouse_ID)]   
    ## join
    new_Alice <- full_join(new_Alice, MT_Tail_Length) %>% select(-c(Tail_length, LCd))
    rm(MT_Tail_Length)
    
## Trap_Date == "Capture"
    MT_Trap_Date <- new_Alice %>% select(Mouse_ID, Capture) %>% pivot_longer(names_to = "Temp",  values_to = "Trap_Date", cols = c(Capture)) %>% 
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Trap_Date) %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Trap_Date$Mouse_ID[duplicated(MT_Trap_Date$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Trap_Date) %>% select(-c(Capture))
    rm(MT_Trap_Date)


## Trichuris == "Trichuris" "Trichuris_muris"
    MT_Trichuris <- new_Alice %>% select(Mouse_ID, Trichuris, Trichuris_muris) %>% pivot_longer(names_to = "Temp",  values_to = "Trichuris", cols = c(Trichuris, Trichuris_muris)) %>% 
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Trichuris) %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Trichuris$Mouse_ID[duplicated(MT_Trichuris$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Trichuris) %>% select(-c(Trichuris_muris))
    rm(MT_Trichuris)
    ## eliminate duplicates in new_Alice
    new_Alice$Mouse_ID[duplicated(new_Alice$Mouse_ID)]
    new_Alice <- new_Alice %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

    
## Year == "Year", "year"
    MT_Year <- new_Alice %>% select(Mouse_ID, Year, year) %>% pivot_longer(names_to = "Temp",  values_to = "Year", cols = c(Year, year)) %>% 
      arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Year)  %>% distinct(Mouse_ID, .keep_all = T) 
    MT_Year$Mouse_ID[duplicated(MT_Year$Mouse_ID)]    
    ## join
    new_Alice <- full_join(new_Alice, MT_Year) %>% select(-c(year))
    rm(MT_Year)
    ## eliminate duplicates in new_Alice
    new_Alice$Mouse_ID[duplicated(new_Alice$Mouse_ID)]
    new_Alice <- new_Alice %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


#### COLUMN INPUT CORRECTION ##################################################
    ## Status
    new_Alice <- new_Alice %>% mutate(Status = replace(Status, Status == "(pregnant)", "pregnant"),
                                      Status = replace(Status, Status == "(young)", "young"))
    ## State
    new_Alice <- new_Alice %>% mutate(State = replace(State, State == "Germany", "DE"),
                                      State = replace(State, State == "D", "DE"),
                                      State = replace(State, State == "Poland", "PL"))
    ## multiple Mice per Box
    new_Alice <- new_Alice %>% mutate(Multiple_Mice_per_Box = ifelse(Mouse_ID %in% c("AA_0514", "AA_0515", "AA_0513","AA_0349", "AA_0454", "ZZ_0037", "ZZ_0038"), TRUE, FALSE))
    
    ## this would be the full data, some columns are cut down, but it's still
    ## pretty busy, so the next step would be direct column selection for 
    ## manageable working Size
    #write.csv(new_Alice, "new_Alice_MiceTable.csv")
    
#### Select Columns ############################################################
    

basics          <- c("Mouse_ID", "Transect", "Code", "Region", "Sex", "Longitude", "Latitude", "Year", "State", "HI", "HI_NLoci")

ILWE_DNA.cols <- c("ILWE_DNA_Content_ng.microliter", "ILWE_DNA_used_up")
    
Crypto.cols     <- c("Tested_by", "Machine", "Measurements", "Ct_mean", "Ct_mean_Ep", "Ct_mean_ABI", "Oocyst_Predict")


gen.loci        <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
                      "Es1", "Gpd1", "Idh1", "Mpi", "Np", "Sod1", "Es1C", "Gpd1C",
                      "Idh1C", "MpiC", "NpC", "Sod1C", "HI_NLoci",
                      "HI", "Zfy2", "SRY1", "Y")

dissection.cols <- c("Body_Weight", "Body_Length", "Tail_Length", "Spleen", 
                     "Left_Testis", "Right_Testis", "Seminal_Vesicles_Weight", 
                     "Sperm", "Left_Epididymis", "Right_Epididymis", 
                     "Arrival", "Wean", "Death", "Dissection_date", "DaysInLab",
                     "Lepid", "Uterus", "Ovaria", "Grav", "Litters", "NN", 
                     "MM", "FF", "Protocol")

parasite.cols   <- c("Aspiculuris_tetraptera", "Syphacia_obvelata", "Aspiculuris_Syphacia",
                     "Trichuris", "Taenia", "Flea", "Mix_Syphacia_Aspiculuris", "Heligmosomoides_polygurus",
                     "Heterakis","Heterakis", "Mastophorus","Hymenolepis",
                     "Ectoparasites", "Ectoparasites_Count", "Ectoparasites_Logical",
                     "Worms_presence")

oocyst.cols     <- c("counter", "Feces_Weight", "Date_count", "N_oocysts_sq1",
                     "N_oocysts_sq2", "N_oocysts_sq3",  "N_oocysts_sq4",
                     "N_oocysts_sq5", "N_oocysts_sq6", "N_oocysts_sq7",
                     "N_oocysts_sq8", "mean_neubauer", "PBS_dil_in_mL", 
                     "OPG", "Ncells")

EqPCR.cols      <- c("delta_ct_ilwe_MminusE", "delta_ct_cewe_MminusE",
                    ## from 2018 on AN IMPORTANT IMPROVEMENT!!!
                    "MC.Eimeria")

EimGeno.cols    <- c("n18S_Seq", "COI_Seq", "ORF470_Seq", "eimeriaSpecies")

    
    
## cut down columns, deselect empty columns
new_Alice <- new_Alice[, colnames(new_Alice) %in% c(basics, ILWE_DNA.cols, Crypto.cols ,gen.loci, dissection.cols, parasite.cols, EimGeno.cols, oocyst.cols), ]      
new_Alice <- new_Alice %>% select(!which(!colSums(!is.na(new_Alice))))
#vis_miss(new_Alice)


#write.csv(new_Alice, "new_Alice_MiceTable_Selected_Cols.csv")


## Select Columns that we have performed qPCR on, cut down columns
Ct_filter_AA <- new_Alice[, colnames(new_Alice) %in% c(basics, ILWE_DNA.cols, Crypto.cols ,gen.loci, dissection.cols, parasite.cols, EimGeno.cols, oocyst.cols), ] %>% 
  filter(Ct_mean != "NA")
#vis_miss(Ct_filter_AA)