library(dplyr)
library(tidyverse)
library(visdat)
library(stringr)
library(data.table)
library(reshape)


Jarda <- read.csv("EmanuelData.csv") # Email from Jarda
      Jarda$HI_NLoci <- gsub(pattern = "HI ", replacement = "", x = Jarda$HI_NLoci)
      setnames(Jarda, old = c("PIN", "X_Longit", "Y_Latit"), new = c("Mouse_ID", "Longitude", "Latitude"), skip_absent = T)
      Jarda$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_", x = Jarda$Mouse_ID)
      Jarda$Mouse_ID <- gsub(pattern = "Sk3173", replacement = "SK_3173", x = Jarda$Mouse_ID)


basics <- c("Mouse_ID", "Sex", "Longitude", "Latitude", "Year")

gen.loci <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
              "Es1", "Gpd1", "Idh1", "Mpi", "Np", "Sod1", "Es1C", "Gpd1C",
              "Idh1C", "MpiC", "NpC", "Sod1C", "Zfy2", "SRY1", "Y", "HI_NLoci",
              "HI")

parasite.cols <- c("Aspiculuris_tetraptera", "Syphacia_obvelata", "Trichuris_muris",
                   "Taenia_taeniformis", "Flea", "Mix_Syphacia_Aspiculuris",
                   "Heterakis_spumosa", "Mastophorus_muris", "Hymenolepis_microstoma",
                   "Catenotaenia_pusilla", "Cysticercus", "Ectoparasites",
                   "Worms_presence", "Hymenolepis_diminiuta", "Taenia_martis",
                   "Heligmosomoides_polygurus", "Taenia", "Aspiculuris_Syphacia",
                   "Trichuris", "Heterakis", "Mastophorus")
ILWE_DNA.cols   <- c("ILWE_DNA_Content_ng.microliter", "ILWE_DNA_used_up", 
                     "DNA_Dilution_ng.microliter", "dest_Water")

qPCR.cols       <- c("Ct_mean", "Ct_mean_Ep", "Ct_mean_ABI", 
                     "Ct_mean_1_ABI", "Ct_mean_2_ABI", "Ct_mean_3_ABI",
                     "Ct_believable", "Machine", "Measurements", "Tested_by", 
                     "qPCR_Date")

Jarda <- Jarda[, colnames(Jarda) %in% c(basics, gen.loci), ]


### add Crypto DNA Extraction Data
    ILWE_DNA_Extraction   <- read.csv("WD_07_22/DNA_Extraction_ILWE_2018_2019.csv")
    ILWE_DNA_Extraction$Mouse_ID[duplicated(ILWE_DNA_Extraction$Mouse_ID)]
    
### add Crypto_qPCR Data and filter out duplicates
    Crypto_qPCR           <- read.csv("WD_07_22/qPCR_Data_MouseID_forJoin.csv") 
    Crypto_qPCR$Mouse_ID  <- gsub(pattern = "SK", replacement = "SK_", x = Crypto_qPCR$Mouse_ID)
    Crypto_qPCR$Mouse_ID  <- gsub(pattern = "SK__", replacement = "SK_", x = Crypto_qPCR$Mouse_ID)
    
### merging qPCR and and Extraction data
    Crypto  <- full_join(Crypto_qPCR[colnames(Crypto_qPCR) %in% c(basics, qPCR.cols, ILWE_DNA.cols)], 
                         ILWE_DNA_Extraction[colnames(ILWE_DNA_Extraction) %in% c(basics, qPCR.cols, ILWE_DNA.cols)], by = "Mouse_ID")

    ## add HI Data with Jarda
    Crypto  <- merge(Crypto[colnames(Crypto) %in% c("Mouse_ID", qPCR.cols, ILWE_DNA.cols)], Jarda[colnames(Jarda) %in% c(basics, "HI")]) %>% filter(Ct_mean >= 0)
    ## includes FULLY available data, all gen.loci, all qPCR data
        #Crypto  <- full_join(Crypto[colnames(Crypto) %in% c("Mouse_ID", qPCR.cols, ILWE_DNA.cols)], Jarda[colnames(Jarda) %in% c(basics, "HI", gen.loci)]) %>% filter(Ct_mean >= 0)
        ## would include data that is missing gen.loci data + HI
    #vis_miss(Crypto, cluster = T)

    write.csv(Crypto, "Crypto_qPCR_Results.csv")
    
    
        
    
    