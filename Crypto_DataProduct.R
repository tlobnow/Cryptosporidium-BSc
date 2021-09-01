library(dplyr)
library(tidyverse)
library(visdat)
library(data.table)

# Email from Jarda (2019)
Jarda <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/EmanuelData.csv") 
      Jarda$HI_NLoci <- gsub(pattern = "HI ", replacement = "", x = Jarda$HI_NLoci)
      setnames(Jarda, old = c("PIN", "Y_Latit", "X_Longit"), new = c("Mouse_ID", "Latitude", "Longitude"), skip_absent = T)

      basics <- c("Mouse_ID", "Sex", "Latitude", "Longitude", "Year")
      
      gen.loci <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
                    "Es1", "Gpd1", "Idh1", "Mpi", "Np", "Sod1", "Es1C", "Gpd1C",
                    "Idh1C", "MpiC", "NpC", "Sod1C", "Zfy2", "SRY1", "Y", "HI_NLoci",
                    "HI")
      
      Crypto_DNA.cols   <- c("ILWE_DNA_Content_ng.microliter", "ILWE_used_up")
      
      Crypto_qPCR.cols  <- c("Ct_mean", "Ct_mean_Ep", "Ct_mean_ABI", 
                             "Machine", "Measurements", "Tested_by", 
                           "qPCR_Date", "Oocyst_Predict", "Crypto_Positive")
      
      Address.cols      <- c("Address", "Longitude")
      

      Jarda <- Jarda[, colnames(Jarda) %in% c(basics, gen.loci), ]
      
### add Crypto DNA Extraction Data
    Crypto_DNA   <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/WD_07_22/DNA_Extraction_ILWE_2018_2019.csv")
    setnames(Crypto_DNA, old = "ILWE_DNA_used_up", new = "ILWE_used_up")

### add Crypto_qPCR Data and filter out duplicates
    Crypto_qPCR           <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/qPCR_Data_MouseID_forJoin.csv") 


### calculate Oocysts with prediction model
    ABI_Best_thSC     <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Crap/ABI_Best_SC.csv")
    ABI_Best_thSC     <-  filter(ABI_Best_thSC, Ct_mean > 0)
    linear_model0     <- lm(log2(Amount_Oocysts) ~ Ct_mean, data = ABI_Best_thSC)
    Oocyst_Predict    <- 2^predict(linear_model0, newdata = Crypto_qPCR)
    Crypto_qPCR <- data.frame(Crypto_qPCR, Oocyst_Predict)
    Crypto_qPCR <- Crypto_qPCR %>%
      mutate(Oocyst_Predict = replace(Oocyst_Predict, Oocyst_Predict == "4292821751815.77", "0"))
    Crypto_qPCR$Oocyst_Predict <- as.integer(Crypto_qPCR$Oocyst_Predict)
    
## add Status (Crypto-positive or negative)
    Crypto_qPCR <- Crypto_qPCR %>% mutate(Crypto_Positive = ifelse(Ct_mean > 0, T, F))
  
### merging qPCR and and Extraction data
    Crypto_Detection  <- full_join(Crypto_qPCR[colnames(Crypto_qPCR) %in% c(basics, Crypto_qPCR.cols, Crypto_DNA.cols)], 
                         Crypto_DNA[colnames(Crypto_DNA) %in% c(basics, Crypto_qPCR.cols, Crypto_DNA.cols)], by = "Mouse_ID")

## add HI Data with Jarda
    Crypto_Detection  <- merge(Crypto_Detection[colnames(Crypto_Detection) %in% c("Mouse_ID", Crypto_qPCR.cols, Crypto_DNA.cols)], Jarda[colnames(Jarda) %in% c(basics, "HI")]) %>% filter(Ct_mean >= 0)

## add Address variable for improved mapping
    Trapping_Data     <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/HZ14-HZ19_Trapping_coodinates.csv")
    Crypto_Detection <- left_join(Crypto_Detection, Trapping_Data[colnames(Trapping_Data) %in% Address.cols])
    Crypto_Detection <- Crypto_Detection %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## add mus_caught per Location and Crypto_mus_caught per Location
    Crypto_pull <- Crypto_Detection %>% count(Longitude)%>% mutate(mus_caught = n) %>% select(Longitude, mus_caught) %>% arrange(mus_caught)
    Crypto_pull_pos <- Crypto_Detection %>% filter(Ct_mean > 0) %>% count(Longitude) %>% mutate(Crypto_mus_caught = n) %>% select(Longitude, Crypto_mus_caught) %>% arrange(Crypto_mus_caught)
    Crypto_Detection_1 <- left_join(Crypto_Detection, Crypto_pull)
    Crypto_Detection <- left_join(Crypto_Detection_1, Crypto_pull_pos)
    Crypto_Detection <- Crypto_Detection %>% replace_na(list(Crypto_mus_caught = 0))
    
## add Infection Rate per Location
    Crypto_Detection <- Crypto_Detection %>%
      mutate(Infection_Rate = Crypto_mus_caught / mus_caught)

## Top Location
    df <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Top_Locations.csv")
    Crypto_Detection <- full_join(Crypto_Detection, df, by = c("Latitude", "Longitude", "mus_caught", "Crypto_mus_caught")) %>% 
      mutate(Top_Location = Crypto_mus_caught >= 3) %>% select(-Address.x)
    setnames(Crypto_Detection, old = "Address.y", new = "Address")
    
## write csv
    write.csv(Crypto_Detection, "Crypto_Detection.csv")
    