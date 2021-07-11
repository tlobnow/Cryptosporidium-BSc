library(dplyr)
library(tidyverse)
library(visdat)
library(ggplot2)
library(stringr)


# Start with all Data ever collected in previous years and collect it in 1 .csv file

setwd("/Users/FinnLo/Documents/Programming/R/HZ_SC_and_Raw_Data")
#setwd("/Users/FinnLo/Documents/Programming/R/HZ_SC_and_Raw_Data/Cryptosporidium-BSc")

# original Data
DNA_Extraction_ILWE_2019  <-  read.csv("./Mouse_Eimeria_Field/data_input/Cryptosporidium/DNA_Extraction_ILWE_2019.csv")
DNA_Extraction_ILWE_2018  <-  read.csv("./Mouse_Eimeria_Field/data_input/Cryptosporidium/cryptoDetection2018.csv")
EmanuelData               <-  read.csv("./Mouse_Eimeria_Field/data_input/Cryptosporidium/EmanuelData.csv")
Gen16                     <-  read.csv("./Mouse_Eimeria_Field/data_input/Mouse_data/HZ16_Genotypes_47-211.csv")
Gen17                     <-  read.csv("./Mouse_Eimeria_Field/data_input/Mouse_data/HZ10_HZ17_Genotypes_47-523_Location.csv")
Gen18                     <-  read.csv("./Mouse_Eimeria_Field/data_input/Mouse_data/HZ18_Genotypes.csv")
Dis19                     <-  read.csv("./Mouse_Eimeria_Field/data_input/Mouse_data/HZ19_Dissections.csv")

f_Gen16 <- Gen16 %>% select(Mouse_ID, Transect, Latitude, Longitude, Year, HI, Sex, HI_NLoci) #, HI, HI_NLoci)
f_Gen17 <- Gen17 %>% select(Mouse_ID, Transect, Latitude, Longitude, Year, HI, Sex, HI_NLoci) #, HI, HI_NLoci)
f_Gen18 <- Gen18 %>% select(Mouse_ID, Transect, Latitude, Longitude, Year, HI, Sex, HI_NLoci) #, HI, HI_NLoci)
f_Dis19 <- Dis19 %>% select(Mouse_ID, Transect, Latitude, Longitude, Year, Sex)

EmanuelData <- EmanuelData %>%
  mutate(Mouse_ID  = PIN) %>%
  select(Mouse_ID, HI, HI_NLoci, Sex)

join   <- bind_rows(f_Gen16, f_Gen17, f_Gen18, f_Dis19)

full_join <- full_join(join, EmanuelData, by = "Mouse_ID") %>%
  distinct(Mouse_ID, .keep_all = T) %>%
  mutate(HI = HI.x,
         HI_NLoci = HI_NLoci.x) %>%
  select(Mouse_ID, Transect, Latitude, Longitude, Sex.x, HI, HI_NLoci)


full_join %>%
  vis_miss()

HI_Emanuel <- EmanuelData %>%
  select(Mouse_ID, HI)
HI_Rest <- full_join %>%
  select(Mouse_ID, HI)

