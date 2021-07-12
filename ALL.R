library(dplyr)
library(tidyverse)
library(visdat)
library(ggplot2)
library(stringr)
library(data.table)


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
join   <- bind_rows(f_Gen16, f_Gen17, f_Gen18, f_Dis19)

EmanuelData <- EmanuelData %>%
  mutate(Mouse_ID  = PIN) %>%
  select(Mouse_ID, HI, HI_NLoci, Sex)


full_join <- full_join(join, EmanuelData, by = "Mouse_ID") %>%
  distinct(Mouse_ID, .keep_all = T)
# joined missing HI, HI_NLoci, Sex data manually
write.csv(full_join, "full_join.csv")
full_join <- read.csv("./Cryptosporidium-BSc/full_join_Corrected.csv") %>%
  select(Mouse_ID, Transect, Latitude, Longitude, Sex, HI.Gen.Jarda, HI_NLoci.Gen.Jarda)

full_join %>% vis_miss()

full_join %>% count(HI.Gen.Jarda)         # 71 unconfirmed, 207 NA / 2186 obs.
full_join %>% count(HI_NLoci.Gen.Jarda)   # 75 unconfirmed, 321 NA /  2186 obs.
# unconfirmed == two different HI or HI_NLoci values in Genotype vs. Jarda data
#### Delta_max_HI = 0.14      range[0,1]
#### Delta_max_HI_NLoci = 6   range[0,14]
# NA == no data available: 
#### neither table had information or Delta could not be calculated due to single value supplied










