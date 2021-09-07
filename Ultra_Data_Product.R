library(tidyverse)
library(data.table)
library(visdat)

#### Select Columns ############################################################
basics          <- c("Mouse_ID", "Transect", "Code", "Region", "Sex", "Longitude", 
                     "Latitude", "Year", "State", "HI", "HI_NLoci")

Crypto_DNA.cols <- c("ILWE_DNA_Content_ng.microliter", "ILWE_DNA_used_up")

Crypto_qPCR.cols <- c("Tested_by", "Machine", "Measurements", "Ct_mean", "Ct_mean_Ep", 
                      "Ct_mean_ABI", "Oocyst_Predict")


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

parasite.cols   <- c("Aspiculuris", "Aspiculuris_tetraptera", "Syphacia_obvelata", "Aspiculuris_Syphacia",
                     "Trichuris", "Taenia", "Flea", "Mix_Syphacia_Aspiculuris", "Heligmosomoides_polygurus",
                     "Heterakis", "Mastophorus","Hymenolepis", "Taenia_martis", "Catenotaenia_pusilla",
                     "Hymenolepis_microstoma", "Hymenolepis_diminuta", 
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


Alice <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/MiceTableMusAliceArticle.csv")
Alice$HI_NLoci <- gsub(pattern = "HI ", replacement = "", x = Alice$HI_NLoci)
Alice$HI_NLoci <- as.integer(Alice$HI_NLoci)
Alice$Mouse_ID <- gsub(pattern = "Sk3173", replacement = "SK_3173", x = Alice$Mouse_ID)
wsh <- c(paste0("AA_000", 1:9), paste0("AA_00", 10:46))
apd <- c("A_0001", "A_0002", "A_0003")
useless <- c(wsh, apd)
Alice <- Alice[!(Alice$Mouse_ID %in% useless),]

Jarda <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/EmanuelData.csv", na.strings=c(""," ","NA"))

## change Column names to standardizes names for merging
setnames(Jarda, old = c("PIN", "X_Longit", "Y_Latit"), new = c("Mouse_ID", "Longitude", "Latitude"), skip_absent = T)
Jarda$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_", x = Jarda$Mouse_ID)
Jarda$Mouse_ID <- gsub(pattern = "Sk3173", replacement = "SK_3173", x = Jarda$Mouse_ID)



new_Alice <- left_join(Alice, Jarda) %>% select(!which(!rowSums(!is.na(Alice)))) %>% select(!which(!colSums(!is.na(Alice))))
new_Alice <- full_join(Alice, Jarda)

## Column Corrections of New_Alice #############################################


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
MT_Ectoparasites_Logical$Ectoparasites_Logical <- as.logical(MT_Ectoparasites_Logical$Ectoparasites_Logical)
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

## correct Year
new_Alice$Year[ new_Alice$Mouse_ID %in% c("SK_2903", "SK_2904")] <- 2014
new_Alice$Year[ new_Alice$Mouse_ID %in% c("AA_0330", "AA_0450", "AA_0451", "AA_0452")] <- 2017


## eliminate duplicates in new_Alice
new_Alice$Mouse_ID[duplicated(new_Alice$Mouse_ID)]
new_Alice <- new_Alice %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 



## Status
new_Alice <- new_Alice %>% mutate(Status = replace(Status, Status == "(pregnant)", "pregnant"),
                                  Status = replace(Status, Status == "(young)", "young"))
## State
new_Alice <- new_Alice %>% mutate(State = replace(State, State == "Germany", "DE"),
                                  State = replace(State, State == "D", "DE"),
                                  State = replace(State, State == "Poland", "PL"))
## multiple Mice per Box
new_Alice <- new_Alice %>% mutate(Multiple_Mice_per_Box = ifelse(Mouse_ID %in% c("AA_0514", "AA_0515", "AA_0513","AA_0349", "AA_0454", "ZZ_0037", "ZZ_0038"), TRUE, FALSE))

###############################################################################
new_Alice <- new_Alice[colnames(new_Alice) %in% c(basics, dissection.cols, EimGeno.cols, EqPCR.cols, gen.loci, oocyst.cols, parasite.cols)]
new_Alice <- new_Alice %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 
#new_Alice$Mouse_ID[duplicated(new_Alice$Mouse_ID)]

## add 2018 Oocyst count data
Eflot2018 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Eimeria_detection/HZ18_Eim_Flotation.csv")
Eflot2018$Ncells <- Eflot2018$Sume
Eflot2018$PBS_dil_in_mL <- Eflot2018$PBS_vol
Eflot2018$Feces_Weight <- Eflot2018$Feces

colnames(Eflot2018)[colnames(Eflot2018)%in%oocyst.cols]
Eflot2018 <- Eflot2018[colnames(Eflot2018)%in%c(basics,oocyst.cols)]

## add 2018 qPCR data
EqPCR2018 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Eimeria_detection/HZ18_qPCR.csv")
colnames(EqPCR2018)[colnames(EqPCR2018)%in%"delta"] <- "delta_ct_cewe_MminusE"
EqPCR2018 <- EqPCR2018[colnames(EqPCR2018)%in%c(basics, EqPCR.cols)]

## add 2018 Eim Genotyping data
EimPCR <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Eimeria_detection/Svenja/table_ct_and_more.csv")
EimPCR$Mouse_ID <- gsub("CEWE_AA_", "AA_0", EimPCR$Name)
EimPCR$eimeriaSpecies <-  gsub("E\\. ", "E_", EimPCR$Eimeria.subspecies)
EimPCR$eimeriaSpecies[EimPCR$eimeriaSpecies%in%c("non infected", "Eimeria sp.")] <- "Negative"
EimPCR$Sex <- NULL
EimPCR <- EimPCR[colnames(EimPCR)%in%c(basics, EimGeno.cols)]

Detection18 <- merge(EimPCR, EqPCR2018)
Detection18 <- merge(Detection18, Eflot2018)

new_Alice <- full_join(new_Alice, Detection18)


EqPCR2019 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Eimeria_detection/HZ19_CEWE_qPCR.csv")
colnames(EqPCR2019)[colnames(EqPCR2019)%in%"delta"] <- "delta_ct_cewe_MminusE"
colnames(EqPCR2019)[colnames(EqPCR2019)%in%"MC"] <- "MC.Eimeria"
EqPCR2019 <- EqPCR2019[colnames(EqPCR2019)%in%c(basics, EqPCR.cols)]

new_Alice <- full_join(new_Alice, EqPCR2019)



## adjust colnames  ###########################################################
## Address #####################################################################
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
MT_Ectoparasites_Logical$Ectoparasites_Logical <- as.logical(MT_Ectoparasites_Logical$Ectoparasites_Logical)
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

## Year == "Year", "year"
MT_Year <- new_Alice %>% select(Mouse_ID, Year) %>% pivot_longer(names_to = "Temp",  values_to = "Year", cols = c(Year)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Year)  %>% distinct(Mouse_ID, .keep_all = T) 
MT_Year$Mouse_ID[duplicated(MT_Year$Mouse_ID)]    
## join
new_Alice <- full_join(new_Alice, MT_Year)
rm(MT_Year)


## adjust col number ##########################################################

## add 2018 dissection data
Dis2018 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ18_Dissections.csv")
Dis2018$Aspiculuris           <- Dis2018$ASP
Dis2018$Syphacia_obvelata     <- Dis2018$SYP
Dis2018$Heterakis             <- Dis2018$HET
Dis2018$Taenia_martis         <- Dis2018$MART
Dis2018$Catenotaenia_pusilla  <- Dis2018$CP
Dis2018$Hymenolepis_microstoma <- Dis2018$HM
Dis2018$Hymenolepis_diminuta  <- Dis2018$HD
Dis2018$Trichuris_muris       <- Dis2018$TM
Dis2018$Mastophorus           <- Dis2018$MM
Dis2018$Left_Epididymis       <- Dis2018$Epididymis
Dis2018$Body_Weight           <- Dis2018$Body_weight
Dis2018$Body_Length           <- Dis2018$Body_length
Dis2018$Tail_Length           <- Dis2018$Tail_length
Dis2018 <- Dis2018[colnames(Dis2018) %in% c("Mouse_ID", basics, dissection.cols, parasite.cols)]




## add 2019 dissection data
Dis2019 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ19_Dissections.csv")
Dis2019$Aspiculuris           <- Dis2019$ASP
Dis2019$Syphacia_obvelata     <- Dis2019$SYP
Dis2019$Heterakis             <- Dis2019$HET
Dis2019$Taenia_martis         <- Dis2019$MART
Dis2019$Catenotaenia_pusilla  <- Dis2019$CP
Dis2019$Hymenolepis_microstoma <- Dis2019$HM
Dis2019$Hymenolepis_diminuta  <- Dis2019$HD
Dis2019$Trichuris_muris       <- Dis2019$TM
Dis2019$Mastophorus           <- Dis2019$MM
Dis2019$Left_Epididymis       <- Dis2019$Epididymis
Dis2019$Body_Weight           <- Dis2019$Body_weight
Dis2019$Body_Length           <- Dis2019$Body_length
Dis2019$Tail_Length           <- Dis2019$Tail_length
Dis2019 <- Dis2019[colnames(Dis2019) %in% c("Mouse_ID", basics, dissection.cols, parasite.cols)]


## merge
new_Alice <- left_join(new_Alice, Dis2018)
new_Alice <- full_join(new_Alice, Dis2019)

## find and remove duplicates
new_Alice <- new_Alice %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## selected passages needed
new_Alice <- new_Alice[colnames(new_Alice) %in% c(basics, dissection.cols, EimGeno.cols, EqPCR.cols, gen.loci, oocyst.cols, parasite.cols)]



## add:
# Gene Expression
# CEWE Elisa
# MES FACS
# Immuno
# Crypto
# other Rodents?

