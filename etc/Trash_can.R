
ALL <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/MiceTable_fullEimeriaInfos_2014to2017.csv")


ALL   <- ALL[,colnames(ALL) %in% c(basics, gen.loci, "HI"), ]
vis_miss(ALL)
# 69.5% present




#### add 2018 Data #############################################################
### The basics from the 2018 dissection
Dis2018 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ18_Dissections.csv")

### The mouse genotyping from Jarda's table
Gen2018 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ18_Genotypes.csv")

## So Jarda's genotyping table acutally has all the main data!  As
## long as we don't want to look at the other parasites we have to do
## nothing to this before the merge
DisGen2018 <- Gen2018[, colnames(Gen2018)%in%c(basics, "HI"), ] %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


#### add 2019 Data #############################################################
### The basics from the 2019 dissection
Dis2019 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ19_Dissections.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)

### extract 2019 Genotyping from Jarda's table
Gen_2010_2019   <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ10-19_Genotypes.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE) %>% filter(Year == "2019")
Gen_2010_2019$EH_ID <- gsub(pattern = "SK", replacement = "SK_", x = Gen_2010_2019$EH_ID)
setnames(Gen_2010_2019, old = c("EH_ID"), new = c("Mouse_ID"), skip_absent = T)

### Join Dis
DisGen2019 = merge(DisGen2019, Jarda[, colnames(Jarda) %in% c(basics,"HI"), ])


### JOINING ###################################################################

### Join 2018 and 2019 Data 
DisGen_18_19          <- full_join(DisGen2018[colnames(DisGen2018)%in%c(basics, "HI")], DisGen2019[colnames(DisGen2019)%in%c(basics, gen.loci)]) %>%
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 



## includes FULLY available data, all gen.loci, all qPCR data
#Crypto  <- full_join(Crypto[colnames(Crypto) %in% c("Mouse_ID", qPCR.cols, ILWE_DNA.cols)], Jarda[colnames(Jarda) %in% c(basics, "HI", gen.loci)]) %>% filter(Ct_mean >= 0)
## would include data that is missing gen.loci data + HI
#vis_miss(Crypto, cluster = T)



parasite.cols <- c("Aspiculuris_tetraptera", "Syphacia_obvelata", "Trichuris_muris",
                   "Taenia_taeniformis", "Flea", "Mix_Syphacia_Aspiculuris",
                   "Heterakis_spumosa", "Mastophorus_muris", "Hymenolepis_microstoma",
                   "Catenotaenia_pusilla", "Cysticercus", "Ectoparasites",
                   "Worms_presence", "Hymenolepis_diminiuta", "Taenia_martis",
                   "Heligmosomoides_polygurus", "Taenia", "Aspiculuris_Syphacia",
                   "Trichuris", "Heterakis", "Mastophorus")
