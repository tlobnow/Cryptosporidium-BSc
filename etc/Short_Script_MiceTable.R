library(dplyr)
library(tidyverse)
library(visdat)


###############################################################################
    
## Load MiceTable Table with selected columns as a basis to add new data onto:
MiceTable   <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/MiceTable_fullEimeriaInfos_2014to2017.csv")
MiceTable   <- MiceTable[colnames(MiceTable) %in% c(basics, "HI")]
  ## get rid of useless Samples 
  wsh <- c(paste0("AA_000", 1:9), paste0("AA_00", 10:46))
  apd <- c("A_0001", "A_0002", "A_0003")
  useless <- c(wsh, apd)
  MiceTable <- MiceTable[!(MiceTable$Mouse_ID %in% useless),]



## 1) add new mice from Jarda (contains new Mouse_IDs and HI)
    new_Jarda <- read.csv("yourLinkHere.csv")
    ## is Mouse_ID == Mouse_ID?, otherwise e.g.:
    ## setnames(Jarda, old = c("PIN), new = c("Mouse_ID"))
    ## join with full_join()
    MiceTable <- full_join(MiceTable, new_Jarda)
    ## find duplicates, fill rows, remove duplicates: 
    MiceTable <- MiceTable %>% 
      arrange(Mouse_ID) %>% 
      group_by(Mouse_ID) %>% 
      fill(c(everything()), .direction = "downup") %>% 
      ungroup() %>% 
      distinct(Mouse_ID, .keep_all = T) 

    
    
## 2) add new Dissection Data
    new_Dis <- read.csv("yourLinkHere.csv")
    ## is Mouse_ID == Mouse_ID?, otherwise use setnames(df, old = c(), new = c("Mouse_ID"))
    ## Join with full_join()
    MiceTable <- full_join(MiceTable, new_Dis)
    ## find duplicates, fill rows, remove duplicates
    MiceTable <- MiceTable %>%  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  distinct(Mouse_ID, .keep_all = T) 

    
    
## 3) add new Genotyping Data
    new_Gen <- read.csv("yourLinkHere.csv")
    ## is Mouse_ID == Mouse_ID?, otherwise use setnames(df, old = c(), new = c("Mouse_ID"))
    ## join with full_join()
    MiceTable <- full_join(MiceTable, new_Gen)
    ## find duplicates, fill rows, remove duplicates
    MiceTable <- MiceTable %>%  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  distinct(Mouse_ID, .keep_all = T) 
    


## 4) add new qPCR, Eimeria, Crypto data
    any_new_Data <- read.csv("yourLinkHere.csv")
    ## join with full_join()
    MiceTable <- full_join(MiceTable, any_new_Data)
    ## find duplicates, fill rows, remove duplicates
    MiceTable <- MiceTable %>%  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  distinct(Mouse_ID, .keep_all = T) 
    


## 5) Column Correction
    ## if there are Columns, that have a different name in the new data, correct them
    ## Colnames should technically be:
    colnames(MiceTable)
    ## if you have discrepancies, e.g. Latitude vs. latitude 
    ## and now you have 2 columns, which are missing different values
    ## you can combine them by applying:
    Col_to_adjust <- MiceTable %>% select(Mouse_ID, Wrong_Name, Correct_Name) %>%pivot_longer(names_to = "Temp", values_to = "Correct_Name", cols = c(Wrong_Name, Correct_Name)) %>% 
      arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>%  
      select(Mouse_ID, Correct_Name)
    ## now you have joined the columns and filled the gaps, time to remove duplicates
    Col_to_adjust <- Col_to_adjust %>% distinct(Mouse_ID, .keep_all = T) 
    ## join, deselect (remove) the wrong column name from popping up later
    MiceTable <- full_join(MiceTable, Col_to_adjust) %>% select(-c(Wrong_Name))
    rm(Col_to_adjust)
    
## Hopefully there aren't many Columns you need to correct, otherwise this might take a minute. 
    ## If you have ensured that the column names are now joined correctly
    ## and there is no new data you can add, you can look at your baby!
    ## Let's visualize the missingness with the vis_miss() function
    vis_miss(MiceTable, cluster = T, sort_miss = T)
 

    

################################################################################
## Interesting Column Names we want to stick to: ###############################
    basics            <- c("Mouse_ID", "Transect", "Code", "Region", "Sex", "Longitude", "Latitude", "Year", "State", "HI", "HI_NLoci")
    
   
    gen.loci          <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
                         "Es1", "Gpd1", "Idh1", "Mpi", "Np", "Sod1", "Es1C", "Gpd1C",
                         "Idh1C", "MpiC", "NpC", "Sod1C", "HI_NLoci",
                         "HI", "Zfy2", "SRY1", "Y")
    
    dissection.cols   <- c("Body_Weight", "Body_Length", "Tail_Length", "Spleen", 
                         "Left_Testis", "Right_Testis", "Seminal_Vesicles_Weight", 
                         "Sperm", "Left_Epididymis", 
                         "Arrival", "Wean", "Death", "Dissection_date", "DaysInLab",
                         "Lepid", "Uterus", "Ovaria", "Grav", "Litters", "NN", 
                         "MM", "FF", "Protocol")
    
    parasite.cols     <- c("Aspiculuris_tetraptera", "Syphacia_obvelata", "Aspiculuris_Syphacia",
                         "Trichuris", "Taenia", "Flea", "Mix_Syphacia_Aspiculuris", "Heligmosomoides_polygurus",
                         "Heterakis","Heterakis", "Mastophorus","Hymenolepis",
                         "Ectoparasites", "Ectoparasites_Count", "Ectoparasites_Logical",
                         "Worms_presence")
    
    oocyst.cols       <- c("counter", "Feces_Weight", "Date_count", "N_oocysts_sq1",
                         "N_oocysts_sq2", "N_oocysts_sq3",  "N_oocysts_sq4",
                         "N_oocysts_sq5", "N_oocysts_sq6", "N_oocysts_sq7",
                         "N_oocysts_sq8", "mean_neubauer", "PBS_dil_in_mL", 
                         "OPG", "Ncells")
    
    EqPCR.cols        <- c("delta_ct_ilwe_MminusE", "delta_ct_cewe_MminusE",
                         ## from 2018 on AN IMPORTANT IMPROVEMENT!!!
                         "MC.Eimeria")
    
    EimGeno.cols      <- c("n18S_Seq", "COI_Seq", "ORF470_Seq", "eimeriaSpecies")
    
    Crypto_DNA.cols   <- c("ILWE_DNA_Content_ng.microliter", "ILWE_DNA_used_up", 
                         "DNA_Dilution_ng.microliter", "dest_Water")
    
    Crypto_qPCR.cols  <- c("Ct_mean", "Ct_mean_Ep", "Ct_mean_ABI", 
                         "Ct_mean_1_ABI", "Ct_mean_2_ABI", "Ct_mean_3_ABI",
                         "Ct_believable", "Machine", "Measurements", "Tested_by", 
                         "qPCR_Date")
    
################################################################################
    
    ## save as csv file
    write.csv(MiceTable, "full_Mouse_Data_2021.csv")
    
    
    
