library(dplyr)
library(tidyverse)
library(visdat)
library(stringr)
library(data.table)
library(reshape)


###############################################################################
    
## Load new_Alice Table with selected columns as a basis to add new data onto:
new_Alice   <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/new_Alice_MiceTable_Selected_Cols.csv")
## make sure that you toss out the Count Column called "X"
new_Alice   <- new_Alice[!names(new_Alice) == "X"]




## 1) add new mice from Jarda (contains new Mouse_IDs and HI)
    new_Jarda <- read.csv("yourLinkHere.csv")
    ## is Mouse_ID == Mouse_ID?, otherwise e.g.:
    ## setnames(Jarda, old = c("PIN), new = c("Mouse_ID"))
    ## join with full_join()
    new_Alice <- full_join(new_Alice, new_Jarda)
    ## find duplicates, fill rows, remove duplicates: 
    new_Alice <- new_Alice %>% 
      arrange(Mouse_ID) %>% 
      group_by(Mouse_ID) %>% 
      fill(c(everything()), .direction = "downup") %>% 
      ungroup() %>% 
      distinct(Mouse_ID, .keep_all = T) 

    
    
## 2) add new Dissection Data
    new_Dis <- read.csv("yourLinkHere.csv")
    ## is Mouse_ID == Mouse_ID?, otherwise use setnames(df, old = c(), new = c("Mouse_ID"))
    ## Join with full_join()
    new_Alice <- full_join(new_Alice, new_Dis)
    ## find duplicates, fill rows, remove duplicates
    new_Alice <- new_Alice %>%  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  distinct(Mouse_ID, .keep_all = T) 

    
    
## 3) add new Genotyping Data
    new_Gen <- read.csv("yourLinkHere.csv")
    ## is Mouse_ID == Mouse_ID?, otherwise use setnames(df, old = c(), new = c("Mouse_ID"))
    ## join with full_join()
    new_Alice <- full_join(new_Alice, new_Gen)
    ## find duplicates, fill rows, remove duplicates
    new_Alice <- new_Alice %>%  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  distinct(Mouse_ID, .keep_all = T) 
    


## 4) add new qPCR, Eimeria, Crypto data
    any_new_Data <- read.csv("yourLinkHere.csv")
    ## join with full_join()
    new_Alice <- full_join(new_Alice, any_new_Data)
    ## find duplicates, fill rows, remove duplicates
    new_Alice <- new_Alice %>%  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  distinct(Mouse_ID, .keep_all = T) 
    


## 5) Column Correction
    ## if there are Columns, that have a different name in the new data, correct them
    ## Colnames should technically be:
    colnames(new_Alice)
    ## if you have discrepancies, e.g. Latitude vs. latitude 
    ## and now you have 2 columns, which are missing different values
    ## you can combine them by applying:
    Col_to_adjust <- new_Alice %>% select(Mouse_ID, Wrong_Name, Correct_Name) %>%pivot_longer(names_to = "Temp", values_to = "Correct_Name", cols = c(Wrong_Name, Correct_Name)) %>% 
      arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>%  
      select(Mouse_ID, Correct_Name)
    ## now you have joined the columns and filled the gaps, time to remove duplicates
    Col_to_adjust <- Col_to_adjust %>% distinct(Mouse_ID, .keep_all = T) 
    ## join, deselect (remove) the wrong column name from popping up later
    new_Alice <- full_join(new_Alice, Col_to_adjust) %>% select(-c(Wrong_Name))
    rm(Col_to_adjust)
    
## Hopefully there aren't many Columns you need to correct, otherwise this might take a minute. 
    ## If you have ensured that the column names are now joined correctly
    ## and there is no new data you can add, you can look at your baby!
    ## Let's visualize the missingness with the vis_miss() function
    vis_miss(new_Alice, cluster = T, sort_miss = T)
    ## are some row entries missing, that can be adjusted manually?
    ## e.g. Location with Lat/Long is given, but there is no Address --> look up and add
    ## e.g. State is sometimes "Germany", "D", or "DEU" --> standardize
    new_Alice <- new_Alice %>% mutate(State = replace(State, State == "Germany", "DE"),
                                      State = replace(State, State == "D", "DE"),
                                      State = replace(State, State == "Poland", "PL"))
    ## sometimes there are Columns that have been filled incorrectly, 
    ## because the type was not taken into account, e.g. Flea & Fleas
        ## was used as numeric (0-12 Fleas were counted)
        ## or as a logical Condition (TRUE, FALSE)
        ## or as a general finding (fleas)
    ## so this should be corrected by either introducing separate Columns 
    ## or adjusting the input as needed, here I separated the input as two columns
    ## which take into account the logical or numeric input with an ifelse() fct:
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

    

################################################################################
## Interesting Column Names we want to stick to: ###############################
    basics          <- c("Mouse_ID", "Transect", "Code", "Region", "Sex", "Longitude", "Latitude", "Year", "State", "HI", "HI_NLoci")
    
    ILWE_DNA.cols   <- c("ILWE_DNA_Content_ng.microliter", "ILWE_DNA_used_up")
    
    Crypto.cols     <- c("Tested_by", "Machine", "Measurements", "Ct_mean", "Ct_mean_Ep", "Ct_mean_ABI", "Oocyst_Predict")
    
    gen.loci        <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
                         "Es1", "Gpd1", "Idh1", "Mpi", "Np", "Sod1", "Es1C", "Gpd1C",
                         "Idh1C", "MpiC", "NpC", "Sod1C", "HI_NLoci",
                         "HI", "Zfy2", "SRY1", "Y")
    
    dissection.cols <- c("Body_Weight", "Body_Length", "Tail_Length", "Spleen", 
                         "Left_Testis", "Right_Testis", "Seminal_Vesicles_Weight", 
                         "Sperm", "Left_Epididymis", 
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
    
################################################################################
    
    ## save as csv file
    write.csv(new_Alice, "full_Mouse_Data_2021.csv")
    
    
    
