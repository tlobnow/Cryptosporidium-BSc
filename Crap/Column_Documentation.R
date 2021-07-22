library(dplyr)
library(tidyverse)
library(visdat)
library(ggplot2)
library(stringr)
library(data.table)
library(reshape)
library(ggmap)
library(stringdist)

## Selected Columns
basics <- c("Mouse_ID", "Transect", "Code", "Region", "Sex", "Longitude", "Latitude", "Year")

gen.loci <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
              "Es1", "Gpd1", "Idh1", "Mpi", "Np", "Sod1", "Es1C", "Gpd1C",
              "Idh1C", "MpiC", "NpC", "Sod1C", "HI_NLoci",
              "HI", "Zfy2", "SRY1", "Y")

dissection.cols <- c("Body_weight", "Body_length", "Tail_length", "Spleen", 
                     "Left_testis", "Right_testis", "SemVes", 
                     "Sperm",
                     "Arrival", "Wean", "Death", "Dissection_date", "DaysInLab",
                     "Testes","Lepid", "Seminal_vesicles", "Uterus", "Ovaria", "Grav", "Litters", "NN", 
                     "MM", "FF", "Protocol")

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



