library(tidyverse)
library(data.table)
library(visdat)
library(viridis)
library(colorRamps)

#### Select Columns ############################################################
basics          <- c("Mouse_ID", "Address", "Sex", "Longitude", "Latitude", "Year", "HI", "HI_NLoci")
gen.loci        <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1", "Es1", "Gpd1", "Idh1", "Mpi", "Np", "Sod1", "Es1C", "Gpd1C", "Idh1C", "MpiC", "NpC", "Sod1C", "HI_NLoci", "HI", "Zfy2", "Y", "Region")
dissection.cols <- c("Body_Weight", "Body_Length", "Tail_Length", "Status", "Spleen", 
                     "Left_Testis", "Right_Testis", "Seminal_Vesicles_Weight", "Liver",
                     "Sperm", "Left_Epididymis", "Right_Epididymis", 
                     "Right_Ovarium_Weight", "Left_Ovarium_Weight",
                     "Left_Embryo", "Right_Embryo", "Fleas", "Ticks", "Ectoparasites_Logical", 
                     "Arrival", "Dissection_Date", "Trap_Date", "Host")
tissue.cols     <- c("SPL1", "SPL2", "ELFO", "LIV", "KI", "LUN", "SG",
                     "MES", "COWE", "COCE", "COCE2", "CEWE", "CECE", 
                     "ILWE", "SICE", "WEOH", "WFOR", "FEC")
initial.worms.cols   <- c("Aspiculuris","Syphacia_obvelata","Trichuris_muris", "Taenia_taeniformis", "Flea", "Ticks", "Fleas",
                          "Ectoparasites",    "Worms_presence", "Syphacia",
                          "Hymenolepis_diminiuta", "Hymenolepis_diminuta", "Taenia_martis",    "Heligmosomoides_polygurus" ,
                          "Heterakis_spumosa","Mastophorus_muris","Hymenolepis_microstoma", "Catenotaenia_pusilla",
                          "Cysticercus", "Hymenolepis", "Taenia", "Aspiculuris_Syphacia",  
                          "Trichuris", "Heterakis", "Mastophorus", "Ectoparasites_Logical", 
                          "Aspiculuris", "Catenotaenia_pusilla")
final.worms.cols <- c("Aspiculuris_sp", "Syphacia_sp", "Trichuris_muris", "Taenia_sp",
                      "Heterakis_sp", "Mastophorus_muris", "Hymenolepis_sp", "Catenotaenia_pusilla",
                      "Heligmosomoides_polygurus", "Worms_presence")
oocyst.cols     <- c("counter", "Feces_Weight", "Date_count", "N_oocysts_sq1",
                     "N_oocysts_sq2", "N_oocysts_sq3",  "N_oocysts_sq4",
                     "N_oocysts_sq5", "N_oocysts_sq6", "N_oocysts_sq7",
                     "N_oocysts_sq8", "mean_neubauer", "PBS_dil_in_mL", 
                     "OPG", "Ncells")
EqPCR.cols      <- c("delta_ct_ilwe_MminusE", "delta_ct_cewe_MminusE", "MC.Eimeria", "Ct.Eimeria", "Ct.Mus")
EimGeno.cols    <- c("n18S_Seq", "COI_Seq", "ORF470_Seq", "eimeriaSpecies")
Gene.Exp.cols   <- c("IFNy", "IL.12", "IRG6", "CXCR3", "IL.6", "GBP2",
                     "IL.10", "IL.13", "IL.10", "IL.13", "IL1RN",
                     "CXCR3", "CASP1", "CXCL9", "GAPDH", 
                     "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", 
                     "NCR1", "PPIB", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")
CellCount.cols <- c( "Treg", "CD4", "Treg17", "Th1", "Th17", "CD8",
                     "Act_CD8", "IFNy_CD4", "IL17A_CD4", "IFNy_CD8")
Crypto_qPCR.cols <- c("Ct_mean", "Oocyst_Predict")

SOTA <- read.csv("SOTA_Data_Product.csv") %>% select(-X)

a <- SOTA[colnames(SOTA)%in% c(basics, Crypto_qPCR.cols)] %>% filter(Oocyst_Predict > 5000) %>% filter(Mouse_ID != "AA_0322")

write.csv(a, "High_Infection_Samples.csv")

