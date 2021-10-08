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


SOTA <- SOTA %>%
  mutate(Eim_Species = ifelse(eimeriaSpecies == "E_falciformis", "E_falciformis",
                              ifelse(eimeriaSpecies == "E_ferrisi", "E_ferrisi",
                                     ifelse(eimeriaSpecies == "Eimeria_alorani", "Eimeria_alorani",
                                            ifelse(eimeriaSpecies == "Eimeria_apionodes", "Eimeria_apionodes",
                                                   ifelse(eimeriaSpecies == "Eimeria_falciformis", "Eimeria_falciformis",
                                                          ifelse(eimeriaSpecies == "Eimeria_sp_Apodemus", "Eimeria_sp_Apodemus",
                                                                 ifelse(eimeriaSpecies == "Eimeria_vermiformis", "Eimeria_vermiformis",
                                                                        ifelse(eimeriaSpecies == "Negative", "Negative",
                                                                               ifelse(NA)))))))))) %>%
  mutate(Tested_for = case_when(Ct_mean >= 0 & is.na(Ct.Eimeria) & is.na(Eim_Species) ~ "Crypto",
                                is.na(Ct_mean) & !is.na(Eim_Species) ~ "Eimeria",
                                is.na(Ct_mean) & (!is.na(Ct.Eimeria) | !is.na(Eim_Species)) ~ "Eimeria",
                                !is.na(Ct_mean) & (!is.na(Ct.Eimeria) | !is.na(Eim_Species)) ~ "Crypto_Eimeria",
                                is.na(Ct_mean) & is.na(Ct.Eimeria) & is.na(Eim_Species) ~ "none",
                                is.na(Ct_mean) & (is.na(Ct.Eimeria) | is.na(Eim_Species)) ~ "none"),
         Crypto_Positive  = case_when(Ct_mean > 0 ~ T,
                                      Ct_mean == 0 ~ F),
         Eimeria_Positive = case_when(Eim_Species != "Negative" | Ct.Eimeria > 0 | OPG > 0 ~ T,
                                      Eim_Species == "Negative" | Ct.Eimeria == 0 | OPG == 0 ~ F),
         Infection = case_when(Tested_for == "Crypto" & Ct_mean > 0                                                          ~ "Positive Crypto",
                               Tested_for == "Crypto" & Ct_mean == 0                                                         ~ "Negative Crypto",
                               Tested_for == "Eimeria" & Eimeria_Positive == T                                               ~ "Positive Eimeria",
                               Tested_for == "Eimeria" & Eimeria_Positive == F                                               ~ "Negative Eimeria",
                               Tested_for == "Crypto_Eimeria" & Ct_mean > 0  & (Ct.Eimeria > 0 |  Eim_Species != "Negative" | OPG > 0) ~ "Positive Crypto_Eimeria",
                               Tested_for == "Crypto_Eimeria" & Ct_mean > 0  & (Ct.Eimeria == 0 | Eim_Species == "Negative"| OPG == 0) ~ "Positive Crypto",
                               Tested_for == "Crypto_Eimeria" & Ct_mean == 0 & (Ct.Eimeria > 0 | Eim_Species != "Negative" | OPG > 0)  ~ "Positive Eimeria",
                               Tested_for == "Crypto_Eimeria" & Ct_mean == 0 & (Ct.Eimeria == 0 | Eim_Species == "Negative" | OPG == 0) ~ "Negative Crypto_Eimeria",
                               Tested_for == "none" ~ "not tested"),
         Infected = case_when(Infection == "not tested" ~ F,
                              Infection == "Negative Crypto" ~ F,
                              Crypto_Positive == T | Eimeria_Positive == T ~ T,
                              Crypto_Positive == F & Eimeria_Positive == F ~ F,
                              is.na(Crypto_Positive) & Eimeria_Positive == F ~ F,
                              Crypto_Positive == F & is.na(Eimeria_Positive) ~ F,
                              ))
                            

b <- SOTA %>% select(Ct_mean, Ct.Eimeria, eimeriaSpecies, Infection, Tested_for)
SOTA %>% select(Ct_mean, Ct.Eimeria, eimeriaSpecies, Infection, Tested_for) %>% count(Infection, Tested_for)

SOTA %>% count(Eimeria_Positive, Crypto_Positive, Infected)

## visualize HI vs. Value, faceted by Marker, colored by Infection (C, E, C+E)
a <- SOTA [colnames(SOTA) %in% c(basics, Gene.Exp.cols, EimGeno.cols, Crypto_qPCR.cols, "Tested_for", "Infection", "Eim_Species", "Infected")]
a <- a %>% pivot_longer(cols = c(10:34, -"IFNy"),
               names_to = "Marker",
               values_to = "Value") %>%
  filter(!is.na(Value))


a %>%
  ggplot(aes(Value, HI, col = Infection)) +
  geom_point() +
  facet_wrap(~ Marker)


a %>%
  filter(Infected == T) %>% 
  ggplot(aes(HI, Value, col = Marker)) +
  geom_point() +
  facet_wrap(~ Infection)


#### investigating Infection profiles
Pos_Crypto <- a %>% filter(Infection == "Positive Crypto")
Pos_Crypto %>%
  ggplot(aes(Value, HI, col = Value)) +
  geom_point() +
  facet_wrap(~ Marker) +
  scale_color_viridis(option = "inferno") +
  labs(title = "Crypto Positive Samples")

Pos_Eimeria <- a %>% filter(Infection == "Positive Eimeria")
Pos_Eimeria %>%
  ggplot(aes(Value, HI, col = Value)) +
  geom_point() +
  facet_wrap(~ Marker) +
  scale_color_viridis(option = "inferno") +
  labs(title = "Eimeria Positive Samples")

Pos_Crypto_Eimeria <- a %>% filter(Infection == "Positive Crypto_Eimeria")
Pos_Crypto_Eimeria %>%
  ggplot(aes(Value, HI, col = Value)) +
  geom_point(stat = "identity") +
  facet_wrap(~ Marker) +
  scale_color_viridis(option = "inferno") +
  labs(title = "Crypto and Eimeria Positive Samples")



