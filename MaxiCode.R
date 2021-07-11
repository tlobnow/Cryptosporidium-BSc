library(ggplot2)
library(dplyr)
library(ggeffects)
library(assertive)
library(stringr)
library(visdat)
library(stringdist)
library(fuzzyjoin)
library(reclin)
library(RColorBrewer)
library(viridis)
library(viridisLite)
# display.brewer.all()
# library(DHARMa) # investigation of models
# library(AICcmodavg) # which model has the best AIC score (complexity vs. accuracy)


ABI_thSC_log2   <- read.csv("ABI_thSC_log2.csv")
ABI_Best_thSC   <- read.csv("ABI_Best_SC.csv")
ABI_thSC_log10  <- read.csv("ABI_thSC_log10.csv")
Eppi_thSC_log10 <- read.csv("Eppi_thSC.csv")
ABI_SC_log2     <- read.csv("ABI_SC_log2.csv") %>%
  select(X, Sample, Ct_mean, Assay, Oocyst_mean, Oocyst_Prediction)
ABI_SC_log10    <- read.csv("ABI_SC_log10.csv") %>%
  select(X, Sample, Ct_mean, Assay, Oocyst_mean, OP_ABI_log10)
Eppi_SC_log10   <- read.csv("Eppi_SC.csv") %>%
  select(Sample, Ct_mean, Oocyst_mean)

Eppi_Samples    <- read.csv("Eppi_Samples.csv") %>%
  select(Sample, Ct_mean, Oocyst_mean)
#HI_Data   <- read.csv("Hybrid_Data.csv") %>%
#  select(PIN, Year, Transect, Sex, X_Longit, Y_Latit, HI) %>%
#  arrange(PIN)
full1     <- nf_full1 %>% filter(Ct_mean > 0)
nf_full1  <- read.csv("nf_full1_Corrected.csv")
f_ABI_Best_thSC <- read.csv("f_ABI_Best_thSC_predicted.csv")
f_ABI_thSC_log2 <- read.csv("f_ABI_thSC_log2_predicted.csv")
f_2016          <- read.csv("f_2016.csv")
f_2017          <- read.csv("f_2017.csv")
f_2018          <- read.csv("f_2018.csv")
f_2019          <- read.csv("f_2019.csv")
nf_2016         <- read.csv("nf_2016.csv")
nf_2017         <- read.csv("nf_2017.csv")
nf_2018         <- read.csv("nf_2018.csv")
nf_2019         <- read.csv("nf_2019.csv")
full_Gen        <- read.csv("full_Gen.csv")
fd_ABI_log2     <- read.csv("fd_ABI_log2.csv")


f_ABI_thSC_log2

   setwd("~/Documents/GitHub/Mouse_Eimeria_Field/data/Field_data/")
Gen16 <-  read.csv("HZ16_Genotypes_47-211.csv")
Gen17 <-  read.csv("HZ10_HZ17_Genotypes_47-523_Location.csv")
Gen18 <-  read.csv("HZ18_Genotypes.csv")
Dis19 <-  read.csv("HZ19_Dissections.csv")
   setwd("~/Documents/Programming/R/HZ_SC_and_Raw_Data")


# filter data frames
f_ABI_thSC_log2   <-  filter(ABI_thSC_log2, Ct_mean > 0, Assay != 5)
f_ABI_thSC_log10  <-  filter(ABI_thSC_log10, Ct_mean > 0, Ct_mean != 36.42)
f_ABI_Best_thSC   <-  filter(ABI_Best_thSC, Ct_mean > 0)
f_ABI_SC_log2     <-  filter(ABI_SC_log2, Ct_mean > 0, Sample != 250000, Assay != 5)
f_ABI_SC_log10    <-  filter(ABI_SC_log10, Ct_mean > 0)
f_Eppi_thSC_log10 <-  filter(Eppi_thSC_log10, Ct_mean > 0)
f_Eppi_SC_log10   <-  filter(Eppi_SC_log10, Ct_mean >0)
f_Eppi_Samples    <-  filter(Eppi_Samples, Ct_mean > 0)

linear_model0 <- lm(log2(Amount_Oocysts) ~ Ct_mean, data = f_ABI_Best_thSC)
linear_model1 <- lm(log2(Amount_Oocysts) ~ Ct_mean, data = f_ABI_thSC_log2)
linear_model2 <- lm(log10(Amount_Oocysts) ~ Ct_mean, data = f_ABI_thSC_log10)
linear_model3 <- lm(log10(Amount_Oocysts) ~ Ct_mean, data = f_Eppi_thSC_log10)


f_ABI_Best_thSC$predicted     <- 2^predict(linear_model0)
f_ABI_thSC_log2$predicted     <- 2^predict(linear_model1)
f_ABI_thSC_log10$predicted    <- 10^predict(linear_model2)
f_Eppi_thSC_log10$predicted   <- 10^predict(linear_model3)

Oocyst_Prediction         <-  2^predict(linear_model0, newdata = full1)
OP_ABI_log10              <-  10^predict(linear_model2, newdata = full1)
OP_Eppi_log10             <-  10^predict(linear_model3, newdata = full1)
OP_2016                   <-  2^predict(linear_model0, newdata = f_2016)
OP_2017                   <-  2^predict(linear_model0, newdata = f_2017)
OP_2018                   <-  2^predict(linear_model0, newdata = f_2018)
OP_2019                   <-  2^predict(linear_model0, newdata = f_2019)
Oocyst_Predict            <-  2^predict(linear_model0, newdata = full1)

nf_Oocyst_Prediction      <-  2^predict(linear_model0, newdata = nf_full1)
nf_OP_ABI_log10           <- 10^predict(linear_model2, newdata = nf_full1)
nf_OP_Eppi_log10          <- 10^predict(linear_model3, newdata = nf_full1)
nf_OP_2016                <-  2^predict(linear_model0, newdata = nf_2016)
nf_OP_2017                <-  2^predict(linear_model0, newdata = nf_2017)
nf_OP_2018                <-  2^predict(linear_model0, newdata = nf_2018)
nf_OP_2019                <-  2^predict(linear_model0, newdata = nf_2019)
fd_ABI_log2               <- data.frame(full1, Oocyst_Prediction)
fd_ABI_log10              <- data.frame(full1, OP_ABI_log10)
fd_Eppi_log10             <- data.frame(full1, OP_Eppi_log10)
nf_fd_ABI_log2            <- data.frame(nf_full1, nf_Oocyst_Prediction)
nf_fd_ABI_log10           <- data.frame(nf_full1, nf_OP_ABI_log10)
nf_fd_Eppi_log10          <- data.frame(nf_full1, nf_OP_Eppi_log10)
all_machines_comparison   <- data.frame(nf_full1, nf_Oocyst_Prediction, nf_OP_ABI_log10, nf_OP_Eppi_log10)
f_all_machines_comparison <- data.frame(full1, Oocyst_Prediction, OP_ABI_log10, OP_Eppi_log10)
#f_ABI_SC_log2            <- data.frame(f_ABI_SC_log2, Oocyst_Prediction)
#ABI_SC_log10             <- data.frame(f_ABI_SC_log10, OP_ABI_log10)
#Eppi_SC_log10            <- data.frame(f_Eppi_SC_log10, OP_Eppi_log10)




# csv files created throughout code
#   write.csv(fd_ABI_log2, "fd_ABI_log2.csv")     # fd_ABI_log2
#   write.csv(fd_ABI_log10, "fd_ABI_log10.csv")   # fd_ABI_log10
#   write.csv(fd_Eppi_log10, "fd_Eppi_log10.csv") # fd_Eppi_log10
#   write.csv(ABI_SC_log2, "ABI_SC_log2.csv")     # ABI_Test
#   write.csv(ABI_SC_log10, "ABI_SC_log10.csv")   # ABI_Test
#   write.csv(Eppi_SC_log10, "Eppi_SC_log10.csv") # ABI_Test
#   write.csv(full_Gen, "full_Gen.csv")           # HZ_Data_Combo
#   write.csv(all_Samples_Complete, "all_Samples_Complete.csv")             # HZ_Data_Combo
#   write.csv(f_all_machines_comparison, "f_all_machines_comparison.csv")   # Eppi_ABI_full_Data_Sheet.R, all pos. Samples
#   write.csv(all_machines_comparison, "all_machines_comparison.csv")       # Eppi_ABI_full_Data_Sheet.R, all Samples, pos and negative
#   write.csv(fd_OP_2016, "2016.csv")   # Eppi_ABI_Years_Comp. 
#   write.csv(fd_OP_2017, "2017.csv")   # Eppi_ABI_Years_Comp. 
#   write.csv(fd_OP_2018, "2018.csv")   # Eppi_ABI_Years_Comp. 
#   write.csv(fd_OP_2019, "2019.csv")   # Eppi_ABI_Years_Comp. 
#   write.csv(full1, "full1.csv")       # Eppi_ABI_Years_Comp. 
#   write.csv(nf_full1, "nf_full1.csv") # All Samples per year with HI
#   write.csv(f_2016, "f_2016.csv")
#   write.csv(f_2017, "f_2017.csv")
#   write.csv(f_2018, "f_2018.csv")
#   write.csv(f_2019, "f_2019.csv")
#   write.csv(f_HI_Data_19, "f_HI_Data_19.csv")
#   write.csv(nf_2016, "nf_2016.csv")
#   write.csv(nf_2017, "nf_2017.csv")
#   write.csv(nf_2018, "nf_2018.csv")
#   write.csv(nf_2019, "nf_2019.csv")
#   write.csv(nf_fd_OP_2016, "nf_fd_OP_2016.csv")
#   write.csv(nf_fd_OP_2017, "nf_fd_OP_2017.csv")
#   write.csv(nf_fd_OP_2018, "nf_fd_OP_2018.csv")
#   write.csv(nf_fd_OP_2019, "nf_fd_OP_2019.csv")
#write.csv(Gen16, "HZ16_Genotypes_47-211.csv")
#write.csv(Gen17, "HZ10_HZ17_Genotypes_47-523_Location.csv")
#write.csv(Gen18, "HZ18_Genotypes.csv")
#write.csv(Dis19, "HZ19_Dissections.csv")


# ABI_log_2.R ##################################################################
# changing it to linMod0 for annotation of the entire data!
linear_model0 <- lm(log2(Amount_Oocysts) ~ Ct_mean, data = f_ABI_Best_thSC)

# prediction 
f_ABI_Best_thSC$predicted   <- 2^predict(linear_model0)


# predict the model on samples (HZ19)
Oocyst_Prediction   <-  2^predict(linear_model0, newdata = full1)
Oocyst_Prediction

fd_ABI_log2   <- data.frame(full1, Oocyst_Prediction)
fd_ABI_log2 %>%
count(Transect)

nf_full1 %>%
  count(Transect)

ggplot(fd_ABI_log2, aes(Ct_mean, log2(Oocyst_Prediction))) +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  geom_smooth(method = "lm")

# plot predicted Oocysts versus Samples
fd_ABI_log2_Plot <- ggplot(fd_ABI_log2, aes(Mouse_ID, Oocyst_Prediction, col = Ct_mean)) +
  geom_point(stat = "summary", fun = "mean") +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  scale_color_gradient(low = "blue", high = " light blue") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust=1))
fd_ABI_log2_Plot

#write.csv(fd_ABI_log2, "fd_ABI_log2.csv")



# Eppi_ABI_full_Data_Sheet.R ##################################################

# DATA FRAMES 
fd_ABI_log2   <- data.frame(full1, Oocyst_Prediction)
fd_ABI_log10  <- data.frame(full1, OP_ABI_log10)
fd_Eppi_log10 <- data.frame(full1, OP_Eppi_log10)

nf_fd_ABI_log2 <- data.frame(nf_full1, nf_Oocyst_Prediction)
nf_fd_ABI_log10 <- data.frame(nf_full1, nf_OP_ABI_log10)
nf_fd_Eppi_log10 <- data.frame(nf_full1, nf_OP_Eppi_log10)

all_machines_comparison <- data.frame(nf_full1, nf_Oocyst_Prediction, nf_OP_ABI_log10, nf_OP_Eppi_log10)
all_machines_comparison

f_all_machines_comparison <- data.frame(full1, Oocyst_Prediction, OP_ABI_log10, OP_Eppi_log10)
f_all_machines_comparison

# PLOT THE DATA FRAMES 

# filtered samples, contains no negative samples
ggplot(fd_ABI_log2, aes(Ct_mean, log2(Oocyst_Prediction))) +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  geom_smooth(method = "lm", se = T)
ggplot(fd_ABI_log10, aes(Ct_mean, log10(OP_ABI_log10))) +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  geom_smooth(method = "lm", se = T)
ggplot(fd_Eppi_log10, aes(Ct_mean, log10(OP_Eppi_log10))) +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  geom_smooth(method = "lm", se = T)

ggp <- ggplot(NULL, aes(Ct_mean, log2(Oocyst_Prediction), color = Year)) +
  geom_point(data = fd_ABI_log2, size = 2, shape = 1) +
  geom_point(data = fd_ABI_log10, size = 2, shape = 1) +
  geom_point(data = fd_Eppi_log10, size = 2, shape = 1) +
  scale_color_gradient(low = "blue", high = "green")
ggp

# plot predicted Oocysts versus Samples 
fd_ABI_log2_Plot <- ggplot(fd_ABI_log2, aes(Mouse_ID, Oocyst_Prediction, col = Ct_mean)) +
  geom_point(stat = "summary", fun = "mean") +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  scale_color_gradient(low = "blue", high = "light blue") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust=1))
fd_ABI_log2_Plot

#write.csv(fd_ABI_log2, "fd_ABI_log2.csv")
#write.csv(fd_ABI_log10, "fd_ABI_log10.csv")
#write.csv(fd_Eppi_log10, "fd_Eppi_log10.csv")


# positive Samples filtered
#write.csv(f_all_machines_comparison, "f_all_machines_comparison.csv")
# all Samples, pos and negative
#write.csv(all_machines_comparison, "all_machines_comparison.csv")

# ABI_Test_Oocysts.R ###########################################################
## This file is for double-checking the efficiency of the Standard curve
## by feeding the linear model with known Standard curve values.
## prediction values should mimic the content for each dilution level of the SC
# plot the filtered data frames
ggplot(f_ABI_SC_log2, aes(log2(Oocyst_mean), Ct_mean)) +
  geom_point(stat = "summary", fun = "mean") +
  geom_smooth(method = "lm") +
  geom_point(position = position_jitter(0), size = 2, shape = 1)
ggplot(f_ABI_thSC_log2, aes(log2(Amount_Oocysts), Ct_mean)) +
  geom_point(stat = "summary", fun = "mean") +
  geom_smooth(method = "lm") +
  geom_point(position = position_jitter(0), size = 2, shape = 1)
ggplot(f_ABI_Best_thSC, aes(log2(Amount_Oocysts), Ct_mean)) +
  geom_point(stat = "summary", fun = "mean") +
  geom_smooth(method = "lm") +
  geom_point(position = position_jitter(0), size = 2, shape = 1)


# Create a linear model to predict the oocysts in unknown samples
linear_model0 <- lm(log2(Amount_Oocysts) ~ Ct_mean, data = f_ABI_Best_thSC)
linear_model1 <- lm(log2(Amount_Oocysts) ~ Ct_mean, data = f_ABI_thSC_log2)

# THIS IS IT!!!!
f_ABI_thSC_log2$predicted <- 2^predict(linear_model1)
f_ABI_thSC_log2
f_ABI_Best_thSC$predicted <- 2^predict(linear_model0)
f_ABI_Best_thSC

Oocyst_Prediction <- 2^predict(linear_model1, newdata = f_ABI_SC_log2)
OP_ABI_log10    <-  2^predict(linear_model0, newdata = f_ABI_SC_log10)
OP_Eppi_log10   <-  2^predict(linear_model0, newdata = f_Eppi_SC_log10)

# data frame
# op predict == predicted with normal thSC, op.predict.1 == predicted w/ best_thSC
f_ABI_SC_log2   <- data.frame(f_ABI_SC_log2, Oocyst_Prediction)
f_ABI_SC_log10  <- data.frame(f_ABI_SC_log10, OP_ABI_log10)
fd_Eppi_log10 <- data.frame(f_Eppi_SC_log10, OP_Eppi_log10)

summary(linear_model1)
summary(linear_model0)
ggpredict(linear_model1)
ggpredict(linear_model0)

#write.csv(f_ABI_thSC_log2, "f_ABI_thSC_log2_predicted.csv")
#write.csv(f_ABI_thSC_log2, "f_ABI_Best_thSC_predicted.csv")



# HZ_Data_Combo.R ##############################################################
# filtered Generations 
f_Gen16 <- Gen16 %>%
  select(Mouse_ID,
         Letter,
         Transect, 
         Latitude, 
         Longitude,
         Year,
         Sex,
         HI, 
         HI_NLoci)
f_Gen16

glimpse(f_Gen16)
glimpse(f_Gen17)

f_Gen17 <- Gen17 %>%
  select(Mouse_ID,
         Transect, 
         Latitude, 
         Longitude,
         Year,
         Sex,
         HI,
         HI_NLoci)

f_Gen17

f_Gen18 <- Gen18 %>%
  select(Mouse_ID, 
         Transect, 
         Latitude, 
         Longitude,
         Year,
         Sex,
         HI, 
         HI_NLoci)
f_Gen18

f_Dis19 <- Dis19 %>%
  select(Mouse_ID, 
         Transect, 
         Latitude, 
         Longitude,
         Year,
         Sex)
f_Dis19

# JOINING DATA FRAMES 
join1 <- full_join(f_Gen16, f_Gen17)
join2 <- full_join(f_Gen18, f_Dis19)
full_Gen <- full_join(join1, join2) %>%
  select(Mouse_ID, 
         Transect, 
         Latitude, 
         Longitude,
         Year,
         Sex,
         HI,
         HI_NLoci)

full_Gen_unique <- full_Gen %>%
  distinct(Mouse_ID, .keep_all = TRUE)
full_Gen_unique %>%
  count(Mouse_ID) %>%
  filter(n > 1)

Emanuel_Data <- read.csv("EmanuelData.csv") %>%
  mutate(Mouse_ID = PIN)
full_Data <- left_join(full_Gen_unique, Emanuel_Data, by = "Mouse_ID", "Transect")
full_Data <- left_join(full_Data, nf_full1, by = "Mouse_ID")
full_Data <- full_Data %>%
  select(Mouse_ID,
         Year.x, Date, Plate, Machine, Tested_by,
         Ct_1_Ep, Ct_2_Ep, Ct_Ep_Consistent,
         Ct_1_ABI, Ct_2_ABI, Ct_3_ABI,
         Ct_4_ABI, Ct_5_ABI, Ct_6_ABI, 
         Ct_7_ABI, Ct_8_ABI, Ct_9_ABI,
         Ct_mean_1_ABI, Ct_mean_2_ABI, Ct_mean_3_ABI,
         Ct_mean_Ep, Ct_mean_ABI, Ct_mean, 
         StDev, Ct_1_Ep_Real_Pred, Ct_2_Ep_Real_Pred, 
         Consistent, 
         HI.x, HI_Level, HI_NLoci.x, Sex.x, 
         Transect.x, Latitude.x, Longitude.x,
         Oocyst_Predict)

glimpse(full_Data)
glimpse(nf_full1)



full_Gen_HI <- read.csv("alles.csv")
#  select(Mouse_ID, Transect, Latitude, Longitude, Year, Sex, HI)
#setdiff(full_Gen_unique, full_Gen_HI)

full_Gen <- left_join(full_Gen_unique, full_Gen_HI) %>%
  select(Mouse_ID, Transect, Latitude, Longitude, Year, Sex, HI)

full_Data <- left_join(full_Gen_unique, full_Gen_HI)
full_Data <- left_join(full_Data, nf_full1)

write.csv(full_Gen, "full_Gen.csv")


# Eppi_ABI_Years_Comparison ####################################################

HI_Data    <- read.csv("Hybrid_Data.csv")

f_HI_Data_19 <- HI_Data %>%
  filter(Year == 2019) %>%
  select(PIN, Year, Transect, Sex, X_Longit, Y_Latit, HI) %>%
  arrange(PIN)

Mouse_ID  <- f_HI_Data_19$PIN
HI        <- f_HI_Data_19$HI

f_HI_Data_19 <- data.frame(Mouse_ID, HI)

# filter the data frames
f_ABI_Best_thSC <-  filter(ABI_Best_thSC, Ct_mean > 0, Assay != 5)
full1   <- nf_full1 %>%
  filter(Ct_mean > 0) %>%
  mutate(ID = row_number()) %>%
  select(ID, Mouse_ID, Year, Ct_mean, Transect, Latitude, Longitude,
         Sex, HI)

# all POSITIVE Samples per year with HI ########################################

f_2016  <- nf_full1 %>%
  filter(Year == 2016, Ct_mean > 0) %>%
  mutate(ID = row_number()) %>%
  select(ID, Mouse_ID, Year, Ct_mean, 
         Transect, Latitude, Longitude,
         Sex, HI)

f_2017  <- nf_full1 %>%
  filter(Year == 2017, Ct_mean > 0) %>%
  mutate(ID = row_number()) %>%
  select(ID, Mouse_ID, Year, Ct_mean, 
         Transect, Latitude, Longitude,
         Sex, HI)
f_2018  <- nf_full1 %>%
  filter(Year == 2018, Ct_mean > 0) %>%
  mutate(ID = row_number()) %>%
  select(ID, Mouse_ID, Year, Ct_mean, 
         Transect, Latitude, Longitude,
         Sex, HI)
f2019  <- nf_full1 %>%
  filter(Year == 2019) %>%
  mutate(ID = row_number()) %>%
  select(ID, Mouse_ID, Year, Ct_mean, 
         Transect, Latitude, Longitude,
         Sex, HI)
f2019

f_2019 <- left_join(f2019, f_HI_Data_19, by = "Mouse_ID")
f_2019




# PLOTS

# plot the filtered data frames
ggplot(f_ABI_Best_thSC, aes(log2(Amount_Oocysts), Ct_mean)) +
  geom_point(stat = "summary", fun = "mean") +
  geom_smooth(method = "lm") +
  geom_point(position = position_jitter(0), size = 2, shape = 1)
ggplot(full1, aes(Mouse_ID, Ct_mean, col = Ct_mean)) +
  geom_point(stat = "summary", fun = "mean") +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



# Create a linear model to predict the oocysts in unknown samples
linear_model0 <- lm(log2(Amount_Oocysts) ~ Ct_mean, data = f_ABI_Best_thSC)
summary(linear_model0)
ggpredict(linear_model0)

# prediction
f_ABI_Best_thSC$predicted   <- 2^predict(linear_model0)

# predict the model on samples (HZ19)

Oocyst_Prediction   <- 2^predict(linear_model0, newdata = full1)
OP_2016   <-  2^predict(linear_model0, newdata = f_2016)
OP_2016
OP_2017   <-  2^predict(linear_model0, newdata = f_2017)
OP_2017
OP_2018   <-  2^predict(linear_model0, newdata = f_2018)
OP_2018
OP_2019   <-  2^predict(linear_model0, newdata = f_2019)
OP_2019

# produce data frames
fd_ABI_log2 <- data.frame(full1, Oocyst_Prediction)
fd_OP_2016  <- data.frame(f_2016, OP_2016) %>%
  mutate(pos_ratio = 41 / 152) 
fd_OP_2017  <- data.frame(f_2017, OP_2017) %>%
  mutate(pos_ratio = 64 / 210) 
fd_OP_2018  <- data.frame(f_2018, OP_2018) %>%
  mutate(pos_ratio = 19 / 155) 
fd_OP_2019  <- data.frame(f_2019, OP_2019) %>%
  mutate(pos_ratio = 17 / 136) 



full <-  union_all(fd_OP_2016, fd_OP_2017)
full0 <- union_all(full, fd_OP_2018)
full1 <- union_all(full0, fd_OP_2019) 
Oocyst_Predict <- 2^predict(linear_model0, newdata = full1)

full1 <- data.frame(full1, Oocyst_Predict) %>%
  select(ID, Mouse_ID, Year, Ct_mean, HI, Transect, Latitude,
         Longitude, Sex, HI, Oocyst_Predict, pos_ratio)
full1

write.csv(full1, "full1.csv")


# PLOTS ----


ggplot(fd_ABI_log2, aes(Ct_mean, log2(Oocyst_Prediction))) +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  geom_smooth(method = "lm")

# plot predicted Oocysts versus Samples
fd_ABI_log2_Plot <- ggplot(fd_ABI_log2, aes(Mouse_ID, Oocyst_Prediction, col = Ct_mean)) +
  geom_point(stat = "summary", fun = "mean") +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  scale_color_gradient(low = "blue", high = " light blue") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust=1))
fd_ABI_log2_Plot


## 2016 Plot
#ggsave("fd16.jpg", width = 8, height = 5)
fd16_Plot <- full1 %>%
  filter(Year == 2016) %>%
  ggplot(aes(Mouse_ID, Oocyst_Predict, col = Ct_mean)) +
  ggtitle("2016 Predicted Oocysts in Positive Samples") +
  geom_point(stat = "summary", fun = "mean") +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  scale_color_gradient(low = "red", high = " yellow") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust=1))
fd16_Plot

## 2017 Plot
#ggsave("fd17.jpg", width = 8, height = 5)
fd17_Plot <- full1 %>%
  filter(Year == 2017) %>%
  ggplot(aes(Mouse_ID, Oocyst_Predict, col = Ct_mean)) +
  ggtitle("2017 Predicted Oocysts in Positive Samples") +
  geom_point(stat = "summary", fun = "mean") +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  scale_color_gradient(low = "red", high = " yellow") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust=1))
fd17_Plot

## 2018 Plot
#ggsave("fd18.jpg", width = 8, height = 5)
fd18_Plot <- full1 %>%
  filter(Year == 2018) %>%
  ggplot(aes(Mouse_ID, Oocyst_Predict, col = Ct_mean)) +
  ggtitle("2018 Predicted Oocysts in Positive Samples") +
  geom_point(stat = "summary", fun = "mean") +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  scale_color_gradient(low = "red", high = " yellow") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust=1))
fd18_Plot

## 2019 Plot
#ggsave("fd19.jpg", width = 8, height = 5)
fd19_Plot <- full1 %>%
  filter(Year == 2019) %>%
  ggplot(aes(Mouse_ID, Oocyst_Predict, col = Ct_mean)) +
  ggtitle("2019 Predicted Oocysts in Positive Samples") +
  geom_point(stat = "summary", fun = "mean") +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  scale_color_gradient(low = "blue", high = " red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust=1))
fd19_Plot



# Comparison of positive cases per year -----

# Normalizing Positives / Samples per Year
# 2016: 41 pos  / 153 Samples
# 2017: 64 pos  / 220 Samples
# 2018: 19 pos  / 179 Samples
# 2019: 17 pos  / 144 Samples

#ggsave("SampleRatio.jpg", width = 8, height = 5)
ggplot(full1, aes(x = Year, y = ratio, fill = Year)) +
  ggtitle("Crypto Positive Samples per Year") +
  geom_bar(stat = "summary", fun = "mean")


# SAVE CSV FILES


write.csv(fd_OP_2016, "2016.csv")
write.csv(fd_OP_2017, "2017.csv")
write.csv(fd_OP_2018, "2018.csv")
write.csv(fd_OP_2019, "2019.csv")

# introducing nf_full1 ----

# all Samples per Year, non filtered
nf_2016  <- nf_full1 %>%
  filter(Year == 2016) %>%
  mutate(ID = row_number()) %>%
  select(ID, Mouse_ID, Year, Ct_mean,
         Transect, Latitude, Longitude,
         Sex, HI)

nf_2017  <- nf_full1 %>%
  filter(Year == 2017) %>%
  mutate(ID = row_number()) %>%
  select(ID, Mouse_ID, Year, Ct_mean,
         Transect, Latitude, Longitude,
         Sex, HI)
nf_2018  <- nf_full1 %>%
  filter(Year == 2018) %>%
  mutate(ID = row_number()) %>%
  select(ID, Mouse_ID, Year, Ct_mean,
         Transect, Latitude, Longitude,
         Sex, HI)
nf2019  <- nf_full1 %>%
  filter(Year == 2019) %>%
  mutate(ID = row_number()) %>%
  select(ID, Mouse_ID, Year, Ct_mean,
         Transect, Latitude, Longitude,
         Sex)
nf2019

nf_2019 <- full_join(nf2019, f_HI_Data_19, by = "Mouse_ID")
nf_2019


Oocyst_Prediction   <- 2^predict(linear_model0, newdata = full1)
nf_OP_2016   <-  2^predict(linear_model0, newdata = nf_2016)
nf_OP_2016
nf_OP_2017   <-  2^predict(linear_model0, newdata = nf_2017)
nf_OP_2017
nf_OP_2018   <-  2^predict(linear_model0, newdata = nf_2018)
nf_OP_2018
nf_OP_2019   <-  2^predict(linear_model0, newdata = nf_2019)
nf_OP_2019

# produce data frames

fd_ABI_log2 <- data.frame(full1, Oocyst_Prediction)
nf_fd_OP_2016  <- data.frame(nf_2016, nf_OP_2016) %>%
  mutate(pos_ratio = 41 / 152)
nf_fd_OP_2017  <- data.frame(nf_2017, nf_OP_2017) %>%
  mutate(pos_ratio = 64 / 210)
nf_fd_OP_2018  <- data.frame(nf_2018, nf_OP_2018) %>%
  mutate(pos_ratio = 19 / 155) 
nf_fd_OP_2019  <- data.frame(nf_2019, nf_OP_2019) %>%
  mutate(pos_ratio = 17 / 136)

 nf_full <-  union_all(nf_fd_OP_2016, nf_fd_OP_2017)
 nf_full0 <- union_all(nf_full, nf_fd_OP_2018)
 nf_full1 <- union_all(nf_full0, nf_fd_OP_2019) 
 Oocyst_Predict <- 2^predict(linear_model0, newdata = nf_full1)



 nf_full1 <- data.frame(nf_full1, Oocyst_Predict) %>%
  select(ID, Mouse_ID, Year, Ct_mean, HI, Transect, Latitude,
         Longitude, Sex, HI, Oocyst_Predict, pos_ratio)
 
 
 # CLEANING DATA ##############################################################
 nf_full1 <- nf_full1 %>%
   filter(!is.na(Year)) %>%
   mutate(Oocyst_Predict = replace(Oocyst_Predict, Oocyst_Predict == '4292821751815.77', '0')) %>%
   arrange(Mouse_ID)
 
 write.csv(nf_full1, "nf_full1_Corrected.csv")

 nf_full1
 
 nf_full1 %>%
   vis_miss()
 
 
 # CSV
 #write.csv(nf_fd_OP_2016, "nf_2016.csv")
 #write.csv(nf_fd_OP_2017, "nf_2017.csv")
 #write.csv(nf_fd_OP_2018, "nf_2018.csv")
 #write.csv(nf_fd_OP_2019, "nf_2019.csv")
 #write.csv(Combo_unique, "nf_full1.csv")
 #write.csv(nf_full1, "nf_full1.csv")
 #nf_full1 <- read.csv("nf_full1.csv")

 


# GGPLOT GRAPHS ###############################################################

# theoretical SC
ggplot(f_ABI_thSC_log2, aes(log2(Amount_Oocysts), Ct_mean)) +
  geom_point(stat = "summary", fun = "mean") +
  geom_smooth(method = "lm") +
  geom_point(position = position_jitter(0), size = 2, shape = 1)
# best SC - ggsave("best_thSC_log2.jpg", width = 8, height = 5)
ggplot(f_ABI_Best_thSC, aes(log2(Amount_Oocysts), Ct_mean)) +
  ggtitle("Standard Curve 1:8 dilution") +
  geom_point(stat = "summary", fun = "mean") +
  geom_smooth(method = "lm") +
  geom_point(position = position_jitter(0), size = 2, shape = 1)
#  all Samples
ggplot(full1, aes(Mouse_ID, Ct_mean, col = Ct_mean)) +
  geom_point(stat = "summary", fun = "mean") +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Comparison of Oocyst predictions across the years
ggp <- ggplot(NULL, aes(Ct_mean, log2(Oocyst_Prediction), color = Year)) +
  geom_point(data = fd_ABI_log2, size = 2, shape = 1) +
  geom_point(data = fd_ABI_log10, size = 2, shape = 1) +
  geom_point(data = fd_Eppi_log10, size = 2, shape = 1) +
  scale_color_gradient(low = "yellow", high = "red")
ggp

# unfiltered samples, contains negative samples!
ggplot(nf_fd_ABI_log2, aes(Ct_mean, log2(nf_Oocyst_Prediction))) +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  geom_smooth(method = "lm", se = T)
ggplot(nf_fd_ABI_log10, aes(Ct_mean, log10(nf_OP_ABI_log10))) +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  geom_smooth(method = "lm", se = T)
ggplot(nf_fd_Eppi_log10, aes(Ct_mean, log10(nf_OP_Eppi_log10))) +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  geom_smooth(method = "lm", se = T)

## 2016 Plot - ggsave("fd16.jpg", width = 8, height = 5)
fd16_Plot <- full1 %>%
  filter(Year == 2016) %>%
  ggplot(aes(Mouse_ID, Oocyst_Predict, col = Ct_mean)) +
  ggtitle("2016 Predicted Oocysts in Positive Samples") +
  geom_point(stat = "summary", fun = "mean") +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  scale_color_gradient(low = "red", high = " yellow") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust=1))
fd16_Plot

## 2017 Plot - ggsave("fd17.jpg", width = 8, height = 5)
fd17_Plot <- full1 %>%
  filter(Year == 2017) %>%
  ggplot(aes(Mouse_ID, Oocyst_Predict, col = Ct_mean)) +
  ggtitle("2017 Predicted Oocysts in Positive Samples") +
  geom_point(stat = "summary", fun = "mean") +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  scale_color_gradient(low = "red", high = " yellow") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust=1))
fd17_Plot

## 2018 Plot - ggsave("fd18.jpg", width = 8, height = 5)
fd18_Plot <- full1 %>%
  filter(Year == 2018) %>%
  ggplot(aes(Mouse_ID, Oocyst_Predict, col = Ct_mean)) +
  ggtitle("2018 Predicted Oocysts in Positive Samples") +
  geom_point(stat = "summary", fun = "mean") +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  scale_color_gradient(low = "red", high = " yellow") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust=1))
fd18_Plot

## 2019 Plot - ggsave("fd19.jpg", width = 8, height = 5)
fd19_Plot <- full1 %>%
  filter(Year == 2019) %>%
  ggplot(aes(Mouse_ID, Oocyst_Predict, col = Ct_mean)) +
  ggtitle("2019 Predicted Oocysts in Positive Samples") +
  geom_point(stat = "summary", fun = "mean") +
  geom_point(position = position_jitter(0), size = 2, shape = 1) +
  scale_color_gradient(low = "blue", high = " red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust=1))
fd19_Plot

# Comparison of positive cases per year
##ggsave("SampleRatio.jpg", width = 8, height = 5)
ggplot(full1, aes(x = Year, y = ratio_pos_per_yr, fill = ratio_pos_per_yr)) +
  ggtitle("Crypto Positive Samples per Year") +
  geom_bar(stat = "summary", fun = "mean")

# mice caught per year
ggplot(full1, aes(x = Year, y = ratio_mus_caught, fill = ratio_mus_caught)) +
  ggtitle("Mus caught per Year") +
  geom_bar(stat = "summary", fun = "mean")






