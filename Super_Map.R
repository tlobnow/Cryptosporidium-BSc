# Load Libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer) #display.brewer.all()
library(viridis)
library(colorRamps)
library(cowplot)
library(uwot)
library(leaflet)
library(leaflet.extras)
library(sp)
library(htmltools)
library(stringr)
library(patchwork)

# Load Data ####

full_Data <- read.csv("full_Data_Corrected.csv")
full_Data <- full_Data %>%
  mutate(HI_State = ifelse(HI == 0.000, "Mmd", 
                           ifelse(HI == 1.000, "Mmm",
                                  ifelse(HI != 0.000 & HI <= 0.500, "Eastern Hybrid",
                                         ifelse(HI != 1.000 & HI > 0.500, "Western Hybrid",
                                                ifelse(NA)))))) %>%
  mutate(Oocyst_Level = ifelse(Oocyst_Predict == 0, "OP = 0",
                               ifelse(Oocyst_Predict <= 100, "OP < 100",
                                      ifelse(Oocyst_Predict <= 10^3, "OP < 1000",
                                             ifelse(Oocyst_Predict <= 10^4, "OP < 10.000",
                                                    ifelse(Oocyst_Predict <= 10^5, "OP < 100.000",
                                                           ifelse(Oocyst_Predict <= 10^6, "OP < 1.000.000")
                                                    ))))))
full_Data$Location <- paste(full_Data$Latitude, full_Data$Longitude, sep = ",")
full_Data_long_pull <- full_Data %>%
  count(Longitude)%>% mutate(mus_caught = n) %>% select(Longitude, mus_caught) %>% arrange(mus_caught)

full_Data_long_pull_Crypto <- full_Data %>%
  filter(Ct_mean > 0) %>% count(Longitude) %>% mutate(Crypto_mus_caught = n) %>% select(Longitude, Crypto_mus_caught) %>% arrange(Crypto_mus_caught)

full_Data1 <- left_join(full_Data, full_Data_long_pull)
full_Data2 <- left_join(full_Data1, full_Data_long_pull_Crypto)
full_Data <- full_Data2
rm(full_Data1)
rm(full_Data2)

glimpse(full_Data)

# Subset Data for work ####
sub_full_Data <- full_Data %>%
  select(Mouse_ID, Transect, Sex, 
         Latitude, Longitude, Location, 
         HI, HI_NLoci, HI_State,
         Year, Date, Plate, Tested_by, Machine, Measurements,
         Ct_believable, Ct_Ep_Consistent,
         Ct_1_Ep, Ct_2_Ep, Ct_1_Ep_Real_Pred, Ct_2_Ep_Real_Pred, Ct_Ep_Consistent,
         Ct_1_ABI, Ct_2_ABI, Ct_3_ABI, Ct_4_ABI, Ct_5_ABI, Ct_6_ABI, Ct_7_ABI, Ct_8_ABI, Ct_9_ABI,
         Ct_mean_1_ABI, Ct_mean_2_ABI, Ct_mean_3_ABI, 
         Ct_mean_Ep, Ct_mean_ABI, Ct_mean, 
         StDev, Oocyst_Predict, Oocyst_Level,
         mus_caught, Crypto_mus_caught
         )

# visualize #### 

#sub_full_Data %>%
#count(Ct_1_Ep_Real_Pred, Ct_2_Ep_Real_Pred)

#   Ct_1_Ep_Real_Pred        Ct_2_Ep_Real_Pred    n
#1             FALSE             FALSE            82
#2             FALSE              TRUE            100
#3              TRUE             FALSE            1
#4              TRUE              TRUE           413
#5                NA                NA           1486

#sub_full_Data %>%
#  filter(Ct_believable == T) %>%
#  ggplot(aes(Mouse_ID, Oocyst_Predict, col = Oocyst_Predict, label = Mouse_ID)) +
#  geom_point() +
#  scale_color_gradient(low = "blue", high = " red") +
#  theme(axis.text.x = element_blank())


# Map Data ####################################################################
map <- map <- full_Data %>%
  leaflet() %>%
  addProviderTiles("CartoDB") %>%
  setView(lat = 52.520007, lng =13.404954, zoom = 7)
map

glimpse(sub_full_Data)

# create Objects and Layers for later ####
# Transects ####
sub_full_Data %>%
  count(Transect)
sub_full_Data_Transect <- sub_full_Data %>%
  filter(Transect  != "NA",
         Longitude != "NA",
         Latitude  != "NA")
data_col_Transect = colorFactor(matlab.like2(6), sub_full_Data_Transect$Transect)
      Transect_Allo_PL <- sub_full_Data %>%
        filter(Transect == "Allo_PL")
      Transect_HZ_BAV <- sub_full_Data %>%
        filter(Transect == "HZ_BAV")
      Transect_HZ_BR <- sub_full_Data %>%
        filter(Transect == "HZ_BR")
      Transect_HZ_CZ <- sub_full_Data %>%
        filter(Transect == "HZ_CZ")
      Transect_HZ_MV <- sub_full_Data %>%
        filter(Transect == "HZ_MV")
      Transect_HZ_SX <- sub_full_Data %>%
        filter(Transect == "HZ_SX")

# Sex ####
sub_full_Data %>%
  count(Sex)
sub_full_Data_Sex <- sub_full_Data %>%
  filter(Sex  != "NA",
         Longitude != "NA",
         Latitude  != "NA")
data_col_Sex = colorFactor(matlab.like2(2), sub_full_Data_Sex$Sex)
     Sex_F <- sub_full_Data %>%
       filter(Sex == "F")
     Sex_M <- sub_full_Data %>%
       filter(Sex == "M")

# HI ####
sub_full_Data %>%
  count(HI)
sub_full_Data_HI <- sub_full_Data %>%
  filter(HI  != "NA",
         Longitude != "NA",
         Latitude  != "NA")
data_col_HI = colorFactor(matlab.like2(6), sub_full_Data_HI$HI)
     HI_0 <- sub_full_Data %>%
       filter(HI == 0)
     HI_bl_0.25 <- sub_full_Data %>%
       filter(HI > 0, HI <= 0.25)
     HI_bl_0.5 <- sub_full_Data %>%
       filter(HI > 0.25, HI <= 0.5)
     HI_bl_0.75 <- sub_full_Data %>%
       filter(HI > 0.5, HI <= 0.75)
     HI_bl_1 <- sub_full_Data %>%
       filter(HI > 0.75, HI < 1)
     HI_1 <- sub_full_Data %>%
       filter(HI == 1)

# HI State ####
sub_full_Data %>%
  count(HI_State)
sub_full_Data$HI_State <- cut(sub_full_Data$HI, c(0, 0.001, 0.5, 0.999, 1), include.lowest = T ,
                              labels = c('pure Mmd', 'Western Hybrid', 'Eastern Hybrid', 'pure Mmm'))
sub_full_Data_HI_State <- sub_full_Data %>%
  filter(HI_State  != "NA",
         Longitude != "NA",
         Latitude  != "NA")
data_col_HI_State = colorFactor(matlab.like2(4), sub_full_Data_HI_State$HI_State)


# Years ####
sub_full_Data %>%
  count(Year)
sub_full_Data_Year <- sub_full_Data %>%
  filter(Year  != "NA",
         Longitude != "NA",
         Latitude  != "NA")
data_col_Year = colorFactor(matlab.like2(9), sub_full_Data_Year$Year)
     Year_2010 <- sub_full_Data %>%
       filter(Year == "2010")
     Year_2011 <- sub_full_Data %>%
       filter(Year == "2011")
     Year_2012 <- sub_full_Data %>%
       filter(Year == "2012")
     Year_2013 <- sub_full_Data %>%
       filter(Year == "2013")
     Year_2014 <- sub_full_Data %>%
       filter(Year == "2014")
     Year_2015 <- sub_full_Data %>%
       filter(Year == "2015")
     Year_2016 <- sub_full_Data %>%
      filter(Year == "2016")
     Year_2017 <- sub_full_Data %>%
      filter(Year == "2017")
     Year_2018 <- sub_full_Data %>%
      filter(Year == "2018")
     Year_2019 <- sub_full_Data %>%
      filter(Year == "2019")

# Oocyst Predictions ####
sub_full_Data %>%
  count(Oocyst_Predict)
sub_full_Data$Oocyst_Level <- cut(sub_full_Data$Oocyst_Predict, c(1, 10^2, 10^3, 10^4, 10^5, 10^6), include.lowest = T,
                                  labels = c('< 100 Oocysts', '< 1000 Oocysts', '< 10.000 Oocysts', '< 100.000 Oocysts', '< 1 Mio. Oocysts'))
     
sub_full_Data_OP_Level <- sub_full_Data %>%
  filter(Oocyst_Level  != "NA",
         Oocyst_Level != "0",
         Longitude != "NA",
         Latitude  != "NA")
data_col_OP_Level = colorFactor(matlab.like2(5), sub_full_Data_OP_Level$Oocyst_Level)

# Tested by ####
sub_full_Data %>%
  count(Tested_by)
sub_full_Data_Tested_by <- sub_full_Data %>%
  filter(Tested_by  != "NA",
         Longitude != "NA",
         Latitude  != "NA")
data_col_Tested_by = colorFactor(matlab.like2(4), sub_full_Data_Tested_by$Tested_by)
    Tes <- sub_full_Data %>%
      filter(Tested_by == "Tes")
    Tes_Tes <- sub_full_Data %>%
      filter(Tested_by == "Tes-Tes")
    YVT <- sub_full_Data %>%
      filter(Tested_by == "YVT")
    YVT_Tes <- sub_full_Data %>%
      filter(Tested_by == "YVT-Tes")


# Machines ####
sub_full_Data %>%
  count(Machine)
sub_full_Data %>%
  count(Machine)
sub_full_Data_Machine <- sub_full_Data %>%
  filter(Machine  != "NA",
         Longitude != "NA",
         Latitude  != "NA")
data_col_Machine = colorFactor(matlab.like2(5), sub_full_Data_Machine$Machine)
    Ep <- sub_full_Data %>%
      filter(Machine == "Eppendorf")
    Ep_nABI <- sub_full_Data %>%
      filter(Machine == "Eppendorf + new ABI")
    nABI_nABI <- sub_full_Data %>%
      filter(Machine == "new ABI + new ABI")
    oABI_nABI <- sub_full_Data %>%
      filter(Machine == "old ABI + new ABI")
    noABI_nABI_nABI <- sub_full_Data %>%
      filter(Machine == "old ABI + new ABI + new ABI")

# Measurements  ####
sub_full_Data %>%
  count(Measurements)
sub_full_Data_Measurements <- sub_full_Data %>%
  filter(Measurements  != "NA",
         Longitude != "NA",
         Latitude  != "NA")
data_col_Measurements = colorFactor(matlab.like2(3), sub_full_Data_Measurements$Measurements)
    Measured_1x <- sub_full_Data %>%
      filter(Measurements == 1)
    Measured_2x <- sub_full_Data %>%
      filter(Measurements == 2)
    Measured_3x <- sub_full_Data %>%
      filter(Measurements == 3)


# Ct_believable  ####
sub_full_Data_Ct_bel <- sub_full_Data %>%
  filter(Ct_believable  != "NA",
         Longitude != "NA",
         Latitude  != "NA")
data_col_Ct_bel = colorFactor(matlab.like2(2), sub_full_Data_Ct_bel$Ct_believable)
    Ct_believable_T <- sub_full_Data %>%
      filter(Ct_believable == T)
    Ct_believable_F <- sub_full_Data %>%
      filter(Ct_believable == F)

# Ct mean ####
sub_full_Data %>%
  count(Ct_mean)
sub_full_Data$Ct_Level <- cut(sub_full_Data$Ct_mean, c(0, 1, 25, 30, 40), include.lowest = F ,
                              labels = c('0', '< 25', '< 30', '< 40'))
    
sub_full_Data_Ct_mean <- sub_full_Data %>%
  filter(Ct_mean  != "NA",
         Longitude != "NA",
         Latitude  != "NA")
data_col_Ct_mean = colorFactor(matlab.like2(4), sub_full_Data_Ct_mean$Ct_Level)
    Ct_bl_25 <- sub_full_Data %>%
      filter(Ct_mean < 25 & Ct_mean != 0)
    Ct_bl_30 <- sub_full_Data %>%
      filter(Ct_mean <30 & Ct_mean >= 25)
    Ct_bl_40 <- sub_full_Data %>%
      filter(Ct_mean <40 & Ct_mean >= 30)
    Ct_mean_Ep <- sub_full_Data %>%
      filter(Ct_mean_Ep != "NA")
    Ct_mean_ABI <- sub_full_Data %>%
      filter(Ct_mean_ABI != "NA")

# mus_caught ####
sub_full_Data %>%
  summarize(Longitude) %>%
  count(Longitude) %>%
  arrange(n)
sub_full_Data_mus_caught <- sub_full_Data %>%
  filter(mus_caught  != "NA",
         Longitude != "NA",
         Latitude  != "NA")
data_col_mus_caught = colorFactor(matlab.like2(5), sub_full_Data_mus_caught$mus_caught)
    mus_caught_10 <- sub_full_Data %>%
      filter(mus_caught <= 10)
    mus_caught_20 <- sub_full_Data %>%
      filter(mus_caught <= 20 & mus_caught > 10)
    mus_caught_30 <- sub_full_Data %>%
      filter(mus_caught <= 30 & mus_caught > 20)
    mus_caught_40 <- sub_full_Data %>%
      filter(mus_caught <= 40 & mus_caught > 30)
    mus_caught_50 <- sub_full_Data %>%
      filter(mus_caught <= 50 & mus_caught > 40)
    mus_caught_60 <- sub_full_Data %>%
      filter(mus_caught <= 60 & mus_caught > 50)

# Crypto mus caught ####
sub_full_Data %>%
  count(Crypto_mus_caught)
sub_full_Data_Crypto_caught <- sub_full_Data %>%
  filter(Crypto_mus_caught  != "NA",
         Longitude != "NA",
         Latitude  != "NA")
data_col_Crypto_caught = colorFactor(matlab.like2(4), sub_full_Data_Crypto_caught$Crypto_mus_caught)


# Map ####
map <- map <- sub_full_Data %>%
  leaflet() %>%
  addProviderTiles("CartoDB") %>%
  setView(lat = 52.520007, lng =13.404954, zoom = 7)
map


map %>%
  addCircleMarkers(data = sub_full_Data_Transect, 
                   col = ~data_col_Transect(Transect),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>HI State:<b>",as.character(HI_State), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Transect") %>%
  addCircleMarkers(data = sub_full_Data_Sex, 
                   col = ~data_col_Sex(Sex),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>HI State:<b>",as.character(HI_State), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Sex") %>%
  addCircleMarkers(data = sub_full_Data_HI_State, 
                   col = ~data_col_HI_State(HI_State),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>HI State:<b>",as.character(HI_State), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "HI_State") %>%
  addCircleMarkers(data = sub_full_Data_Year, 
                   col = ~data_col_Year(Year),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>HI State:<b>",as.character(HI_State), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Year") %>%
  addCircleMarkers(data = sub_full_Data_OP_Level, 
                   col = ~data_col_OP_Level(Oocyst_Level),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>HI State:<b>",as.character(HI_State), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Oocyst_Level",
                   opacity = 1) %>%
  addCircleMarkers(data = sub_full_Data_Tested_by, 
                   col = ~data_col_Tested_by(Tested_by),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>HI State:<b>",as.character(HI_State), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Tested_by",
                   opacity = 1) %>%
  addCircleMarkers(data = sub_full_Data_Machine, 
                   col = ~data_col_Machine(Machine),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>HI State:<b>",as.character(HI_State), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Machine",
                   opacity = 1) %>%
  addCircleMarkers(data = sub_full_Data_Measurements, 
                   col = ~data_col_Measurements(Measurements),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>HI State:<b>",as.character(HI_State), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Measurements",
                   opacity = 1) %>%
  addCircleMarkers(data = sub_full_Data_Ct_bel, 
                   col = ~data_col_Ct_bel(Ct_believable),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>HI State:<b>",as.character(HI_State), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Ct_believable",
                   opacity = 1) %>%
  addCircleMarkers(data = sub_full_Data_Ct_mean, 
                   col = ~data_col_Ct_mean(Ct_Level),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>HI State:<b>",as.character(HI_State), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Ct_mean",
                   opacity = 0.5) %>%
  addCircleMarkers(data = sub_full_Data_mus_caught, 
                   col = ~data_col_mus_caught(mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>HI State:<b>",as.character(HI_State), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(mus_caught),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "mus_caught",
                   opacity = 0.5) %>%
  addCircleMarkers(data = sub_full_Data_Crypto_caught, 
                   col = ~data_col_Crypto_caught(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>HI State:<b>",as.character(HI_State), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Crypto_mus_caught),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Crypto_caught",
                   opacity = 0.5) %>%
  addLegend("bottomleft", 
            pal = data_col_Transect, 
            title = "Transect",
            values = sub_full_Data_Transect$Transect, 
            group = "Transect") %>%
  addLegend("bottomleft", 
            pal = data_col_Sex, 
            title = "Sex",
            values = sub_full_Data_Sex$Sex, 
            group = "Sex") %>%
  addLegend("bottomleft", 
            pal = data_col_HI_State, 
            title = "HI_State",
            values = sub_full_Data_HI_State$HI_State,  
            group = "HI_State") %>%
  addLegend("bottomleft",
            pal = data_col_Year, 
            title = "Year",
            values = sub_full_Data_Year$Year,
            group = "Year") %>%
  addLegend("bottomleft", 
            pal = data_col_OP_Level, 
            title = "Oocyst Prediction",
            values = sub_full_Data_OP_Level$Oocyst_Level,  
            group = "Oocyst_Level",
            opacity = 1) %>%
  addLegend("bottomleft", 
            pal = data_col_Tested_by, 
            title = "Tested by",
            values = sub_full_Data_Tested_by$Tested_by,  
            group = "Tested_by",
            opacity = 1) %>%
  addLegend("bottomleft", 
            pal = data_col_Machine, 
            title = "Machine",
            values = sub_full_Data_Machine$Machine,  
            group = "Machine",
            opacity = 1) %>%
  addLegend("bottomleft", 
            pal = data_col_Measurements, 
            title = "Measurements",
            values = sub_full_Data_Measurements$Measurements,  
            group = "Measurements",
            opacity = 1) %>%
  addLegend("bottomleft", 
            pal = data_col_Ct_bel, 
            title = "Ct believable",
            values = sub_full_Data_Ct_bel$Ct_believable,  
            group = "Ct_believable",
            opacity = 1) %>%
  addLegend("bottomleft", 
            pal = data_col_Ct_mean, 
            title = "Ct_mean",
            values = sub_full_Data_Ct_mean$Ct_Level,  
            group = "Ct_mean",
            opacity = 1) %>%
  addLegend("bottomleft", 
            pal = data_col_mus_caught, 
            title = "mus caught",
            values = sub_full_Data_mus_caught$mus_caught,  
            group = "mus_caught",
            opacity = 1) %>%
  addLegend("bottomleft", 
            pal = data_col_Crypto_caught, 
            title = "Crypto caught",
            values = sub_full_Data_Crypto_caught$Crypto_mus_caught,  
            group = "Crypto_caught",
            opacity = 1) %>%
  addLayersControl(#baseGroups = c("Transect", "Sex", "HI_State", "Year", "Oocyst_Level", "Tested_by", "Machine"),
    overlayGroups = c("Transect","Sex", "HI_State", "Year", "Oocyst_Level", "Tested_by", "Machine", "Measurements", "Ct_believable", "Ct_mean", "mus_caught", "Crypto_caught"),
    options = layersControlOptions(collapsed = F))


  












