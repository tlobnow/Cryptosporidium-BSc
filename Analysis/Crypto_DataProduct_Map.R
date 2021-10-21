library(tidyverse)
library(data.table)
library(leaflet)
library(ggplot2)
library(RColorBrewer) #display.brewer.all()
library(colorRamps)
library(cowplot)
library(uwot)
library(leaflet.extras)
library(sp)
library(htmltools)
library(stringr)
library(patchwork)

# Load Palette ####
r <- c(0, 64, 128, 179, 217, 255)
g <- c(0, 12, 25, 25, 12,  0)
b <- c(255, 249, 243, 191,  95,   0)
  
beach <- function (n, name = c("beach.colors")) 
  {
    beach.colors = rgb(r,g,b,maxColorValue = 255)
    name = match.arg(name)
    orig = eval(parse(text = name))
    rgb = t(col2rgb(orig))
    temp = matrix(NA, ncol = 3, nrow = n)
    x = seq(0, 1, , length(orig))
    xg = seq(0, 1, , n)
    for (k in 1:3) {
      hold = spline(x, rgb[, k], n = n)$y
      hold[hold < 0] = 0
      hold[hold > 255] = 255
      temp[, k] = round(hold)
    }
    palette = rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
    palette
  }


# load Data, add necessary columns ####
Crypto_Detection <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/Crypto_Detection.csv") %>% select(-X)

Crypto_pull <- Crypto_Detection %>% count(Longitude)%>% mutate(mus_caught = n) %>% select(Longitude, mus_caught) %>% arrange(mus_caught)
Crypto_pull_pos <- Crypto_Detection %>% filter(Ct_mean > 0) %>% count(Longitude) %>% mutate(Crypto_mus_caught = n) %>% select(Longitude, Crypto_mus_caught) %>% arrange(Crypto_mus_caught)
Crypto_Detection_1 <- left_join(Crypto_Detection, Crypto_pull)
Crypto_Detection_2 <- left_join(Crypto_Detection_1, Crypto_pull_pos)
Crypto_Detection <- Crypto_Detection_2
rm(Crypto_Detection_1)
rm(Crypto_Detection_2)
rm(Crypto_pull)
rm(Crypto_pull_pos)

# Map Data ####
map <- Crypto_Detection %>%
  leaflet() %>%
  addProviderTiles("CartoDB") %>%
  setView(lat = 52.520007, lng =13.404954, zoom = 7)
map

# for now this excludes HZ21 due to missing genotyping info
Crypto_Detection <- Crypto_Detection %>% filter(HI  != "NA", Longitude != "NA", Latitude  != "NA")

Crypto_Detection <- Crypto_Detection %>% replace_na(list(Crypto_mus_caught = 0))
Crypto_Detection[,'HI']=format(round(Crypto_Detection[,'HI'],2), nsmall=2)

Crypto_Detection$mus_Level <-  cut(Crypto_Detection$mus_caught, c(0, 1, 5, 10, 15, 20, 21), include.lowest = T ,
                                       labels = c('1 mouse', 'up to 5 mice', 'up to 10 mice', 'up to 15 mice', 'up to 20 mice', '21 mice'))

data_col_mus        = colorFactor(matlab.like(6), Crypto_Detection$mus_caught)
data_col_mus_Level  = colorFactor(matlab.like(6), Crypto_Detection$mus_Level)

mus_1 <- Crypto_Detection %>% filter(mus_caught == 1)
mus_5 <- Crypto_Detection %>% filter(mus_caught > 1 & mus_caught <= 5)
mus_10 <- Crypto_Detection %>% filter(mus_caught > 5 & mus_caught <= 10)
mus_15 <- Crypto_Detection %>% filter(mus_caught > 10 & mus_caught <= 15)
mus_20 <- Crypto_Detection %>% filter(mus_caught > 15 & mus_caught <= 20)
mus_21 <- Crypto_Detection %>% filter(mus_caught == 21)

# Mus Map ####
map %>%
  addCircleMarkers(data = mus_1, 
                   col = ~data_col_mus(mus_caught),
                   label = ~htmlEscape(mus_caught),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3, 
                   opacity = 0.5,
                   group = "mus_1") %>%
  addCircleMarkers(data = mus_5, 
                   col = ~data_col_mus(mus_caught),
                   label = ~htmlEscape(mus_caught),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 0.5,
                   group = "mus_5") %>%
  addCircleMarkers(data = mus_10, 
                   col = ~data_col_mus(mus_caught),
                   label = ~htmlEscape(mus_caught),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "mus_10") %>%
  addCircleMarkers(data = mus_15, 
                   col = ~data_col_mus(mus_caught),
                   label = ~htmlEscape(mus_caught),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "mus_15") %>%
  addCircleMarkers(data = mus_20, 
                   col = ~data_col_mus(mus_caught),
                   label = ~htmlEscape(mus_caught),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "mus_20") %>%
  addCircleMarkers(data = mus_21, 
                   col = ~data_col_mus(mus_caught),
                   label = ~htmlEscape(mus_caught),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "mus_21") %>%
  addLegend("bottomleft", 
            pal = data_col_mus_Level, 
            title = "Number of mice caught",
            values = Crypto_Detection$mus_Level, 
            group = c('1 mouse', 'up to 5 mice', 'up to 10 mice', 'up to 15 mice', 'up to 20 mice', '21 mice'),
            opacity = 1) %>%
  addLayersControl(overlayGroups = c('mus_1', 'mus_5', 'mus_10', 'mus_15', 'mus_20', 'mus_21'),
                   options = layersControlOptions(collapsed = F))


# Crypto_Map ####
data_col_mus_crypto = colorFactor(matlab.like(8), Crypto_Detection$Crypto_mus_caught)

Crypto_Infections_1 <- Crypto_Detection %>% filter(Crypto_mus_caught == 1)
Crypto_Infections_2 <- Crypto_Detection %>% filter(Crypto_mus_caught == 2)
Crypto_Infections_3 <- Crypto_Detection %>% filter(Crypto_mus_caught == 3)
Crypto_Infections_4 <- Crypto_Detection %>% filter(Crypto_mus_caught == 4)
Crypto_Infections_5 <- Crypto_Detection %>% filter(Crypto_mus_caught == 5)
Crypto_Infections_6 <- Crypto_Detection %>% filter(Crypto_mus_caught == 6)
Crypto_Infections_7 <- Crypto_Detection %>% filter(Crypto_mus_caught == 7)
Crypto_Infections_8 <- Crypto_Detection %>% filter(Crypto_mus_caught == 8)

map %>%
  addCircleMarkers(data = Crypto_Infections_1, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Crypto_Infections_1") %>%
  addCircleMarkers(data = Crypto_Infections_2, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Crypto_Infections_2") %>%
  addCircleMarkers(data = Crypto_Infections_3, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Crypto_Infections_3") %>%
  addCircleMarkers(data = Crypto_Infections_4, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Crypto_Infections_4") %>%
  addCircleMarkers(data = Crypto_Infections_5, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Crypto_Infections_5") %>%
  addCircleMarkers(data = Crypto_Infections_6, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Crypto_Infections_6") %>%
  addCircleMarkers(data = Crypto_Infections_7, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Crypto_Infections_7") %>%
  addCircleMarkers(data = Crypto_Infections_8, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Crypto_Infections_8") %>%
  addLegend("bottomright", 
            pal = data_col_mus_crypto, 
            title = "Crypto caught",
            values = Crypto_Detection$Crypto_mus_caught, 
            group = c('Crypto_Infections_1', 
                      'Crypto_Infections_2', 
                      'Crypto_Infections_3', 
                      'Crypto_Infections_4', 
                      'Crypto_Infections_5',
                      'Crypto_Infections_6',
                      'Crypto_Infections_7',
                      'Crypto_Infections_8'),
            opacity = 1) %>%
  addLayersControl(overlayGroups = c('Crypto_Infections_1', 
                                     'Crypto_Infections_2', 
                                     'Crypto_Infections_3', 
                                     'Crypto_Infections_4', 
                                     'Crypto_Infections_5',
                                     'Crypto_Infections_6',
                                     'Crypto_Infections_7',
                                     'Crypto_Infections_8'),
                   options = layersControlOptions(collapsed = T))

# HI Map ####

Crypto_Detection$HI <- as.numeric(Crypto_Detection$HI)
data_col_HI       = colorFactor(beach(6), Crypto_Detection$HI)
data_col_HI_Level = colorFactor(beach(6), Crypto_Detection$HI_Level)
data_col_Sequenced = colorFactor(matlab.like(2), Crypto_Detection$Sequenced)


Crypto_Detection$HI_Level <-  cut(Crypto_Detection$HI, c(0, 0.001, 0.250, 0.500, 0.750, 0.999, 1), include.lowest = T ,
                                  labels = c('HI = 0', 'HI < 0.25', 'HI < 0.5', 'HI < 0.75', 'HI < 1', 'HI = 1'))


HI_0 <- Crypto_Detection %>% filter(HI < 0.01)
HI_below_0.25 <- Crypto_Detection %>% filter(HI > 0.01, HI <= 0.25)
HI_below_0.5 <- Crypto_Detection %>% filter(HI > 0.25, HI <= 0.5)
HI_below_0.75 <- Crypto_Detection %>% filter(HI > 0.5, HI <= 0.75)
HI_below_1 <- Crypto_Detection %>% filter(HI > 0.75, HI < 1)
HI_equal_1 <- Crypto_Detection %>% filter(HI > 0.99)

map %>%
  addCircleMarkers(data = HI_0, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "HI_0") %>%
  addCircleMarkers(data = HI_below_0.25, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "HI_below_0.25") %>%
  addCircleMarkers(data = HI_below_0.5, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "HI_below_0.5") %>%
  addCircleMarkers(data = HI_below_0.75, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "HI_below_0.75") %>%
  addCircleMarkers(data = HI_below_1, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "HI_below_1") %>%
  addCircleMarkers(data = HI_equal_1, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "HI_equal_1") %>%
  addLegend("bottomleft", 
            pal = data_col_HI_Level, 
            title = "HI",
            values = Crypto_Detection$HI_Level, 
            group = c('HI = 0', 'HI < 0.25', 'HI < 0.5', 'HI < 0.75', 'HI < 1', 'HI = 1'),
            opacity = 1) #%>%
  #addLayersControl(overlayGroups = c("HI_0", 
 # "HI_below_0.25", 
#  "HI_below_0.5", 
#  "HI_below_0.75", 
#  "HI_below_1", 
#  "HI_equal_1"),
  #                 options = layersControlOptions(collapsed = F))


# Map Crypto Infections per mice caught ####

Crypto_Detection <- Crypto_Detection %>%
  mutate(Infection_Rate = Crypto_mus_caught / mus_caught)
Crypto_Detection[,'Infection_Rate']=format(round(Crypto_Detection[,'Infection_Rate'],2),nsmall=2)
Crypto_Detection$Infection_Rate <- as.numeric(Crypto_Detection$Infection_Rate)

Crypto_Detection$Infection_Level <-  cut(Crypto_Detection$Infection_Rate, c(0, 0.01, 0.25, 0.5, 0.75, 1), include.lowest = T ,
                                      labels = c('0 %', '< 25 %', '< 50 %', '< 75 %', '100 %'))

data_col_mus_infection_rate = colorFactor(matlab.like(6), Crypto_Detection$Infection_Level)

Infection_Rate_0    <- Crypto_Detection %>% filter(Infection_Level == '0 %')
Infection_Rate_25   <- Crypto_Detection %>% filter(Infection_Level == '< 25 %')
Infection_Rate_50   <- Crypto_Detection %>% filter(Infection_Level == '< 50 %')
Infection_Rate_75   <- Crypto_Detection %>% filter(Infection_Level == '< 75 %')
Infection_Rate_equal_100   <- Crypto_Detection %>% filter(Infection_Level == '100 %')

map <- Crypto_Detection %>%
  leaflet() %>%
  addProviderTiles("CartoDB") %>%
  setView(lat = 52.520007, lng =13.404954, zoom = 7) 



map %>%
  addCircleMarkers(data = Crypto_Detection, 
                   col = ~data_col_mus_infection_rate(Infection_Level),
                   label = ~htmlEscape(mus_caught),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   opacity = 0.3,
                   radius = 3,
                   group = "Crypto_Detection",
                   clusterOptions = markerClusterOptions()) %>%
  addCircleMarkers(data = Infection_Rate_0, 
                   col = ~data_col_mus_infection_rate(Infection_Level),
                   label = ~htmlEscape(mus_caught),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   opacity = 0.3,
                   radius = 3,
                   group = "Infection_Rate_0") %>%
  addCircleMarkers(data = Infection_Rate_25, 
                   col = ~data_col_mus_infection_rate(Infection_Level),
                   label = ~htmlEscape(mus_caught),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Infection_Rate_25") %>%
  addCircleMarkers(data = Infection_Rate_50, 
                   col = ~data_col_mus_infection_rate(Infection_Level),
                   label = ~htmlEscape(mus_caught),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Infection_Rate_50") %>%
  addCircleMarkers(data = Infection_Rate_75, 
                   col = ~data_col_mus_infection_rate(Infection_Level),
                   label = ~htmlEscape(mus_caught),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Infection_Rate_75") %>%
  addCircleMarkers(data = Infection_Rate_equal_100, 
                   col = ~data_col_mus_infection_rate(Infection_Level),
                   label = ~htmlEscape(mus_caught),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "Infection_Rate_equal_100") %>%
  addLegend("bottomleft", 
            pal = data_col_mus_infection_rate, 
            title = "Infection Rate",
            values = Crypto_Detection$Infection_Level, 
            group = c('Infection_Rate_0',
                      'Infection_Rate_25',
                      'Infection_Rate_50',
                      'Infection_Rate_75',
                      'Infection_Rate_equal_100'),
            opacity = 1) %>%
  addLayersControl(overlayGroups = c('Crypto_Detection',
                                     'Infection_Rate_0',
                                     'Infection_Rate_25',
                                     'Infection_Rate_50',
                                     'Infection_Rate_75',
                                     'Infection_Rate_equal_100'),
                   options = layersControlOptions(collapsed = F))



################################################################################
## Visualizing sequenced Samples

Crypto_Detection <- Crypto_Detection %>%
  mutate(Sequenced = Mouse_ID %in% c("AA_0144", "AA_0325", "AA_0689", "AA_0209",
                                     "AA_0282", "AA_0793", "AA_0667", "AA_0805"))

data_col_Sequenced = colorFactor(brewer.pal(3, "Spectral"), Crypto_Detection$Sequenced)

#display.brewer.all()


HI_0 <- Crypto_Detection %>% filter(HI < 0.01)
HI_below_0.25 <- Crypto_Detection %>% filter(HI > 0.01, HI <= 0.25)
HI_below_0.5 <- Crypto_Detection %>% filter(HI > 0.25, HI <= 0.5)
HI_below_0.75 <- Crypto_Detection %>% filter(HI > 0.5, HI <= 0.75)
HI_below_1 <- Crypto_Detection %>% filter(HI > 0.75, HI < 1)
HI_equal_1 <- Crypto_Detection %>% filter(HI > 0.99)
Sequenced <- Crypto_Detection %>% filter(Sequenced == T)
Not_Sequenced <- Crypto_Detection %>% filter(Sequenced == F)


map %>%
  addCircleMarkers(data = HI_0, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   opacity = 0.05,
                   radius = 3,
                   group = "HI_0") %>%
  addCircleMarkers(data = HI_below_0.25, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 0.05,
                   group = "HI_below_0.25") %>%
  addCircleMarkers(data = HI_below_0.5, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 0.05,
                   group = "HI_below_0.5") %>%
  addCircleMarkers(data = HI_below_0.75, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 0.05,
                   group = "HI_below_0.75") %>%
  addCircleMarkers(data = HI_below_1, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 0.05,
                   group = "HI_below_1") %>%
  addCircleMarkers(data = HI_equal_1, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 0.05,
                   group = "HI_equal_1") %>%
  addCircleMarkers(data = Sequenced, 
                   col = ~data_col_Sequenced(Sequenced),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   opacity = 3,
                   radius = 3,
                   group = "Sequenced") %>%
  addCircleMarkers(data = Not_Sequenced, 
                   col = ~data_col_Sequenced(Sequenced),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 3,
                   group = "Not_Sequenced") %>%
  addLegend("bottomleft", 
            pal = data_col_HI_Level, 
            title = "HI",
            values = Crypto_Detection$HI_Level, 
            group = c('HI = 0', 'HI < 0.25', 'HI < 0.5', 'HI < 0.75', 'HI < 1', 'HI = 1'),
            opacity = 1) %>%
addLayersControl("bottomright",
                 baseGroups = c("Sequenced",
                                "Not_Sequenced"),
                 overlayGroups = c("HI_0", 
                                   "HI_below_0.25", 
                                   "HI_below_0.5", 
                                   "HI_below_0.75", 
                                   "HI_below_1", 
                                   "HI_equal_1"),
                 options = layersControlOptions(collapsed = T))




###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Ct Map
Crypto_Detection <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/Crypto_Detection.csv") %>% select(-X)

Crypto_pull <- Crypto_Detection %>% count(Longitude)%>% mutate(mus_caught = n) %>% select(Longitude, mus_caught) %>% arrange(mus_caught)
Crypto_pull_pos <- Crypto_Detection %>% filter(Ct_mean > 0) %>% count(Longitude) %>% mutate(Crypto_mus_caught = n) %>% select(Longitude, Crypto_mus_caught) %>% arrange(Crypto_mus_caught)
Crypto_Detection_1 <- left_join(Crypto_Detection, Crypto_pull)
Crypto_Detection_2 <- left_join(Crypto_Detection_1, Crypto_pull_pos)
Crypto_Detection <- Crypto_Detection_2
rm(Crypto_Detection_1)
rm(Crypto_Detection_2)
rm(Crypto_pull)
rm(Crypto_pull_pos)

Crypto_Detection <- Crypto_Detection %>% filter(Year == 2021,
                                                Ct_mean > 0)

map <- Crypto_Detection %>%
  leaflet() %>%
  addProviderTiles("CartoDB") %>%
  setView(lat = 52.520007, lng =13.404954, zoom = 8) 
Crypto_Detection$Ct_mean <- round(Crypto_Detection$Ct_mean, digits = 2)

Ct_bl_30 <- Crypto_Detection %>% filter(Ct_mean <= 30 & Ct_mean > 0)
Ct_bl_34 <- Crypto_Detection %>% filter(Ct_mean <= 34 & Ct_mean > 30)
Ct_bl_36 <- Crypto_Detection %>% filter(Ct_mean <= 36 & Ct_mean > 34)
Ct_over_36 <- Crypto_Detection %>% filter(Ct_mean > 36)

#Crypto_Detection$Ct_mean <- round(Crypto_Detection$Ct_mean, digits = 2)
Crypto_Detection$Ct_Level  <-  cut(Crypto_Detection$Ct_mean,   c(0, 1, 28, 30, 32, 34, 36, 40), include.lowest = T , labels = c('0',' < 28  ', '28 - 30', '30 - 32', '32 - 34', '34 - 36', ' > 36'))
Crypto_Detection <- Crypto_Detection %>% mutate(Ct_Level = case_when((Ct_mean <= 30 & Ct_mean > 0) ~ "<= 30",
                                                                     (Ct_mean <= 34 & Ct_mean > 30) ~ "30 - 34",
                                                                     (Ct_mean <= 36 & Ct_mean > 34) ~ "34 - 36",
                                                                     (Ct_mean > 36) ~ "> 36"))

data_col_Ct_Level = colorFactor(brewer.pal(4, "YlOrRd" ), Crypto_Detection$Ct_mean, reverse = T)


map %>%
  addCircleMarkers(data = HI_0, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   opacity = 0.05,
                   radius = 3,
                   group = "HI_0") %>%
  addCircleMarkers(data = HI_below_0.25, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 0.05,
                   group = "HI_below_0.25") %>%
  addCircleMarkers(data = HI_below_0.5, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 0.05,
                   group = "HI_below_0.5") %>%
  addCircleMarkers(data = HI_below_0.75, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 0.05,
                   group = "HI_below_0.75") %>%
  addCircleMarkers(data = HI_below_1, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 0.05,
                   group = "HI_below_1") %>%
  addCircleMarkers(data = HI_equal_1, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 0.05,
                   group = "HI_equal_1") %>%
  addCircleMarkers(data = Ct_over_36, 
                   color = ~data_col_Ct_Level(Ct_mean),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>", as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  "<b>Oocyst Prediction:<b>",  as.character(Oocyst_Predict), "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 5,
                   group = "Ct_over_36") %>%
  addCircleMarkers(data = Ct_over_36, 
                   color = ~data_col_Ct_Level(Ct_mean),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>", as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  "<b>Oocyst Prediction:<b>",  as.character(Oocyst_Predict), "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 5,
                   group = "Ct_over_36") %>%
  addCircleMarkers(data = Ct_bl_36, 
                   color = ~data_col_Ct_Level(Ct_mean),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>", as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  "<b>Oocyst Prediction:<b>",  as.character(Oocyst_Predict), "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 5,
                   group = "Ct_bl_36") %>%
  addCircleMarkers(data = Ct_bl_34, 
                   label = ~htmlEscape(Mouse_ID),
                   color = ~data_col_Ct_Level(Ct_mean),
                   popup = ~paste("<b>Mouse_ID:<b>", as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  "<b>Oocyst Prediction:<b>",  as.character(Oocyst_Predict), "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 5,
                   group = "Ct_bl_34") %>%
  addCircleMarkers(data = Ct_bl_30, 
                   label = ~htmlEscape(Mouse_ID),
                   color = ~data_col_Ct_Level(Ct_mean),
                   popup = ~paste("<b>Mouse_ID:<b>", as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  "<b>Oocyst Prediction:<b>",  as.character(Oocyst_Predict), "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 5,
                   group = "Ct_bl_30") %>%
  #addLegend("bottomleft", 
  #          pal = data_col_Ct_Level, 
  #          title = "Measured Ct",
  #          values = Crypto_Detection$Ct_mean, 
  #          group = c("Ct_bl_30", "Ct_bl_34", "Ct_bl_36", "Ct_over_36"),
  #          opacity = 1) %>%
  addLayersControl(overlayGroups = c("Ct_bl_30", "Ct_bl_34", "Ct_bl_36", "Ct_over_36"),
                   options = layersControlOptions(collapsed = T))





################################################################################
################################################################################
################################################################################
################################################################################
## Oocyst Predictions




Crypto_Detection <- Crypto_Detection %>% filter(Year == 2021)
map <- Crypto_Detection %>%
  leaflet() %>%
  addProviderTiles("CartoDB") %>%
  setView(lat = 52.520007, lng =13.404954, zoom = 8) 

Crypto_Detection$Ct_Level <- cut(Crypto_Detection$Ct_mean, c(0, 1, 26, 28, 30, 32, 34, 36, 38, 40), include.lowest = F , labels = c('0', '< 26', ' 26 - 28', '28 - 30', '30 - 32', ' 32 - 34', '34 - 36', '36 - 38', ' 38 - 40'))
data_col_Oocysts = colorFactor(matlab.like(5), Crypto_Detection$Oocyst_Predict)

Opred_100 <- Crypto_Detection %>% filter(Oocyst_Predict <= 100 & Oocyst_Predict > 0)
Opred_1000 <- Crypto_Detection %>% filter(Oocyst_Predict <= 1000 & Oocyst_Predict > 100)
Opred_10.000 <- Crypto_Detection %>% filter(Oocyst_Predict <= 10000 & Oocyst_Predict > 1000)
Opred_100.000 <- Crypto_Detection %>% filter(Oocyst_Predict <= 100000 & Oocyst_Predict > 10000)


SOTA <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv") %>% select(-X)


  
  