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

Crypto_Detection <- read.csv("Crypto_Detection.csv")

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

Crypto_Detection_mus <- Crypto_Detection %>% filter(HI  != "NA", Longitude != "NA", Latitude  != "NA")

Crypto_Detection <- Crypto_Detection %>% replace_na(list(Crypto_mus_caught = 0))
Crypto_Detection[,'HI']=format(round(Crypto_Detection[,'HI'],2),nsmall=2)

Crypto_Detection_mus$mus_Level <-  cut(Crypto_Detection_mus$mus_caught, c(0, 1, 5, 10, 15, 20, 21), include.lowest = T ,
                                       labels = c('1 mouse', '< 5 mice', '< 10 mice', '< 15 mice', '< 20 mice', '21 mice'))


data_col_mus        = colorFactor(matlab.like(6), Crypto_Detection_mus$mus_caught)
data_col_mus_Level  = colorFactor(matlab.like(6), Crypto_Detection_mus$mus_Level)
data_col_mus_crypto = colorFactor(matlab.like(8), Crypto_Detection_mus$Crypto_mus_caught)


mus_1 <- Crypto_Detection_mus %>% filter(mus_caught == 1)
mus_5 <- Crypto_Detection_mus %>% filter(mus_caught > 1 & mus_caught <= 5)
mus_10 <- Crypto_Detection_mus %>% filter(mus_caught > 5 & mus_caught <= 10)
mus_15 <- Crypto_Detection_mus %>% filter(mus_caught > 10 & mus_caught <= 15)
mus_20 <- Crypto_Detection_mus %>% filter(mus_caught > 15 & mus_caught <= 20)
mus_21 <- Crypto_Detection_mus %>% filter(mus_caught == 21)

# Mus Map ####
map %>%
  addCircleMarkers(data = mus_1, 
                   col = ~data_col_mus(mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3, 
                   opacity = 0.5,
                   group = "mus_1") %>%
  addCircleMarkers(data = mus_5, 
                   col = ~data_col_mus(mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 0.5,
                   group = "mus_5") %>%
  addCircleMarkers(data = mus_10, 
                   col = ~data_col_mus(mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "mus_10") %>%
  addCircleMarkers(data = mus_15, 
                   col = ~data_col_mus(mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "mus_15") %>%
  addCircleMarkers(data = mus_20, 
                   col = ~data_col_mus(mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "mus_20") %>%
  addCircleMarkers(data = mus_21, 
                   col = ~data_col_mus(mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                  "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "mus_21") %>%
  addLegend("bottomleft", 
            pal = data_col_mus_Level, 
            title = "Number of mice caught",
            values = Crypto_Detection_mus$mus_Level, 
            group = c('1 mouse', '< 5 mice', '< 10 mice', '< 15 mice', '< 20 mice', '21 mice'),
            opacity = 1) %>%
  addLayersControl(overlayGroups = c('mus_1', 'mus_5', 'mus_10', 'mus_15', 'mus_20', 'mus_21'),
                   options = layersControlOptions(collapsed = F))


# Crypto_Map ####

cryp_1 <- Crypto_Detection_mus %>% filter(Crypto_mus_caught == 1)
cryp_2 <- Crypto_Detection_mus %>% filter(Crypto_mus_caught == 2)
cryp_3 <- Crypto_Detection_mus %>% filter(Crypto_mus_caught == 3)
cryp_4 <- Crypto_Detection_mus %>% filter(Crypto_mus_caught == 4)
cryp_5 <- Crypto_Detection_mus %>% filter(Crypto_mus_caught == 5)
cryp_6 <- Crypto_Detection_mus %>% filter(Crypto_mus_caught == 6)
cryp_7 <- Crypto_Detection_mus %>% filter(Crypto_mus_caught == 7)
cryp_8 <- Crypto_Detection_mus %>% filter(Crypto_mus_caught == 8)

map %>%
  addCircleMarkers(data = cryp_1, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Longitude:<b>", as.character(Longitude), "<br>",
                                  "<b>Latitude:<b>",  as.character(Latitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "cryp_1") %>%
  addCircleMarkers(data = cryp_2, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Longitude:<b>", as.character(Longitude), "<br>",
                                  "<b>Latitude:<b>",  as.character(Latitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "cryp_2") %>%
  addCircleMarkers(data = cryp_3, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Longitude:<b>", as.character(Longitude), "<br>",
                                  "<b>Latitude:<b>",  as.character(Latitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "cryp_3") %>%
  addCircleMarkers(data = cryp_4, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Longitude:<b>", as.character(Longitude), "<br>",
                                  "<b>Latitude:<b>",  as.character(Latitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "cryp_4") %>%
  addCircleMarkers(data = cryp_5, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Longitude:<b>", as.character(Longitude), "<br>",
                                  "<b>Latitude:<b>",  as.character(Latitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "cryp_5") %>%
  addCircleMarkers(data = cryp_6, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Longitude:<b>", as.character(Longitude), "<br>",
                                  "<b>Latitude:<b>",  as.character(Latitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "cryp_6") %>%
  addCircleMarkers(data = cryp_7, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Longitude:<b>", as.character(Longitude), "<br>",
                                  "<b>Latitude:<b>",  as.character(Latitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "cryp_7") %>%
  addCircleMarkers(data = cryp_8, 
                   col = ~data_col_mus_crypto(Crypto_mus_caught),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Longitude:<b>", as.character(Longitude), "<br>",
                                  "<b>Latitude:<b>",  as.character(Latitude), "<br>",
                                  "<b>Crypto Infections:<b>",  as.character(Crypto_mus_caught), "<b>/<b>", as.character(mus_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "cryp_8") %>%
  addLegend("bottomleft", 
            pal = data_col_mus_crypto, 
            title = "Crypto caught",
            values = Crypto_Detection_mus$Crypto_mus_caught, 
            group = c('cryp_1', 
                      'cryp_2', 
                      'cryp_3', 
                      'cryp_4', 
                      'cryp_5',
                      'cryp_6',
                      'cryp_7',
                      'cryp_8'),
            opacity = 1) %>%
  addLayersControl(overlayGroups = c('cryp_1', 
                                     'cryp_2', 
                                     'cryp_3', 
                                     'cryp_4', 
                                     'cryp_5',
                                     'cryp_6',
                                     'cryp_7',
                                     'cryp_8'),
                   options = layersControlOptions(collapsed = F))

# HI Map ####

Crypto_Detection_mus$HI <- as.numeric(Crypto_Detection_mus$HI)
data_col_HI       = colorFactor(beach(6), Crypto_Detection_mus$HI)
data_col_HI_Level = colorFactor(beach(6), Crypto_Detection_mus$HI_Level)


Crypto_Detection_mus$HI_Level <-  cut(Crypto_Detection_mus$HI, c(0, 0.001, 0.250, 0.500, 0.750, 0.999, 1), include.lowest = T ,
                                  labels = c('HI = 0', 'HI < 0.25', 'HI < 0.5', 'HI < 0.75', 'HI < 1', 'HI = 1'))



HI_0 <- Crypto_Detection_mus %>% filter(HI < 0.01)
HI_below_0.25 <- Crypto_Detection_mus %>% filter(HI > 0.01, HI <= 0.25)
HI_below_0.5 <- Crypto_Detection_mus %>% filter(HI > 0.25, HI <= 0.5)
HI_below_0.75 <- Crypto_Detection_mus %>% filter(HI > 0.5, HI <= 0.75)
HI_below_1 <- Crypto_Detection_mus %>% filter(HI > 0.75, HI < 1)
HI_equal_1 <- Crypto_Detection_mus %>% filter(HI > 0.99)

map %>%
  addCircleMarkers(data = HI_0, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>HI State:<b>",as.character(HI), "<br>",
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
                                  "<b>HI State:<b>",as.character(HI), "<br>",
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
                                  "<b>HI State:<b>",as.character(HI), "<br>",
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
                                  "<b>HI State:<b>",as.character(HI), "<br>",
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
                                  "<b>HI State:<b>",as.character(HI), "<br>",
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
                                  "<b>HI State:<b>",as.character(HI), "<br>",
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
            values = Crypto_Detection_mus$HI_Level, 
            group = c('HI = 0', 'HI < 0.25', 'HI < 0.5', 'HI < 0.75', 'HI < 1', 'HI = 1'),
            opacity = 1) #%>%
  #addLayersControl(overlayGroups = c("HI_0", 
  #                                   "HI_below_0.25", 
  #                                   "HI_below_0.5", 
  #                                   "HI_below_0.75", 
  #                                   "HI_below_1", 
  #                                   "HI_equal_1"),
  #                 options = layersControlOptions(collapsed = F))



