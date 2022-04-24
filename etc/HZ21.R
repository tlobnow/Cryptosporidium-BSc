library(tidyverse)
library(leaflet)
library(colorRamps)
library(htmltools)
library(data.table)




Trap21 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ21_Trap.csv")

SOTA <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv")

#Crypto_Detection <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/Crypto_Detection.csv")
#Crypto21 <- Crypto_Detection %>% filter(Year == 2021)
#Crypto21 <- Crypto21[Crypto21$Mouse_ID %like% "AA_", ]
#Crypto21 <- Crypto21 %>% filter(ILWE_DNA_Content_ng.microliter >= 0)

#write.csv(Crypto21, "HZ21_Crypto_qPCR_List.csv")


Trap21 %>% count(Mus_caught) %>% filter(!is.na(Mus_caught)) %>% sum()



map <- Trap21 %>%
  leaflet() %>%
  addProviderTiles("CartoDB") %>%
  setView(lat = 52.520007, lng =13.404954, zoom = 6)

Trap21$mus_Level <-  cut(Trap21$Mus_caught, c(0, 1, 5, 10, 15, 20, 30), include.lowest = T , labels = c('1 mouse', 'up to 5 mice', 'up to 10 mice', 'up to 15 mice', 'up to 20 mice', ' up to 30 mice'))

data_col_mus        = colorFactor(matlab.like(6), Trap21$Mus_caught)
data_col_mus_Level  = colorFactor(matlab.like(6), Trap21$mus_Level)

mus_1 <- Trap21 %>% filter(Mus_caught == 1)
mus_5 <- Trap21 %>% filter(Mus_caught > 1 & Mus_caught <= 5)
mus_10 <- Trap21 %>% filter(Mus_caught > 5 & Mus_caught <= 10)
mus_15 <- Trap21 %>% filter(Mus_caught > 10 & Mus_caught <= 15)
mus_20 <- Trap21 %>% filter(Mus_caught > 15 & Mus_caught <= 20)
mus_30 <- Trap21 %>% filter(Mus_caught > 20 & Mus_caught <= 30)

map %>%
  addCircleMarkers(data = mus_1, 
                   col = ~data_col_mus(Mus_caught),
                   label = ~htmlEscape(Address),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Address:<b>", as.character(Address), "<br>",
                                  "<b>Mus caught:<b>",  as.character(Mus_caught), "<br>",
                                  "<b>Rodents caught:<b>",  as.character(Rodents_caught), "<br>",
                                  sep=" "),
                   radius = 3, 
                   opacity = 1,
                   group = "mus_1") %>%
  addCircleMarkers(data = mus_5, 
                   col = ~data_col_mus(Mus_caught),
                   label = ~htmlEscape(Address),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Address:<b>", as.character(Address), "<br>",
                                  "<b>Mus caught:<b>",  as.character(Mus_caught), "<br>",
                                  "<b>Rodents caught:<b>",  as.character(Rodents_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   opacity = 1,
                   group = "mus_5") %>%
  addCircleMarkers(data = mus_10, 
                   col = ~data_col_mus(Mus_caught),
                   label = ~htmlEscape(Address),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Address:<b>", as.character(Address), "<br>",
                                  "<b>Mus caught:<b>",  as.character(Mus_caught), "<br>",
                                  "<b>Rodents caught:<b>",  as.character(Rodents_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "mus_10") %>%
  addCircleMarkers(data = mus_15, 
                   col = ~data_col_mus(Mus_caught),
                   label = ~htmlEscape(Address),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Address:<b>", as.character(Address), "<br>",
                                  "<b>Mus caught:<b>",  as.character(Mus_caught), "<br>",
                                  "<b>Rodents caught:<b>",  as.character(Rodents_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "mus_15") %>%
  addCircleMarkers(data = mus_20, 
                   col = ~data_col_mus(Mus_caught),
                   label = ~htmlEscape(Address),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Address:<b>", as.character(Address), "<br>",
                                  "<b>Mus caught:<b>",  as.character(Mus_caught), "<br>",
                                  "<b>Rodents caught:<b>",  as.character(Rodents_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "mus_20") %>%
  addCircleMarkers(data = mus_30, 
                   col = ~data_col_mus(Mus_caught),
                   label = ~htmlEscape(Address),
                   popup = ~paste("<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>Address:<b>", as.character(Address), "<br>",
                                  "<b>Mus caught:<b>",  as.character(Mus_caught), "<br>",
                                  "<b>Rodents caught:<b>",  as.character(Rodents_caught), "<br>",
                                  sep=" "),
                   radius = 3,
                   group = "mus_30") %>%
    addLegend("bottomleft", 
              pal = data_col_mus_Level, 
              title = "Number of mice caught",
              values = Trap21$mus_Level, 
              group = c('1 mouse', 'up to 5 mice', 'up to 10 mice', 'up to 15 mice', 'up to 20 mice', 'up to 30 mice'),
             opacity = 1) %>%
  addLayersControl(overlayGroups = c('mus_1', 
                                     'mus_5', 
                                     'mus_10', 
                                     'mus_15', 
                                     'mus_20', 
                                     'mus_30'),
                   options = layersControlOptions(collapsed = T))





