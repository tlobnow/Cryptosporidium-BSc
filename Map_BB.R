library(ggplot2)
library(dplyr)
library(leaflet)
library(sp)
library(RColorBrewer)
display.brewer.all()

#### FULL1 DATASET, contains filtered Samples, (Crypto-positive) ##############
data <- read.csv("nf_full1_Corrected.csv") %>%
  filter(Year == "2019")
glimpse(data)
data <- data[complete.cases(data),]
data$Latitude   <- as.numeric(data$Latitude)
data$Longitude  <- as.numeric(data$Longitude)
data.SP         <- SpatialPointsDataFrame(data[,c(9,10)], data[,-c(9,10)]) 
data$HI_Lvl <- cut(data$HI, c(0, 0.25, 0.5, 0.75, 1), include.lowest = T ,
                   labels = c('0 to 0.25', '0.25 to 0.5', '0.75 to 1', '1'))
data_col = colorFactor(palette = 'RdYlGn', data$HI_Lvl)#puts it into a spatial data frame (negatives of those columns as well)

map <- leaflet() %>%
        addTiles() %>%
        addCircleMarkers(data = data, 
                         lng = ~Longitude, 
                         lat = ~Latitude,
                         color = ~data_col(HI_Lvl),
                         popup= ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                       "<b>HI:<b>",      as.character(HI), "<br>",
                                       "<b>Year:<b>",    as.character(Year),"<br>",
                                       "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                       "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                       "<b>Sex:<b>", Sex, "<br>",
                                       sep=" ")) %>%
        addLegend('bottomleft', pal = data_col, values = data$HI_Lvl,
                  title = 'HI Index<br>of Samples<br>all Years',
                  opacity = 1)
map



#### NON-FILTERED DATASET, all Samples ########################################
data1 <- read.csv("nf_full1_Corrected.csv")
glimpse(data1)
data1 <- data1[complete.cases(data1),]
data1$Latitude    <- as.numeric(data1$Latitude)
data1$Longitude   <- as.numeric(data1$Longitude)
data.SP           <- SpatialPointsDataFrame(data1[,c(9,10)], data1[,-c(9,10)])
data1$HI_Lvl <- cut(data1$HI, c(0, 0.25, 0.5, 0.75, 1), include.lowest = T ,
                    labels = c('0 to 0.25', '0.25 to 0.5', '0.75 to 1', '1'))
data1_col = colorFactor(palette = 'RdYlGn', data$HI_Lvl)

nf_map <- leaflet() %>%
          addTiles() %>%
          addCircleMarkers(data = data1, 
                     lng = ~Longitude, 
                     lat = ~Latitude,
                     color = ~data1_col(HI_Lvl),
                     radius = data$HI*10,
                     popup= ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                   "<b>HI:<b>",      as.character(HI), "<br>",
                                   "<b>Year:<b>",    as.character(Year),"<br>",
                                   "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                   "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                   "<b>Sex:<b>", Sex, "<br>",
                                   sep=" ")) %>%
  addLegend('bottomleft', pal = data_col, values = data$HI_Lvl,
            title = 'HI Index<br>of Samples<br>all Years',
            opacity = 1)
nf_map
             

#### FULL1 DATASET, POSITIVE SAMPLES, BY OOCYST PREDICTION LEVEL ###############

data <- read.csv("full1.csv")
glimpse(data)
data <- data[complete.cases(data),]
data$Latitude   <- as.numeric(data$Latitude)
data$Longitude  <- as.numeric(data$Longitude)
data.SP         <- SpatialPointsDataFrame(data[,c(11,12)], data[,-c(11,12)]) 
data$Oocyst_Predict_Lvl <- cut(data$Oocyst_Predict, c(0, 10^3, 10^4, 10^5, 10^6, 10^7, 10^14), include.lowest = T ,
                   labels = c('0 to 10^3', '10^3 to 10^4', '10^4 to 10^5', '10^5 to 10^6', '10^6 to 10^7', '> 10^7'))
data_col = colorFactor(palette = 'RdYlGn', data$Oocyst_Predict_Lvl)#puts it into a spatial data frame (negatives of those columns as well)

map_Oocyst <- leaflet() %>%
  addTiles() %>%
  addCircleMarkers(data = data, 
                   lng = ~Longitude, 
                   lat = ~Latitude,
                   color = ~data_col(Oocyst_Predict_Lvl),
                   popup= ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                 "<b>HI:<b>",      as.character(HI), "<br>",
                                 "<b>Year:<b>",    as.character(Year),"<br>",
                                 "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                 "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                 "<b>Sex:<b>", Sex, "<br>",
                                 sep=" ")) %>%
  addLegend('bottomleft', pal = data_col, values = data$Oocyst_Predict_Lvl,
            title = 'Oocyst Prediction<br>of Samples<br>',
            opacity = 20)
map_Oocyst


#### FULL1 DATASET, POSITIVE SAMPLES, BY YEAR ##############
data <- read.csv("full1.csv")
data %>%
  filter(Year == 2016)
glimpse(data)
data <- data[complete.cases(data),]
data$Latitude   <- as.numeric(data$Latitude)
data$Longitude  <- as.numeric(data$Longitude)
data.SP         <- SpatialPointsDataFrame(data[,c(11,12)], data[,-c(11,12)]) 

data_col = colorFactor(palette = 'Spectral', data$Year)   #puts it into a spatial data frame (negatives of those columns as well)
# RdYlBu
# Accent
# Greens


map_Year <- leaflet() %>%
  addTiles() %>%
  addCircleMarkers(data = data, 
                   lng = ~Longitude, 
                   lat = ~Latitude,
                   color = ~data_col(Year),
                   popup= ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                 "<b>HI:<b>",      as.character(HI), "<br>",
                                 "<b>Year:<b>",    as.character(Year),"<br>",
                                 "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                 "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                 "<b>Sex:<b>", Sex, "<br>",
                                 sep=" ")) %>%
  addLegend('bottomleft', pal = data_col, values = data$Year,
            title = 'Year',
            opacity = 1)
map_Year

# FILTERED DATA FRAME AS MAP BASE
data <- read.csv("full1.csv")
f_data <- data %>%
  filter(Year == 2016)
glimpse(f_data)
f_data <- f_data[complete.cases(f_data),]
f_data$Latitude   <- as.numeric(f_data$Latitude)
f_data$Longitude  <- as.numeric(f_data$Longitude)
f_data.SP         <- SpatialPointsDataFrame(f_data[,c(11,12)], f_data[,-c(11,12)]) 

f_data_col = colorFactor(palette = 'RdYlGn', f_data$Year)   #puts it into a spatial f_data frame (negatives of those columns as well)

map_Year <- leaflet() %>%
  addTiles() %>%
  addCircleMarkers(data = f_data, 
                   lng = ~Longitude, 
                   lat = ~Latitude,
                   color = ~f_data_col(Year),
                   popup= ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                 "<b>HI:<b>",      as.character(HI), "<br>",
                                 "<b>Year:<b>",    as.character(Year),"<br>",
                                 "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                 "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                 "<b>Sex:<b>", Sex, "<br>",
                                 sep=" ")) %>%
  addLegend('bottomleft', pal = f_data_col, values = f_data$Year,
            title = 'Year',
            opacity = 1)
map_Year



#### ANNIE LEAFLET TIPS #######################################################
#Mittelpunkt_Europa <- data.frame(lat=numeric(), long= numeric())
#Mittelpunkt_Europa[1,] <- c(53.5511,9.9937)
#coordinates(Mittelpunkt_Europa) <- ~long + lat
#Mittelpunkt_Europa


#visualize the buffer on a map
#mapHamburg_Buffer <- leaflet() %>%
 # addTiles() %>%  # Add default OpenStreetMap map tiles
  #addPolygons(data=Hamburg_Buffer)
#mapHamburg_Buffer