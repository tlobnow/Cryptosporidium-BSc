library(ggplot2)
library(dplyr)
library(ggeffects)
library(RColorBrewer) #display.brewer.all()
library(viridis)
library(cowplot)
library(uwot)
library(leaflet)
library(leaflet.extras)
library(sp)
library(htmltools)
library(stringr)


ABI_Best_thSC <- read.csv("ABI_Best_SC.csv")
full1         <- read.csv("full1.csv")
nf_full1      <- read.csv("nf_full1_Corrected.csv")
#all_Samples_adj   <- read.csv("all_Samples_adj.csv")
#all_Samples_adj_trim   <- read.csv("all_Samples_adj_trim.csv")
fd_ABI_log2       <- read.csv("fd_ABI_log2.csv")
Plate6_Candidates <- read.csv("Plate6_Candidates.csv")
new_Plate6        <- read.csv("new_Plate6.csv")
f_ABI_Best_thSC   <-  filter(ABI_Best_thSC, Ct_mean > 0)
linear_model0     <- lm(log2(Amount_Oocysts) ~ Ct_mean, data = f_ABI_Best_thSC)
ggpredict(linear_model0)
summary(linear_model0)



number_ticks <- function(n) {function(limits) pretty(limits, n)}

f_ABI_Best_thSC %>%
  ggplot(aes(log2(Amount_Oocysts), Ct_mean, label = Amount_Oocysts)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_label(nudge_y = 1) +
  theme() +
  scale_x_continuous(breaks=number_ticks(5)) +
  ggtitle("Standard Curve  1:8 Dilution")

#### HI_LEVEL BAR CHARTS DIFF. SAMPLE SIZES ####
full_Gen_HI <- read.csv("alles.csv")
nf_full1 <- read.csv("nf_full1_Corrected.csv")
positive_Samples <- read.csv("positive_Samples.csv")

full_Gen_HI$HI_Level = factor(full_Gen_HI$HI_Level, 
                              levels=c("HI == 0",
                                       "HI < 0.25",
                                       "HI < 0.5",
                                       "HI < 0.75",
                                       "HI < 1",
                                       "HI == 1"))
full_Gen_HI %>%
  ggplot(aes(HI_Level, fill = HI_Level)) +
  geom_histogram(stat = "count")

nf_full1$HI_Level = factor(nf_full1$HI_Level, 
                              levels=c("HI == 0",
                                       "HI < 0.25",
                                       "HI < 0.5",
                                       "HI < 0.75",
                                       "HI < 1",
                                       "HI == 1"))
nf_full1 %>%
  ggplot(aes(HI_Level, fill = HI_Level)) +
  geom_histogram(stat = "count")



positive_Samples$HI_Level = factor(positive_Samples$HI_Level, 
                                   levels=c("HI == 0",
                                            "HI < 0.25",
                                            "HI < 0.5",
                                            "HI < 0.75",
                                            "HI < 1",
                                            "HI == 1"))
positive_Samples %>%
  ggplot(aes(HI_Level, fill = HI_Level)) +
  geom_histogram(stat = "count")


positive_Samples %>%
  ggplot(aes(OP_10_x)) +
  geom_histogram(stat = "count")



#write.csv(Samples_2019_sel, "Samples_2019.csv")    
###############################################################################
#### FULL_GEN_HI DATASET ####
full_Gen_HI <- read.csv("alles.csv")

map <- full_Gen_HI %>%
  leaflet() %>%
  addProviderTiles("OpenStreetMap.DE") %>%
  setView(lat = 52.520007, lng =13.404954, zoom = 7)
map

# objects & layers full_Gen_HI
Year_2016 <- full_Gen_HI %>%
  filter(Year == "2016")
Year_2017 <- full_Gen_HI %>%
  filter(Year == "2017")
Year_2018 <- full_Gen_HI %>%
  filter(Year == "2018")
Year_2019 <- full_Gen_HI %>%
  filter(Year == "2019")
HI_0 <- full_Gen_HI %>%
  filter(HI == 0)
HI_bl_0.25 <- full_Gen_HI %>%
  filter(HI > 0, HI <= 0.25)
HI_bl_0.5 <- full_Gen_HI %>%
  filter(HI > 0.25, HI <= 0.5)
HI_bl_0.75 <- full_Gen_HI %>%
  filter(HI > 0.5, HI <= 0.75)
HI_bl_1 <- full_Gen_HI %>%
  filter(HI > 0.75, HI < 1)
HI_1 <- full_Gen_HI %>%
  filter(HI == 1)

# map YEARS full_Gen_HI
map %>%
  addCircleMarkers(data = Year_2016,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FFD919",
                   group = "2016") %>%
  addCircleMarkers(data = Year_2017,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FFB319",
                   group = "2017") %>%
  addCircleMarkers(data = Year_2018,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF8C19",
                   group = "2018") %>%
  addCircleMarkers(data = Year_2019,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF6619",
                   group = "2019") %>%
  addCircleMarkers(data = positive_Samples,
                   radius = 5,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#19FFFF",
                   group = "Crypto-Positive Samples",
                   opacity = .3) %>%
  addLayersControl(overlayGroups = c("2016","2017","2018","2019", "Crypto-Positive Samples")) %>%
  addLegend(colors = c("#FFD919", "#FFB319", "#FF8C19", "#FF6619", "#19FFFF"),
            labels = c("2016","2017","2018","2019", "Crypto-Positive Samples"),
            opacity = 1, 
            position = "bottomleft")

# map HI  full_Gen_HI

map <- full_Gen_HI %>%
  leaflet() %>%
  addProviderTiles("OpenStreetMap.DE") %>%
  setView(lat = 52.520007, lng =13.404954, zoom = 6)

map %>%
  #addCircleMarkers(data = full_Gen_HI,
  #                 radius = 2,
  #                 label = ~htmlEscape(Mouse_ID),
  #                 color = "yellow",
  #                 group = "full_Gen_HI",
  #                 clusterOptions = markerClusterOptions(showCoverageOnHover = T)) %>%
  addCircleMarkers(data = HI_0,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#1919FF",
                   group = "HI_0") %>%
  addCircleMarkers(data = HI_bl_0.25,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#4019FF",
                   group = "HI_bl_0.25") %>%
  addCircleMarkers(data = HI_bl_0.5,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#6619FF",
                   group = "HI_bl_0.5") %>%
  addCircleMarkers(data = HI_bl_0.75,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#B319FF",
                   group = "HI_bl_0.75") %>%
  addCircleMarkers(data = HI_bl_1,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF1966",
                   group = "HI_bl_1") %>%
  addCircleMarkers(data = HI_1,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF1919",
                   group = "HI_1") %>%
  addLayersControl(overlayGroups = c(#"full_Gen_HI",
                                     "HI_0", 
                                     "HI_bl_0.25", 
                                     "HI_bl_0.5", 
                                     "HI_bl_0.75", 
                                     "HI_bl_1", 
                                     "HI_1"
  )) %>%
  addLegend(colors = c("#1919FF", 
                       "#4019FF", 
                       "#6619FF", 
                       "#B319FF", 
                       "#FF1966", 
                       "#FF1919"
  ),
  labels = c("HI = 0", 
             "HI_bl_0.25", 
             "HI_bl_0.5", 
             "HI_bl_0.75", 
             "HI_bl_1", 
             "HI = 1"),
  opacity = 1, 
  position = "bottomleft")

###############################################################################
#### nf_FULL1 DATASET, contains filtered Samples, (Crypto-positive) ###########

nf_full1        <- read.csv("nf_full1_Corrected.csv")

map <- nf_full1 %>%
  leaflet() %>%
  addProviderTiles("OpenStreetMap.DE") %>%
  setView(lat = 52.520007, lng =13.404954, zoom = 7)
map

# objects and layers
Year_2016 <- nf_full1 %>%
  filter(Year == "2016")
Year_2017 <- nf_full1 %>%
  filter(Year == "2017")
Year_2018 <- nf_full1 %>%
  filter(Year == "2018")
Year_2019 <- nf_full1 %>%
  filter(Year == "2019")
positive_Samples <- nf_full1 %>%
  filter(Ct_mean > 0)
Oocyst_0 <- nf_full1 %>%
  filter(Oocyst_Predict == 0)
Oocyst_10_2 <- nf_full1 %>%
  filter(Oocyst_Predict > 0) %>%
  filter(Oocyst_Predict <= 10^2)
Oocyst_10_3 <- nf_full1 %>%
  filter(Oocyst_Predict > 10^2) %>%
  filter(Oocyst_Predict <= 10^3)
Oocyst_10_4 <- nf_full1 %>%
  filter(Oocyst_Predict > 10^3) %>%
  filter(Oocyst_Predict <= 10^4)
Oocyst_10_5 <- nf_full1 %>%
  filter(Oocyst_Predict > 10^4) %>%
  filter(Oocyst_Predict <= 10^6)

HI_0 <- nf_full1 %>%
  filter(HI == 0)
HI_bl_0.25 <- nf_full1 %>%
  filter(HI > 0, HI <= 0.25)
HI_bl_0.5 <- nf_full1 %>%
  filter(HI > 0.25, HI <= 0.5)
HI_bl_0.75 <- nf_full1 %>%
  filter(HI > 0.5, HI <= 0.75)
HI_bl_1 <- nf_full1 %>%
  filter(HI > 0.75, HI < 1)
HI_1 <- nf_full1 %>%
  filter(HI == 1)


# mapping Samples per Year --> better ALL Samples
map %>%
  addCircleMarkers(data = Year_2016,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FFD919",
                   group = "2016") %>%
  addCircleMarkers(data = Year_2017,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FFB319",
                   group = "2017") %>%
  addCircleMarkers(data = Year_2018,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF8C19",
                   group = "2018") %>%
  addCircleMarkers(data = Year_2019,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF6619",
                   group = "2019") %>%
  addCircleMarkers(data = positive_Samples,
                   radius = 5,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF1940",
                   group = "Crypto-Positive Samples",
                   opacity = .3) %>%
  addLayersControl(overlayGroups = c("2016","2017","2018","2019", "Crypto-Positive Samples")) %>%
  addLegend(colors = c("#FFD919", "#FFB319", "#FF8C19", "#FF6619", "#FF1940"),
            labels = c("2016","2017","2018","2019", "Crypto-Positive Samples"),
            opacity = 1, 
            position = "bottomleft")



# mapping HI
map <- nf_full1 %>%
  leaflet() %>%
  addProviderTiles("CartoDB") %>%
  setView(lat = 52.520007, lng =13.404954, zoom = 6)
map

map %>%
  addCircleMarkers(data = nf_full1,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "yellow",
                   group = "nf_full1",
                   clusterOptions = markerClusterOptions(showCoverageOnHover = T)) %>%
  addCircleMarkers(data = HI_0,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#1919FF",
                   group = "HI_0") %>%
  addCircleMarkers(data = HI_bl_0.25,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#4019FF",
                   group = "HI_bl_0.25") %>%
  addCircleMarkers(data = HI_bl_0.5,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#6619FF",
                   group = "HI_bl_0.5") %>%
  addCircleMarkers(data = HI_bl_0.75,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#B319FF",
                   group = "HI_bl_0.75") %>%
  addCircleMarkers(data = HI_bl_1,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF1966",
                   group = "HI_bl_1") %>%
  addCircleMarkers(data = HI_1,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF1919",
                   group = "HI_1") %>%
  addLayersControl(overlayGroups = c("nf_full1",
                                     "HI_0", 
                                     "HI_bl_0.25", 
                                     "HI_bl_0.5", 
                                     "HI_bl_0.75", 
                                     "HI_bl_1", 
                                     "HI_1"
  )) %>%
  addLegend(colors = c("#1919FF", 
                       "#4019FF", 
                       "#6619FF", 
                       "#B319FF", 
                       "#FF1966", 
                       "#FF1919"
  ),
  labels = c("HI = 0", 
             "HI_bl_0.25", 
             "HI_bl_0.5", 
             "HI_bl_0.75", 
             "HI_bl_1", 
             "HI = 1"),
  opacity = 1, 
  position = "bottomleft")

# mapping Oocyst predictions
map %>%
  addCircleMarkers(data = Oocyst_0,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FFD919", 
                   opacity = 0.1,
                   group = "Oocyst_0") %>%
  addCircleMarkers(data = Oocyst_10_2,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF8C19",
                   group = "Oocyst_10_2") %>%
  addCircleMarkers(data = Oocyst_10_3,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF6619",
                   group = "Oocyst_10_3") %>%
  addCircleMarkers(data = Oocyst_10_4,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF4019",
                   group = "Oocyst_10_4") %>%
  addCircleMarkers(data = Oocyst_10_5,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF1919",
                   group = "Oocyst_10_5") %>%
  addLayersControl(overlayGroups = c("Oocyst_0", 
                                     "Oocyst_10_2", 
                                     "Oocyst_10_3", 
                                     "Oocyst_10_4", 
                                     "Oocyst_10_5")) %>%
  addLegend(colors = c("#FFD919", 
                       "#FFB319", 
                       "#FF8C19", 
                       "#FF6619", 
                       "#FF4019"), 
                       #"#FF1919"),
            labels = c("negative",
                       "< 100 Oocysts",
                       "< 1000 Oocysts",
                       "< 10.000 Oocysts", 
                       "< 100.000 Oocysts"), 
                       #"< 1.000.000 Oocysts"),
            opacity = 1, 
            position = "bottomleft")


################################################################################
#### Samples_2019 ####

Samples_2019 <- nf_full1 %>%
  filter(Year == "2019") %>%
  arrange(Mouse_ID) %>%
  select(Mouse_ID, Year, Ct_mean, HI, HI_Level, 
         Transect, Latitude, Longitude, Sex, 
         pos_ratio, overall_ratio, Oocyst_Predict)

Samples_2019 %>%
  filter(Ct_mean > 0)

Samples_2019_sel <- Samples_2019 %>%
  select(Mouse_ID, Sex, Latitude, Longitude, Ct_mean, HI, Oocyst_Predict)

################################################################################
#### Crypto-positive Samples ####
nf_full1 <- read.csv("nf_full1_Corrected.csv")
positive_Samples <- nf_full1 %>%
  filter(Ct_mean > 0)
write.csv(positive_Samples, "positive_Samples.csv")
positive_Samples <- read.csv("positive_Samples.csv")

positive_Samples$HI_Level = factor(positive_Samples$HI_Level, 
                                   levels=c("HI == 0",
                                               "HI < 0.25",
                                               "HI < 0.5",
                                               "HI < 0.75",
                                               "HI < 1",
                                               "HI == 1"))

positive_Samples$Oocyst_Level = factor(positive_Samples$Oocyst_Level, 
                                   levels=c("OP = 0",
                                            "OP < 10^2",
                                            "OP < 10^3",
                                            "OP < 10^4",
                                            "OP < 10^5",
                                            "OP < 10^6"))

positive_Samples %>%
  ggplot(aes(HI_Level, fill = HI_Level)) +
  geom_bar(stat = "count")

positive_Samples %>%
  ggplot(aes(Oocyst_Level, fill = Oocyst_Level)) +
  geom_bar(stat = "count")
  


HI_0 <- positive_Samples %>%
  filter(HI == 0)
HI_bl_0.25 <- positive_Samples %>%
  filter(HI > 0, HI <= 0.25)
HI_bl_0.5 <- positive_Samples %>%
  filter(HI > 0.25, HI <= 0.5)
HI_bl_0.75 <- positive_Samples %>%
  filter(HI > 0.5, HI <= 0.75)
HI_bl_1 <- positive_Samples %>%
  filter(HI > 0.75, HI < 1)
HI_1 <- positive_Samples %>%
  filter(HI == 1)

Oocyst_0 <- positive_Samples %>%
  filter(Oocyst_Predict == 0)
Oocyst_10_2 <- positive_Samples %>%
  filter(Oocyst_Predict > 0) %>%
  filter(Oocyst_Predict <= 10^2)
Oocyst_10_3 <- positive_Samples %>%
  filter(Oocyst_Predict > 10^2) %>%
  filter(Oocyst_Predict <= 10^3)
Oocyst_10_4 <- positive_Samples %>%
  filter(Oocyst_Predict > 10^3) %>%
  filter(Oocyst_Predict <= 10^4)
Oocyst_10_5 <- positive_Samples %>%
  filter(Oocyst_Predict > 10^4) %>%
  filter(Oocyst_Predict <= 10^5)
Oocyst_10_6 <- positive_Samples %>%
  filter(Oocyst_Predict > 10^5) %>%
  filter(Oocyst_Predict <= 10^6)

map %>%
  addCircleMarkers(data = HI_0,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#1919FF",
                   group = "HI_0") %>%
  addCircleMarkers(data = HI_bl_0.25,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#4019FF",
                   group = "HI_bl_0.25") %>%
  addCircleMarkers(data = HI_bl_0.5,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#6619FF",
                   group = "HI_bl_0.5") %>%
  addCircleMarkers(data = HI_bl_0.75,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#B319FF",
                   group = "HI_bl_0.75") %>%
  addCircleMarkers(data = HI_bl_1,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF1966",
                   group = "HI_bl_1") %>%
  addCircleMarkers(data = HI_1,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF1919",
                   group = "HI_1") %>%
  addLayersControl(overlayGroups = c("HI_0", 
                                     "HI_bl_0.25", 
                                     "HI_bl_0.5", 
                                     "HI_bl_0.75", 
                                     "HI_bl_1", 
                                     "HI_1"
  )) %>%
  addLegend(colors = c("#1919FF", 
                       "#4019FF", 
                       "#6619FF", 
                       "#B319FF", 
                       "#FF1966", 
                       "#FF1919"
  ),
  labels = c("HI = 0", 
             "HI_bl_0.25", 
             "HI_bl_0.5", 
             "HI_bl_0.75", 
             "HI_bl_1", 
             "HI = 1"),
  opacity = 1, 
  position = "bottomleft")


map %>%
  addCircleMarkers(data = Oocyst_10_2,
                   radius = 3,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF8095",
                   group = "Oocyst_10_2") %>%
  addCircleMarkers(data = Oocyst_10_3,
                   radius = 3,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#95B300",
                   group = "Oocyst_10_3") %>%
  addCircleMarkers(data = Oocyst_10_4,
                   radius = 3,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#00CC66",
                   group = "Oocyst_10_4") %>%
  addCircleMarkers(data = Oocyst_10_5,
                   radius = 3,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#198CFF",
                   group = "Oocyst_10_5") %>%
  addCircleMarkers(data = Oocyst_10_6,
                   radius = 3,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#BF00E6",
                   group = "Oocyst_10_6") %>%
  addLayersControl(overlayGroups = c("Oocyst_10_2", 
                                     "Oocyst_10_3", 
                                     "Oocyst_10_4", 
                                     "Oocyst_10_5",
                                     "Oocyst_10_6")) %>%
  addLegend(colors = c("#FF8095", 
                       "#95B300", 
                       "#00CC66", 
                       "#198CFF",
                       "#BF00E6"), 
            #"#FF1919"),
            labels = c("< 100 Oocysts",
                       "< 1000 Oocysts",
                       "< 10.000 Oocysts", 
                       "< 100.000 Oocysts",
                       "< 1.000.000 Oocysts"), 
            #"< 1.000.000 Oocysts"),
            opacity = 1, 
            position = "bottomleft")



###############################################################################
#### HI Analysis ####

full_Gen_HI %>%
  ggplot(aes(HI, Mouse_ID, col = HI)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = " red") +
  theme(axis.text.y = element_blank())

nf_full1 %>%
  ggplot(aes(HI, Mouse_ID, col = HI)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = " red") +
  theme(axis.text.y = element_blank())

positive_Samples %>%
  ggplot(aes(HI, Mouse_ID, col = HI)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = " red") +
  theme(axis.text.y = element_blank())


#### Oocyst Analysis ####
positive_Samples %>%
  ggplot(aes(Mouse_ID, Ct_mean, col = Ct_mean)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = " red") +
  theme(axis.text.x = element_blank())

Plate6_Candidates <- read.csv("Plate6_Candidates.csv")
Plate6_Candidates %>%
  ggplot(aes(Mouse_ID, Ct_mean, col = Ct_mean)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = " red") +
  theme(axis.text.x = element_blank())

Pos_Plate6 <- left_join(Plate6_Candidates, positive_Samples, by = "Mouse_ID")
Pos_Plate6 %>%
  ggplot(aes(x = Mouse_ID)) +
  geom_point(aes(y =  Ct_mean.x), shape = 2) +
  geom_point(aes(y =  Ct_mean.y)) +
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40),
                     limits = c(0, 40)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### old Oocyst Predictions ####
nf_full1 <- read.csv("nf_full1_Corrected.csv")
all_Samples_old <- read.csv("all_Samples_old.csv")
all_Samples_old <- all_Samples_old %>%
  mutate(Mouse_ID = Sample.Name) %>%
  mutate(Oocyst_Predict = Oocyst_Predict) #%>%
  select(Mouse_ID, Ct_mean, Oocyst_Predict)
all_Samples_old_join <- left_join(all_Samples_old, nf_full1, by = "Mouse_ID")


ABI_Best_thSC     <- read.csv("ABI_Best_SC.csv")
f_ABI_Best_thSC   <-  filter(ABI_Best_thSC, Ct_mean > 0)
linear_model0     <- lm(log2(Amount_Oocysts) ~ Ct_mean, data = f_ABI_Best_thSC)
Oocyst_Predict    <- 2^predict(linear_model0, newdata = all_Samples_old)


all_Samples_old <- data.frame(all_Samples_old, Oocyst_Predict)
all_Samples_old <- all_Samples_old %>%
  mutate(Oocyst_Predict = replace(Oocyst_Predict, Oocyst_Predict == '4292821751815.77', '0')) %>%
  select(Mouse_ID, Ct_1, Ct_2, Ct_3, Ct_mean, Oocyst_1, Oocyst_2, Oocyst_mean, Oocyst_calculated, Fnal_Oocyst_calculated, Oocyst_Predict)
write.csv(all_Samples_old, "all_Samples_old1.csv")
write.csv(all_Samples_old_join, "all_Samples_old_join.csv")
all_Samples_old_join <- read.csv("all_Samples_old_join.csv")

# variables ####
Oocyst_0 <- all_Samples_old_join %>%
  filter(Oocyst_Predict.x == 0)
Oocyst_10_2 <- all_Samples_old_join %>%
  filter(Oocyst_Predict.x > 0) %>%
  filter(Oocyst_Predict.x <= 10^2)
Oocyst_10_3 <- all_Samples_old_join %>%
  filter(Oocyst_Predict.x > 10^2) %>%
  filter(Oocyst_Predict.x <= 10^3)
Oocyst_10_4 <- all_Samples_old_join %>%
  filter(Oocyst_Predict.x > 10^3) %>%
  filter(Oocyst_Predict.x <= 10^4)
Oocyst_10_5 <- all_Samples_old_join %>%
  filter(Oocyst_Predict.x > 10^4) %>%
  filter(Oocyst_Predict.x <= 10^5)
Oocyst_10_6 <- all_Samples_old_join %>%
  filter(Oocyst_Predict.x > 10^5) %>%
  filter(Oocyst_Predict.x <= 10^6)
Oocyst_10_7 <- all_Samples_old_join %>%
  filter(Oocyst_Predict.x > 10^6)


map <- nf_full1 %>%
  leaflet() %>%
  addProviderTiles("OpenStreetMap.DE") %>%
  setView(lat = 52.520007, lng =13.404954, zoom = 6)



# mapping Oocyst predictions ####
map %>%
  addCircleMarkers(data = Oocyst_0,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FFD919", 
                   opacity = 0.1,
                   group = "Oocyst_0") %>%
  addCircleMarkers(data = Oocyst_10_2,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF8C19",
                   group = "Oocyst_10_2") %>%
  addCircleMarkers(data = Oocyst_10_3,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF6619",
                   group = "Oocyst_10_3") %>%
  addCircleMarkers(data = Oocyst_10_4,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF4019",
                   group = "Oocyst_10_4") %>%
  addCircleMarkers(data = Oocyst_10_5,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF1919",
                   group = "Oocyst_10_5") %>%
  addCircleMarkers(data = Oocyst_10_6,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "green",
                   group = "Oocyst_10_6") %>%
  addCircleMarkers(data = Oocyst_10_7,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "blue",
                   group = "Oocyst_10_7") %>%
  addLayersControl(overlayGroups = c("Oocyst_0", 
                                     "Oocyst_10_2", 
                                     "Oocyst_10_3", 
                                     "Oocyst_10_4", 
                                     "Oocyst_10_5",
                                     "Oocyst_10_6",
                                     "Oocyst_10_7"
                                     )) %>%
  addLegend(colors = c("#FFD919", 
                       "#FFB319", 
                       "#FF8C19", 
                       "#FF6619", 
                       "#FF4019",
                       "green",
                       "blue"
                       ), 
            #"#FF1919"),
            labels = c("negative",
                       "< 100 Oocysts",
                       "< 1000 Oocysts",
                       "< 10.000 Oocysts", 
                       "< 100.000 Oocysts",
                       "< 1 Mio Oocysts",
                       "< 10 Mio Oocyts"), 
            #"< 1.000.000 Oocysts"),
            opacity = 1, 
            position = "bottomleft")


#### TRUSTING IN DATA ####
all_Samples_old <- read.csv("all_Samples_old1.csv")
nf_full1 <- read.csv("nf_full1_Corrected.csv")

all_Samples_DF1 <- all_Samples_old %>%
  filter(X <= 515) %>%
  mutate(Delta_OP_mean_vs_calc = sqrt((Oocyst_mean - Oocyst_calculated)^2)) %>%
  mutate(Ct0 = Ct_1 == 0 | Ct_2 == 0)
  

all_Samples_DF2 <- all_Samples_old %>%
  filter(X > 515) %>%
  mutate(Delta_OP_mean_vs_calc = sqrt((Oocyst_mean - Oocyst_calculated)^2)) %>%
  mutate(Ct0 = Ct_1 == 0 | Ct_2 == 0 | Ct_3 == 0)
 
all_Samples_old <- union(all_Samples_DF1, all_Samples_DF2)
all_Samples_old
write.csv(all_Samples_old, "all_Samples_old1.csv")

# CT_FRENZY DATA ####
Ct_frenzy <- left_join(all_Samples_old, nf_full1, by = "Mouse_ID")

Ct_frenzy <- Ct_frenzy %>%
  mutate(Ct_Diff = sqrt((Ct_mean.x - Ct_mean.y)^2),
         Delta_Ct = sqrt((Ct_1 - Ct_2)^2),
         old_goofy = Delta_Ct > 3,
         new_goofy = Ct_Diff > 3,
         Sex = replace(Sex, Sex == 'female', 'F'),
         Sex = replace(Sex, Sex == 'male', 'M'),
         Consistent = Ct_1 + Ct_2)

write.csv(Ct_frenzy, "Ct_frenzy.csv")

Ct_frenzy <- Ct_frenzy %>%
  select(Mouse_ID, Ct_1, Ct_2, Ct_3, Delta_Ct, 
         Ct_mean.x, Oocyst_Predict.x,
         Ct_mean.y, Oocyst_Predict.y, 
         Ct_Diff, Ct0,
         Oocyst_1, Oocyst_2, Oocyst_mean, Oocyst_calculated, 
         old_goofy, new_goofy, 
         Sex, HI, HI_Level, Year,
         Transect, Latitude, Longitude,
         Tested, Double_Tested)

count(Ct_frenzy, old_goofy)
count(Ct_frenzy, new_goofy)


Ct_frenzy %>%
  filter(Double_Tested == T) %>%
  ggplot(aes(Ct_mean.x, Ct_mean.y, label = Mouse_ID)) +
  geom_point() +
  geom_label(nudge_y = 1.5, alpha = 0.5) +
  annotate(geom = "line",
           x = c(0,30),
           y = c(0,30))


# Ct = 0 excluded
# if Year not known, it wasn't included in cleaned data
CryptoAll_Pos <- read.csv("CryptoAll_Pos.csv")
CryptoAll_Pos <- CryptoAll_Pos %>%
  #rename(Ct_mean = Ct) %>%
  distinct()


  mutate(Oocyst_Predict = replace(Oocyst_Predict, Oocyst_Predict == '4292821751815.77', '0')) %>%
  
Crypto_nf_full1 <- left_join(CryptoAll_Pos_u, nf_full1, by = "Mouse_ID")
  
setdiff(CryptoAll_Pos, CryptoAll_Pos_u)


Ct_frenzy

Ct_frenzy %>%
  filter(Ct_mean.y != 0) %>%
  #filter(Ct_mean.y > 100) %>%
  #ggplot(aes(Mouse_ID, Oocyst_Predict.x, col = Oocyst_Predict.x)) +
  ggplot(aes(Mouse_ID, Oocyst_Predict.y, col = Oocyst_Predict.y)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_blank())


Ct_frenzy %>%
  ggplot(aes(Ct_mean.x, Ct_mean.y)) +
  geom_point()

Ct_frenzy %>%
  filter(Ct0 == F) %>%
  filter(!is.na(Year)) %>%
  #filter(Ct_mean.y != 0) %>%
  ggplot(aes(x = Mouse_ID, col = Ct_Diff)) +
  geom_point(aes(y = Ct_mean.x), shape = 1) + # see through circles = old qPCRs
  geom_point(aes(y = Ct_mean.y)) + # filled circles = new qPCRs
  scale_color_gradient(low = "blue", # highly similar Ct values = blue
                       high = "red") + # highly divergent Ct values = red
  theme(axis.text.x = element_blank()) +
  labs(y = "Ct values") +
  ggtitle("Comparison of qPCR results ", 
          subtitle = "old   qPCRs   =   see-through circles 
new qPCRs   =   filled circles" )


  
Ct_frenzy %>%  
  #filter(Ct0 == F) %>%
  filter(!is.na(Year)) %>%
  #ggplot(aes(Ct_mean.x, Ct_mean.y, col = Oocyst_Predict.x)) + # shows extremely high predictions
  ggplot(aes(Ct_mean.x, Ct_mean.y, col = Oocyst_Predict.y)) + # shows more reasonable prediction
  geom_point() +
  labs(x = "old Ct values",
       y = "new Ct values") +
  ggtitle("Difference of Ct_means") +
  scale_color_gradient(low = "blue",
                       high = "red")

count(all_Samples_old, Delta_OP_mean_vs_calc)
Ct_frenzy %>%  
  filter(Ct0 == F) %>%
  filter(!is.na(Year)) %>%
  #filter(Oocyst_calculated != 0) %>%
  #filter(Ct_mean.x != 0) %>%
  #filter(Ct_mean.y != 0) %>%
  count()
  

  
  # create col whether it's actually infected
  # if no oocysts & no Ct --> prob F
  # check model prediction with oocyst prediction
  # added data should improve accuracy of the model
