# Load Libraries
library(ggplot2)
library(dplyr)
library(ggeffects)
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

Data_Initial <- read.csv("nf_full1_Corrected.csv")
glimpse(Data_Initial)

Data_work <- Data_Initial


# 77 Samples were measured at twice, at least once by me
# Consistent Data applies to Ct values that are roughly (Ct +- 3) similar

Data_work %>%  
  ggplot(aes(x = Ct_mean_1_ABI, col = Machine, shape = Tested_by)) +
  geom_point(aes(y = Ct_mean_Ep)) +
  geom_point(aes(y = Ct_mean_2_ABI)) +
  ggtitle("Difference of Ct-Values in double-tested Samples") +
  #labs(y = "Ct") +
  scale_color_brewer(palette = "Paired")

# mean Ep values
Data_work %>%
ggplot(aes(x = Ct_mean_1_ABI, col = Machine, shape = Tested_by)) +
  geom_point(aes(y = Ct_mean_Ep)) +
  geom_point(aes(y = Ct_mean_2_ABI)) +
  scale_y_continuous(name = "Eppendorf Ct",
                     sec.axis =  sec_axis(~., name = "old ABI Ct")) +
  scale_color_brewer(palette = "Paired") +
  ggtitle("Difference of Ct-Values in double-tested Samples") +
  theme(
    axis.title.y = element_text(color = "#3399FF", size=13),
    axis.title.y.right = element_text(color = "#00CC22", size=13)
  )
 
# Consistent Ep values 
Data_work %>%
  ggplot(aes(x = Ct_mean_1_ABI, col = Machine, shape = Tested_by)) +
  geom_point(aes(y = Ct_Ep_Consistent)) +
  geom_point(aes(y = Ct_mean_2_ABI)) +
  scale_y_continuous(name = "Eppendorf Ct",
                     sec.axis =  sec_axis(~., name = "old ABI Ct")) +
  scale_color_brewer(palette = "Paired") +
  ggtitle("Difference of Ct-Values in double-tested Samples") +
  theme(
    axis.title.y = element_text(color = "#3399FF", size=13),
    axis.title.y.right = element_text(color = "#00CC22", size=13)
  )



# all mean values
Data_work %>%
  ggplot(aes(x = Ct_mean, col = Machine, shape = Tested_by)) +
  geom_point(aes(y = Ct_mean_Ep)) +
  geom_point(aes(y = Ct_mean_ABI)) +
  scale_y_continuous(name = "Eppendorf Ct",
                     sec.axis =  sec_axis(~., name = "old ABI Ct")) +
  scale_color_brewer(palette = "Paired") +
  theme(
    axis.title.y = element_text(color = "#3399FF", size=13),
    axis.title.y.right = element_text(color = "#00CC22", size=13)
  )


Data_work %>%  
  ggplot(aes(x = Ct_mean, col = Machine, shape = Tested_by, label = Tested_by)) +
  geom_point(aes(y = Ct_mean_Ep)) +
  geom_point(aes(y = Ct_mean_ABI)) +
  #geom_label(aes(y = Ct_mean_Ep), alpha = 0.5) +
  ggtitle("Difference of Ct-Values in double-tested Samples") +
  labs(y = "Ct") +
  scale_color_brewer(palette = "Paired")
  
#Data_work <- Data_work %>%
#  mutate(Ct_Diff_Ep = Ct_mean_1_ABI - Ct_mean_Ep,
#         Ct_Diff_ABI = Ct_mean_2_ABI - Ct_mean_1_ABI)

Data_work %>%
  filter(Plate == "7") %>%
  ggplot(aes(Mouse_ID, col = Consistent)) +
  geom_bar(stat = "identity", aes(y = -(Ct_mean_Ep))) +
  geom_bar(stat = "identity", aes(y = Ct_mean_1_ABI)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Difference of Ct-Values in double-tested Samples")

Data_work %>%
  filter(Plate == "7") %>%
  ggplot(aes(Mouse_ID, col = Consistent)) +
  geom_bar(stat = "identity", aes(y = -(Ct_Ep_Consistent))) +
  geom_bar(stat = "identity", aes(y = Ct_mean_1_ABI)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Difference of Ct-Values in double-tested Samples")

# predict Oocysts for positive Samples ####
ABI_Best_thSC     <- read.csv("ABI_Best_SC.csv")
f_ABI_Best_thSC   <-  filter(ABI_Best_thSC, Ct_mean > 0)
linear_model0     <- lm(log2(Amount_Oocysts) ~ Ct_mean, data = f_ABI_Best_thSC)
Oocyst_Predict    <- 2^predict(linear_model0, newdata = Data_work)


Data_work <- data.frame(Data_work, Oocyst_Predict)
Data_work <- Data_work %>%
  mutate(Oocyst_Predict = replace(Oocyst_Predict, Oocyst_Predict == '4292821751815.77', '0'))



write.csv(Data_work, "Data.csv")
Data_work <- read.csv("Data.csv")
Data_work

Data_Initial <- read.csv("nf_full1_Corrected.csv")

Data_work <- Data_Initial %>%
  select(Mouse_ID,
         Ct_1_Ep, Ct_2_Ep, Ct_Ep_Consistent,
         Ct_1_ABI, Ct_2_ABI, Ct_3_ABI,
         Ct_4_ABI, Ct_5_ABI, Ct_6_ABI, 
         Ct_7_ABI, Ct_8_ABI, Ct_9_ABI,
         Ct_mean_1_ABI, Ct_mean_2_ABI, Ct_mean_3_ABI,
         Ct_mean_Ep, Ct_mean_ABI, Ct_mean, 
         StDev, Measurements,
         Tested_Me, Consistent,
         HI, HI_Level, Sex,
         Transect, Latitude, Longitude,
         Year, Date, Plate, Machine, Tested_by,
         Oocyst_Predict
  )


Data_work %>%
  filter(Ct_mean > 0) %>%
  #filter(Mouse_ID != c("AA_0322")) %>%
  # log10 might improve visuals
  ggplot(aes(Mouse_ID, log10(Oocyst_Predict), col = Tested_by, shape = Machine)) +
  geom_point() +
  theme(axis.text.x = element_blank())
  #scale_color_brewer(palette = "Paired") +
  #ggtitle("Oocyst Predictions excluding AA_0322")

# predict Oocysts, only calculating with old, wrong Ct-values
#pure_Ep           <- read.csv("pure_old_Ep_Ct_for_OP.csv")
#pure_Ep <- pure_Ep %>%
#  select(Mouse_ID, Tested_Me, Ct_mean)
#Oocyst_Predict    <- 2^predict(linear_model0, newdata = pure_Ep)
#pure_Ep <- data.frame(pure_Ep, Oocyst_Predict)

#pure_Ep <- pure_Ep %>%
#  mutate(Oocyst_Predict = replace(Oocyst_Predict, Oocyst_Predict == '4292821751815.77', '0'))


#write.csv(pure_Ep, "pure_Ep.csv")
#pure_Ep <- read.csv("pure_Ep.csv")
  
#pure_Ep %>%
#  ggplot(aes(Mouse_ID, log10(Oocyst_Predict), col = Tested_Me, label = Mouse_ID)) +
#  geom_point() +
  #geom_label() +
#  theme(axis.text.x = element_blank()) +
  #scale_color_brewer(palette = "Paired") +
#  ggtitle("Oocyst Predictions without Data Revision")

# use HI Loci instead of HI_Level ####

#full_Gen <- read.csv("full_Gen.csv")
#full_Gen_unique <- full_Gen %>%
#  distinct(Mouse_ID, .keep_all = TRUE)
#full_Gen_unique %>%
#  count(Mouse_ID) %>%
#  filter(n > 1)

#nf_full1        <- read.csv("nf_full1_Corrected.csv")
#Emanuel_Data    <- read.csv("EmanuelData.csv") %>%
#  mutate(Mouse_ID = PIN)
#full_Data <- left_join(full_Gen_unique, Emanuel_Data, by = "Mouse_ID", "Transect")
#full_Data <- left_join(full_Data, nf_full1, by = "Mouse_ID")
#write.csv(full_Data, "full_Data.csv")

#full_Data <- read.csv("full_Data_Corrected.csv")
#full_Data <- left_join(full_Data, nf_full1)


full_Data <- full_Data %>%
  mutate(HI_State = ifelse(HI == 1, "pure Mmm", 
                           ifelse(HI == 0, "pure Mmd",
                                  ifelse(HI_NLoci > 0 & HI_NLoci <= 0.5, "Western Hybrid",
                                         ifelse(HI_NLoci < 1 & HI_NLoci > 0.5, "Eastern Hybrid",
                                                ifelse(NA)))))) %>%
  mutate(Oocyst_Level = ifelse(Oocyst_Predict == 0, "OP = 0",
                               ifelse(Oocyst_Predict < 100, "OP < 100",
                                      ifelse(Oocyst_Predict < 10^3, "OP < 1000",
                                             ifelse(Oocyst_Predict < 10^4, "OP < 10.000",
                                                    ifelse(Oocyst_Predict < 10^5, "OP < 100.000",
                                                           ifelse(Oocyst_Predict < 10^6, "OP < 1.000.000")
                                                    ))))))
full_Data$Location <- paste(full_Data$Latitude, full_Data$Longitude, sep = ",")

# introduce mus_caught and mus_Crypto_caught 
full_Data_long_pull <- full_Data %>%
  count(Longitude)%>% mutate(mus_caught = n) %>% select(Longitude, mus_caught) %>% arrange(mus_caught)
  
full_Data_long_pull_Crypto <- full_Data %>%
  filter(Ct_mean > 0) %>% count(Longitude) %>% mutate(Crypto_mus_caught = n) %>% select(Longitude, Crypto_mus_caught) %>% arrange(Crypto_mus_caught)

full_Data1 <- left_join(full_Data, full_Data_long_pull)
full_Data2 <- left_join(full_Data1, full_Data_long_pull_Crypto)
full_Data <- full_Data2
write.csv(full_Data, "full_Data_Corrected.csv")


  



##############################################################################
# mapping full_Data Dataset ####################################################

full_Data <- read.csv("full_Data_Corrected.csv")

  
map <- full_Data %>%
  leaflet() %>%
  addProviderTiles("OpenStreetMap.DE") %>%
  setView(lat = 52.520007, lng =13.404954, zoom = 7)
map

# objects and layers
Oocyst_0 <- full_Data %>%
  filter(Oocyst_Predict == 0)
Oocyst_10_2 <- full_Data %>%
  filter(Oocyst_Predict > 0) %>%
  filter(Oocyst_Predict <= 10^2)
Oocyst_10_3 <- full_Data %>%
  filter(Oocyst_Predict > 10^2) %>%
  filter(Oocyst_Predict <= 10^3)
Oocyst_10_4 <- full_Data %>%
  filter(Oocyst_Predict > 10^3) %>%
  filter(Oocyst_Predict <= 10^4)
Oocyst_10_5 <- full_Data %>%
  filter(Oocyst_Predict > 10^4) %>%
  filter(Oocyst_Predict <= 10^5)
Oocyst_10_6 <- full_Data %>%
  filter(Oocyst_Predict > 10^5) %>%
  filter(Oocyst_Predict <= 10^6)


interesting_for_revisit <- full_Data %>%
  filter(mus_caught > 1, Oocyst_Predict > 1000)

map %>%
  addCircleMarkers(data = interesting_for_revisit,
                   radius = 2, 
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF1940",
                   opacity = 0.3,
                   group = "interesting_for_revisit",
                   popup= ~paste("<b>Mouse_ID:<b>", as.character(Mouse_ID), "<br>",
                                 "<b>Location:<b>", as.character(Location), "<br>",
                                 "<b>mus caught:<b>", as.character(mus_caught), "<br>",
                                 "<b>Crypto positive:<b>", as.character(Crypto_mus_caught), "<br>",
                                 "<b>HI:<b>",      as.character(HI), "<br>",
                                 "<b>Year:<b>",    as.character(Year),"<br>",
                                 "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                 "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                 "<b>Sex:<b>", Sex, "<br>",
                                 sep=" ")) %>%
  addLayersControl(overlayGroups = c("interesting_for_revisit")) %>%
  addLegend(colors = "#FF1940",
            labels = c("interesting_for_revisit"),
            opacity = 1, 
            position = "bottomleft")


#### HI easy

full_Data_sel <- full_Data %>%
  filter(HI != "NA")
data_col = colorFactor(palette = 'Blues', full_Data_sel$HI)#puts it into a spatial data frame (negatives of those columns as well)

data_col_HI       = colorFactor(matlab.like2(14), full_Data_sel$HI)
data_col_HI_NLoci = colorFactor(matlab.like2(4), full_Data_sel$HI_NLoci)

data_col_Oocysts  = colorRampPalette(c("red", "purple", "blue"),
                                                    space = "Lab")


map %>%
  addCircleMarkers(data = full_Data_sel,  color = ~data_col_HI(HI), radius = 3, group = "HI") %>%
  #addCircleMarkers(data = full_Data_sel, coloer = ~data_col_HI()) %>%
  addLegend('bottomleft', 
            pal = data_col_HI_NLoci, 
            values = full_Data_sel$HI_NLoci,
            opacity = 1) %>%
  addLayersControl()
map





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
  addLayersControl(overlayGroups = c("2016","2017","2018","2019")) %>%
  addLegend(colors = c("#FFD919", "#FFB319", "#FF8C19", "#FF6619"),
            labels = c("2016","2017","2018","2019"),
            opacity = 1, 
            position = "bottomleft")

map %>%
  addCircleMarkers(data = Oocyst_0,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FFD919", 
                   opacity = 0.1,
                   group = "Oocyst_0",
                   popup= ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                 "<b>Location:<b>",as.character(Location), "<br>",
                                 "<b>mus caught:<b>",as.character(n), "<br>",
                                 "<b>HI:<b>",      as.character(HI), "<br>",
                                 "<b>Year:<b>",    as.character(Year),"<br>",
                                 "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                 "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                 "<b>Sex:<b>", Sex, "<br>",
                                 sep=" ")) %>%
  addCircleMarkers(data = Oocyst_10_2,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF8C19",
                   opacity = 0.1,
                   group = "Oocyst_10_2",
                   popup= ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                 "<b>Location:<b>",as.character(Location), "<br>",
                                 "<b>mus caught:<b>",as.character(n), "<br>",
                                 "<b>HI:<b>",      as.character(HI), "<br>",
                                 "<b>Year:<b>",    as.character(Year),"<br>",
                                 "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                 "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                 "<b>Sex:<b>", Sex, "<br>",
                                 sep=" ")) %>%
  addCircleMarkers(data = Oocyst_10_3,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF6619",
                   opacity = 0.1,
                   group = "Oocyst_10_3",
                   popup= ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                 "<b>Location:<b>",as.character(Location), "<br>",
                                 "<b>mus caught:<b>",as.character(n), "<br>",
                                 "<b>HI:<b>",      as.character(HI), "<br>",
                                 "<b>Year:<b>",    as.character(Year),"<br>",
                                 "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                 "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                 "<b>Sex:<b>", Sex, "<br>",
                                 sep=" ")) %>%
  addCircleMarkers(data = Oocyst_10_4,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF4019",
                   opacity = 0.3,
                   group = "Oocyst_10_4",
                   popup= ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                 "<b>Location:<b>",as.character(Location), "<br>",
                                 "<b>mus caught:<b>",as.character(n), "<br>",
                                 "<b>HI:<b>",      as.character(HI), "<br>",
                                 "<b>Year:<b>",    as.character(Year),"<br>",
                                 "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                 "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                 "<b>Sex:<b>", Sex, "<br>",
                                 sep=" ")) %>%
  addCircleMarkers(data = Oocyst_10_5,
                   radius = 2,
                   label = ~htmlEscape(Mouse_ID),
                   color = "#FF1940",
                   opacity = 0.5,
                   group = "Oocyst_10_5",
                   popup= ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                 "<b>Location:<b>",as.character(Location), "<br>",
                                 "<b>mus caught:<b>",as.character(n), "<br>",
                                 "<b>HI:<b>",      as.character(HI), "<br>",
                                 "<b>Year:<b>",    as.character(Year),"<br>",
                                 "<b>Ct:<b>",      as.character(Ct_mean),"<br>",
                                 "<b>Oocysts:<b>", as.character(Oocyst_Predict), "<br>",
                                 "<b>Sex:<b>", Sex, "<br>",
                                 sep=" ")) %>%
  addLayersControl(overlayGroups = c("Oocyst_0", 
                                     "Oocyst_10_2", 
                                     "Oocyst_10_3", 
                                     "Oocyst_10_4", 
                                     "Oocyst_10_5")) %>%
  addLegend(colors = c("#FFD919", 
                       "#FFB319", 
                       "#FF8C19", 
                       "#FF6619", 
                       "#FF1940"), 
            #"#FF1919"),
            labels = c("negative",
                       "< 100 Oocysts",
                       "< 1000 Oocysts",
                       "< 10.000 Oocysts", 
                       "< 100.000 Oocysts"), 
            #"< 1.000.000 Oocysts"),
            opacity = 1, 
            position = "bottomleft")





