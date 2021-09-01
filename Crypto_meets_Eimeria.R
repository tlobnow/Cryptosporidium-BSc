library(tidyverse)
library(skimr)
library(data.table)
library(visdat)
library(viridis)
library(RColorBrewer)

Crypto_Detection      <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Crypto_Detection.csv")
Eimeria_Detection     <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/a1d58d91589849ea4dbabbbf82e6c896ef7ec95c/data_products/Eimeria_Detection.csv")
HZ_CEWE_ELISA         <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/ef8575f03ff3a167fda850fd0e93c5de13843aa6/data_input/HZ19_CEWE_ELISA.csv")
Eimeria_Summary       <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Eimeria_detection/Summary_eimeria.csv")

# HZ19
HZ19_wild_immuno_long <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/ef8575f03ff3a167fda850fd0e93c5de13843aa6/data_input/HZ19_wild_immuno_long.csv")
HZ19_immuno_long      <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/ef8575f03ff3a167fda850fd0e93c5de13843aa6/data_input/HZ19_immuno_long.csv")
HZ19_MES_FACS         <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/ef8575f03ff3a167fda850fd0e93c5de13843aa6/data_input/HZ19_MES_FACS.csv")
HZ19_CEWE_qPCR        <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/ef8575f03ff3a167fda850fd0e93c5de13843aa6/data_input/Eimeria_detection/HZ19_CEWE_qPCR.csv")
HZ19_immuno           <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/ef8575f03ff3a167fda850fd0e93c5de13843aa6/data_input/HZ19_immuno.csv")
HZ19_Dissections      <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/ef8575f03ff3a167fda850fd0e93c5de13843aa6/data_input/Mouse_data/HZ19_Dissections.csv")
HZ10_19_Genotypes     <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/ef8575f03ff3a167fda850fd0e93c5de13843aa6/data_input/Mouse_data/HZ10-19_Genotypes.csv")

Eimeria_Detection_merge <- full_join(Eimeria_Detection, HZ_CEWE_ELISA)
Eimeria_Detection_merge <- full_join(Eimeria_Detection, Eimeria_Summary)
Eimeria_Detection_merge <- full_join(Eimeria_Detection_merge, Eimeria_Summary)
Eimeria_Detection_merge <- full_join(Eimeria_Detection_merge, HZ19_wild_immuno_long)
Eimeria_Detection_merge <- full_join(Eimeria_Detection_merge, HZ19_immuno_long)
Eimeria_Detection_merge <- full_join(Eimeria_Detection_merge, HZ19_MES_FACS)
Eimeria_Detection_merge <- full_join(Eimeria_Detection_merge, HZ19_CEWE_qPCR)
Eimeria_Detection_merge <- full_join(Eimeria_Detection_merge, HZ19_immuno)
Eimeria_Detection_merge <- full_join(Eimeria_Detection_merge, HZ19_Dissections)
Eimeria_Detection_merge <- Eimeria_Detection_merge %>% select(-X)
#Eimeria_Detection_merge %>% duplicated(Eimeria_Detection_merge$Mouse_ID, incomparables = F, fromLast = T )
Eimeria_Detection_merge <- Eimeria_Detection_merge %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 
Eimeria_Detection_merge$OPG <- as.integer(Eimeria_Detection_merge$OPG)
vis_miss(Cryp_meets_Eim, sort = T, cluster = T)
Cryp_meets_Eim <- full_join(Crypto_Detection, Eimeria_Detection_merge) %>% 
  select(-X) %>%
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) %>% 
  mutate(Eimeria_Positive = case_when(eimeriaSpecies == "Negative" ~ F,
                                      eimeriaSpecies == "Double" ~ T,
                                      eimeriaSpecies == "Double_ferrisi_vermiformis" ~ T,
                                      eimeriaSpecies == "Double_tbd" ~ T,
                                      eimeriaSpecies == "E_falciformis" ~ T,
                                      eimeriaSpecies == "E_ferrisi" ~ T,
                                      eimeriaSpecies == "E_vermiformis" ~ T,
                                      eimeriaSpecies == "Other" ~ T,
                                      Ct.Eimeria > 0 ~ T,
                                      OPG >0 ~ T,
                                      MC.Eimeria == T ~ T),
         Cryp_Eim_Positive = case_when(Crypto_Positive  == F ~ F,
                                       Eimeria_Positive == F ~ F,
                                       Crypto_Positive  == T & Eimeria_Positive == T ~ T,
                                       Crypto_Positive  == T & is.na(Eimeria_Positive) ~ F,
                                       is.na(Crypto_Positive) & Eimeria_Positive == T ~ F),
         Infection = case_when(Crypto_Positive  == T & Eimeria_Positive == T ~ "Cryp_Eim_POS",
                               Crypto_Positive  == T ~ "Crypto_Pos",
                               Eimeria_Positive == T ~ "Eimeria_Pos",
                               Crypto_Positive  == F & Eimeria_Positive == F ~ "Cryp_Eim_NEG",
                               Eimeria_Positive == F ~ "Eimeria_Neg",
                               Crypto_Positive  == F ~ "Crypto_Neg"))



## Filter and exclude NAs, extreme cut down due to Marker values missing for most samples (2019 only)
Cryp_meets_Eim <- Cryp_meets_Eim %>%  
  select(Mouse_ID, HI, Year, Oocyst_Predict, Ct_mean, Ct.Eimeria, OPG, MC.Eimeria, eimeriaSpecies, Crypto_Positive, 
         IFNy, CD4, Treg, Div_Treg, Treg17, Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IL17A_CD4, IFNy_CD8,
         Eimeria_Positive, Cryp_Eim_Positive, Infection) %>%
  pivot_longer(cols = c(IFNy, CD4, Treg, Div_Treg, Treg17, Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IL17A_CD4, IFNy_CD8),
               names_to = "Marker",
               values_to = "Value") %>%
  filter(HI >= 0)
         #!is.na(Infection),
         #!is.na(Value))
         #!is.na(MC.Eimeria))

vis_miss(Cryp_meets_Eim, sort = T)






## how is the distribution of Eimeria_pos vs. Crypto_pos? do we have double-infected Samples?
Cryp_meets_Eim %>% 
  count(Eimeria_Positive, Crypto_Positive, Infection)
  

#Cryp_meets_Eim %>% count(Mouse_ID) %>% mutate(count = case_when(n == 14 ~ n/14,
#                                                                n == 13 ~ n/13,
#                                                                n == 1 ~ n))

# Eimeria_Positive Crypto_Positive    Infection         count
# <lgl>            <lgl>             <chr>              <dbl>
# 1 FALSE            FALSE           Cryp_Eim_NEG       346
# 2 FALSE            TRUE            Crypto_Pos          61
# 3 FALSE            NA              Eimeria_Neg        627
# 4 TRUE             FALSE           Eimeria_Pos        157
# 5 TRUE             TRUE            Cryp_Eim_POS        39
# 6 TRUE             NA              Eimeria_Pos         87
# 7 NA               FALSE           Crypto_Neg          18
# 8 NA               TRUE            Crypto_Pos           1
# 9 NA               NA              NA                 171


Cryp_meets_Eim %>% 
  count(Infection) %>%
  mutate(count = n /14) %>% 
  select(-n) %>%
  ggplot(aes(x = Infection, y = count, fill = Infection)) + 
  geom_col() +
  geom_text(aes(label = count), vjust = 1.4, col = "white")


## faceted by Markers
Cryp_meets_Eim %>%
  filter(Ct_mean > 0 & Ct.Eimeria > 0 | Ct_mean > 0 | Ct.Eimeria > 0) %>%
  ggplot(aes(y = Value, col = Infection)) +
  geom_point(aes(x = Ct_mean)) +
  geom_point(aes(x = Ct.Eimeria)) +
  facet_wrap(~Marker, scales = "free_y", nrow = 3)

## faceted by Infection
Cryp_meets_Eim %>%
  #filter(Ct_mean > 0 & Ct.Eimeria > 0 | Ct_mean > 0 | Ct.Eimeria >0) %>%
  filter(Infection == c("Cryp_Eim_POS", "Eimeria_Pos", "Crypto_Neg")) %>% 
  ggplot(aes(y = Value, col = Marker)) +
  geom_point(aes(x = Ct_mean)) +
  geom_point(aes(x = Ct.Eimeria)) +
  facet_wrap(~Infection, scales = "free_y", nrow = 3)









# Principal Component Analysis ################################################

pca.plot.data <- Cryp_meets_Eim %>%
  pivot_wider(values_from = "Value", names_from = "Marker") %>%
  filter(HI >= 0 & Crypto_Positive == T & Eimeria_Positive == T, CD4 >= 0)


pca <- Cryp_meets_Eim %>%
  pivot_wider(values_from = "Value", names_from = "Marker") %>%
  filter(HI >= 0 & Crypto_Positive == T & Eimeria_Positive == T, CD4 >= 0) %>%
  select(c(CD4, Treg, Div_Treg, Treg17, Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IL17A_CD4, IFNy_CD8)) %>%
  scale() %>%
  prcomp()

pca.plot.data <- pca$x %>%
  as_tibble() %>%
  bind_cols(pca.plot.data, .)

#### normal plot
pca.plot.data %>%
  ggplot(aes(PC1, PC2)) +
  #scale_x_reverse() +
  geom_point(size = 1)

pca.plot.data %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(size = 1)





#### color the dots in your PCA using the values of different signalling markers
pca_col_data <- pca.plot.data %>%
  pivot_longer(values_to = "Value",
               names_to = "Marker",
               cols = c(CD4, Treg, Div_Treg, Treg17, Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IL17A_CD4, IFNy_CD8)) %>%
  dplyr::select(PC1, PC2, Value, Marker)

pca_col_data %>%
  ggplot(aes(PC1, PC2, col = Value)) +
  geom_point(size = 1) +
  scale_color_viridis(option = "inferno") +
  facet_wrap(~ Marker)


# Clustering
#### k Means
#### cluster into 3 clusters color code the pca plot by cluster
#### cluster cells multiple times, what happens?

fitk <- pca.plot.data %>%
  select(c(IFNy, CD4, Treg, Div_Treg, Treg17, Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IL17A_CD4, IFNy_CD8)) %>%
  scale() %>%
  kmeans(centers = 3)

pca.plot.data <- pca.plot.data %>%
  mutate(cluster = as.factor(fitk$cluster))

#display.brewer.all()
pca.plot.data %>%
  ggplot(aes(PC1, PC2, col = cluster)) +
  geom_point(size = 1) +
  scale_color_brewer(palette = "Dark2",
                     direction = 1)


#Heatmap representation of perturbation data
#### visualize changes in signalling by cluster, treatment, cell line w/ HEATMAP
#### group by cluster, treatment, cell line
#### calculate mean signal values within groups for markers of choice
#### plot w/ geom_tile & facet_grid

heatmap_data <- pca.plot.data %>%
  pivot_longer(c("CD4", "Treg", "Div_Treg", "Treg17", "Th1", "Th17", "Div_Th17", "CD8", "Act_CD8", "Div_Act_CD8", "IFNy_CD4", "IL17A_CD4", "IFNy_CD8"),
               names_to = "Marker",
               values_to = "Value") %>%
  group_by(cluster, Infection, Marker) %>%
  summarise(Value = mean(Value))

heatmap_data %>%
  ggplot(aes(Marker,Infection, fill = Value)) +
  geom_tile() +
  facet_grid(cluster ~ Infection) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis(option = "inferno")


# PCA II not limited to Crypto-Eim-Positive Infections ################################################

pca.plot.data <- Cryp_meets_Eim %>%
  pivot_wider(values_from = "Value", names_from = "Marker") %>%
  filter(HI >= 0 & CD4 >= 0)


pca <- Cryp_meets_Eim %>%
  pivot_wider(values_from = "Value", names_from = "Marker") %>%
  filter(HI >= 0 & CD4 >= 0) %>%
  select(c(CD4, Treg, Div_Treg, Treg17, Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IL17A_CD4, IFNy_CD8)) %>%
  scale() %>%
  prcomp()

pca.plot.data <- pca$x %>%
  as_tibble() %>%
  bind_cols(pca.plot.data, .)

#### normal plot
pca.plot.data %>%
  ggplot(aes(PC1, PC2)) +
  #scale_x_reverse() +
  geom_point(size = 1)

pca.plot.data %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(size = 1)


pca.var <- pca$sdev^2
pca.var.per <- round(pca.var / sum(pca.var) * 100, 1)
barplot(pca.var.per, main = "Scree Plot", xlab = "Principal component",
        ylab = "Percent Variation")

loading_scores <- pca$rotation[,1]
scores <- abs(loading_scores)
score_ranked <- sort(scores, decreasing = T)
top10 <- names(score_ranked[1:13])
pca$rotation[top10,1]
pca.var.per.df <- cbind(pca.var.per, top10)
pca.var.per.df <- data.frame(pca.var.per.df)
pca.var.per.df$pca.var.per <- as.numeric(as.character(pca.var.per.df$pca.var.per))
pca.var.per.df$top10 <- factor(pca.var.per.df$top10,levels = top10)
ggplot(pca.var.per.df, aes(top10, pca.var.per, fill = top10)) +
  geom_col() +
  theme(strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold")) +
  scale_x_discrete(name ="Marker") +
  scale_y_discrete(name ="% of variation explained") +
  labs(fill = "") +
  
  geom_text(aes(label = pca.var.per)) +
  theme_bw()


#### color the dots in your PCA using the values of different signalling markers
pca_col_data <- pca.plot.data %>%
  pivot_longer(values_to = "Value",
               names_to = "Marker",
               cols = c(CD4, Treg, Div_Treg, Treg17, Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IL17A_CD4, IFNy_CD8)) %>%
  dplyr::select(PC1, PC2, Value, Marker)

pca_col_data %>%
  ggplot(aes(PC1, PC2, col = Value)) +
  geom_point(size = 1) +
  scale_color_viridis(option = "inferno") +
  facet_wrap(~ Marker)


# Clustering
#### k Means
#### cluster into 3 clusters color code the pca plot by cluster
#### cluster cells multiple times, what happens?

fitk <- pca.plot.data %>%
  select(c(CD4, Treg, Div_Treg, Treg17, Th1, Th17, Div_Th17, CD8, Act_CD8, Div_Act_CD8, IFNy_CD4, IL17A_CD4, IFNy_CD8)) %>%
  scale() %>%
  kmeans(centers = 3)

pca.plot.data <- pca.plot.data %>%
  mutate(cluster = as.factor(fitk$cluster))

#display.brewer.all()
pca.plot.data %>%
  ggplot(aes(PC1, PC2, col = cluster)) +
  geom_point(size = 1) +
  scale_color_brewer(palette = "Dark2",
                     direction = 1)


#Heatmap representation of perturbation data
#### visualize changes in signalling by cluster, treatment, cell line w/ HEATMAP
#### group by cluster, treatment, cell line
#### calculate mean signal values within groups for markers of choice
#### plot w/ geom_tile & facet_grid

heatmap_data <- pca.plot.data %>%
  pivot_longer(c("CD4", "Treg", "Div_Treg", "Treg17", "Th1", "Th17", "Div_Th17", "CD8", "Act_CD8", "Div_Act_CD8", "IFNy_CD4", "IL17A_CD4", "IFNy_CD8"),
               names_to = "Marker",
               values_to = "Value") %>%
  group_by(cluster, Infection, OPG, Marker) %>%
  summarise(Value = mean(Value))

heatmap_data %>%
  ggplot(aes(Marker, Infection, fill = Value)) +
  geom_tile() +
  facet_grid(cluster ~ OPG) +
  #facet_grid(~ cluster) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis(option = "inferno")

