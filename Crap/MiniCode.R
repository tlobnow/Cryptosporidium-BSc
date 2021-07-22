library(ggplot2)
library(dplyr)
library(ggeffects)
#library(DHARMa) # investigation of models
#library(AICcmodavg) # which model has the best AIC score (complexity vs. accuracy)
library(visdat)
library(RColorBrewer) # display.brewer.all()
library(mosaic)
library(pheatmap)
library(tidyverse)
library(assertive)
library(stringr)
library(stringdist)
library(fuzzyjoin)
library(reclin)
library(viridisLite)
library(viridis)
library(cowplot)
library(remotes)
library(uwot)

#ABI_thSC_log2 <- read.csv("ABI_thSC_log2.csv")
#all_Samples   <- read.csv("all_Samples_Complete.csv")

# remove duplicates:
#all_Samples_unique <- all_Samples %>%
#  distinct(Mouse_ID, .keep_all = TRUE)
#all_Samples_unique %>%
#  count(Mouse_ID) %>%
#  filter(n > 1)
#write.csv(all_Samples_unique, "all_Samples_Complete.csv")

getwd()
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
summary(linear_model0)
ggpredict(linear_model0)



#nf_full1 %>%
#  ggplot(aes(Mouse_ID, Oocyst_Predict, col = Year)) +
#  geom_bar(stat = "identity")

#models <- list(linear_model0, linear_model1)
#model.names <- c("linMod0", "linMod1")
#anova(linear_model0, linear_model1)
#aictab(cand.set = models, modnames = model.names, sort = T)
#Model selection based on AICc:
#           K      AICc    Delta_AICc    AICcWt   Cum.Wt     LL
#linMod1    3     12.27         0.00      1      1          2.86
#linMod0    3     47.50         35.23     0      1        -20.12

# predict function
f_ABI_Best_thSC$predicted   <- 2^predict(linear_model0)
# predict model on samples (all Samples)
Oocyst_Predict   <-  2^predict(linear_model0, newdata = nf_full1)
nf_full1 <- data.frame(nf_full1, Oocyst_Predict)
nf_full1 <- nf_full1 %>%
  mutate(Oocyst_Predict = replace(Oocyst_Predict, Oocyst_Predict == '4292821751815.77', '0'))
write.csv(nf_full1, "nf_full1_Corrected.csv")

#full1 <- nf_full1 %>%
#  filter(Ct_mean > 0)
#write.csv(full1, "full1.csv")



#fd_ABI_log2   <- data.frame(fd_ABI_log2, Oocyst_Prediction) %>%
#  mutate(ID = row_number()) %>%
#  mutate(Delta_Ct = sqrt((Ct_1 - Ct_2)^2)) %>%
#  select(ID, Mouse_ID, Year, Ct_1, Ct_2, Ct_3, Ct_mean, Delta_Ct, 
#         Oocyst_Prediction, Transect, Latitude, Longitude, Sex, HI) %>%
#  mutate(goofy = Delta_Ct > 3)

# csv write.csv(fd_ABI_log2, "fd_ABI_log2.csv")

  
#new_Oocyst_Predict   <-  2^predict(linear_model0, newdata = new_Plate6)
#new_Plate6 <- data.frame(new_Plate6, new_Oocyst_Predict)
#new_Plate6 <- new_Plate6 %>%
#  mutate(new_Oocyst_Predict = replace(new_Oocyst_Predict, new_Oocyst_Predict == '4292821751815.77', '0'))
new_Plate6 <- read.csv("new_Plate6.csv")
#write.csv(new_Plate6, "new_Plate6.csv")


#adj, includes negative Samples
new_Plate6 %>%
  #filter(Ct_mean != 0) %>%
  ggplot(aes(x = Mouse_ID, fill = Consistent)) +
  geom_bar(aes(y = Ct_Diff_Adj), stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Difference of Ct-Values in double-tested Samples") +
  ylim(c(-28,26)) +
  labs(y = "Ct Difference")
 
# adj0, excludes negative Samples
new_Plate6 %>%
  ggplot(aes(x = Mouse_ID, fill = Consistent)) +
  geom_bar(aes(y = Ct_Diff_Adj0), stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Adjusted old Ct_values vs. new Ct_values, negative Samples excluded") +
  ylim(c(-28,26))

# vis the HI Level
nf_full1 %>%
  filter(Ct_mean > 0) %>%
  ggplot(aes(Mouse_ID, Ct_mean, col = HI_Level)) +
  geom_point() +
  #scale_color_viridis(option = "inferno") +
  facet_grid(~ HI_Level, nrow(1)) +
  #facet_wrap(~Year, nrow = 1) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())


full1 <- nf_full1 %>%
  filter(Ct_mean > 0)

min(full1$Ct_mean)
max(full1$Ct_mean)

# visualize HI Level
a <- nf_full1 %>%
  filter(Ct_mean > 0) %>%
  ggplot(aes(Mouse_ID, Ct_mean, col = Ct_mean)) +
  geom_point() +
  scale_color_viridis(option = "inferno") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("Ct-Values of Crypto-Positive Samples")
save_plot(a, filename = "Ct_vs_MouseID.jpg", base_height = 6, base_width = 9)


################################################################################
# Visualize HI Data ############################################################
      # HI = 1 = Mmm
      # Hi = 0 = Mmd

# HI vs. Ct - ggsave("HybridIndex_HI_Ct_all_positive_Samples.jpg", width = 8, height = 5)
full1 %>%
  ggplot(aes(HI, Ct_mean, col = HI)) +
  geom_point(position = position_jitter(0), size = 2, alpha = 0.5) +
  scale_color_gradient(low = "blue", high = " red") +
  ggtitle("Crypto Samples HI vs. Ct")

# HI vs.Ct - ggsave("HybridIndex_HI_Ct_all_Samples.jpg", width = 8, height = 5)
nf_full1 %>%
  ggplot(aes(HI, Ct_mean, col = HI)) +
  geom_point(position = position_jitter(0), size = 2, alpha = 0.5) +
  scale_color_gradient(low = "blue", high = " red") +
  ggtitle("All Samples HI vs. Ct")


nf_full1 %>%
  ggplot(aes(HI, Ct_mean, col = HI)) +
  geom_tile() +
  geom_point() +
  facet_wrap(~Year) +
  scale_color_gradient(low = "blue", high = " red")


nf_full1 <- nf_full1 %>%
  filter(!is.na(Year)) %>%
  mutate(Oocyst_Predict = replace(Oocyst_Predict, Oocyst_Predict == '4292821751815.77', '0')) %>%
  arrange(Mouse_ID)

plot(full1$Year, full1$Ct_mean, 
     main = "Ct_means per Year",
     xgap.axis = 1,
     xlab = "Year",
     ylab = "Ct mean")

plot(nf_full1$Year, nf_full1$Ct_mean, 
     main = "Ct_means per Year",
     xlab = "Year",
     ylab = "Ct mean")
plot(nf_full1$Year, nf_full1$HI,
     main = "HI per Year",
     xlab = "Year",
     ylab = "HI")

full1 %>%
ggplot(aes(Year, HI, col = HI)) +
  geom_point()

nf_full1

nf_full1 %>%
  ggplot(aes(Mouse_ID, HI, col = HI)) +
  geom_point() +
  facet_grid(~Year) +
  scale_color_gradient(low = "blue", high = " red") +
  ggtitle("all Samples vs. Hybrid Index")+
  theme(axis.text.x = element_blank())

nf_full1 %>%
  filter(Oocyst_Predict != 0) %>%
  ggplot(aes(HI, Ct_mean, col = HI)) +
  geom_point() +
  #facet_grid(~Year) +
  scale_color_gradient(low = "blue", high = " red") +
  ggtitle("Hybrid index vs. Ct across the years, positive Samples")

HI_pull <- nf_full1 %>%
  filter(!is.na(HI)) %>%
  pull(HI)

HI_R_pull <- nf_full1 %>%
  filter(!is.na(HI_Ratio)) %>%
  pull(HI_Ratio)

Oocyst_pull <- nf_full1 %>%
  filter(!is.na(Oocyst_Predict)) %>%
  pull(Oocyst_Predict)

Ct_pull <- nf_full1 %>%
  filter(!is.na(Ct_mean)) %>%
  pull(Ct_mean)

t.test(HI_pull, HI_R_pull)
# p-value < 2.2e-16

t.test(HI_pull, Oocyst_pull)
# p-value = 0.02231

t.test(Ct_pull, Oocyst_pull)
# p-value = 0.02231

t.test(HI_pull, Ct_pull)
# p-value < 2.2e-16

mus_pull <- nf_full1 %>%
  filter(!is.na(ratio_mus_caught)) %>%
  pull(ratio_mus_caught)
