library(dplyr)
library(tidyverse)
library(visdat)
library(ggplot2)
library(leaflet)

Crypto_Detection <- read.csv("Crypto_qPCR_Results.csv") %>% select(-Oocyst_Predict)
Crypto_Detection  <- Crypto_Detection[!names(Crypto_Detection) == "X"]
Crypto_Detection  <- Crypto_Detection[!names(Crypto_Detection) == "X.1"]

Crypto_Detection %>%
  count(Machine)

## Predict Oocysts
    ABI_Best_thSC     <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Crap/ABI_Best_SC.csv")
    f_ABI_Best_thSC   <-  filter(ABI_Best_thSC, Ct_mean > 0)
    linear_model0     <- lm(log2(Amount_Oocysts) ~ Ct_mean, data = f_ABI_Best_thSC)
    Oocyst_Predict    <- 2^predict(linear_model0, newdata = Crypto_Detection)
    
    Crypto_Detection <- data.frame(Crypto_Detection, Oocyst_Predict)
    Crypto_Detection <- Crypto_Detection %>%
      mutate(Oocyst_Predict = replace(Oocyst_Predict, Oocyst_Predict == "4292821751815.77", "0"))
    Crypto_Detection$Oocyst_Predict <- as.numeric(Crypto_Detection$Oocyst_Predict)

## write csv
    write.csv(Crypto_Detection, "Crypto_Detection.csv")
    
    
## Visualize
    

    
    
    
    
    
        
