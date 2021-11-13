## STATISTICS MANN-WHITHNEY-U TEST



## GP60
Clades <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Clade_Memberships.csv")


IXa                 <- Clades %>% filter(GP60 == "IXa")# 10 samples
IXa_Subtype_family  <- IXa[IXa$GP60 == "IXa", "HI"]
IXb                 <- Clades %>% filter(GP60 %in% c("IXb.1", "IXb.2", "IXb.3")) # 43 samples
IXb_Subtype_family  <- IXb[IXb$GP60 %in% c("IXb.1", "IXb.2", "IXb.3"), "HI"]

t.test(IXa_Subtype_family, IXb_Subtype_family)
# p-value = 3.317e-16

wilcox.test(IXa_Subtype_family, IXb_Subtype_family)
# p-value = 8.926e-07




IXb.2 <- Clades %>% filter(GP60 %in% c("IXb.2")) # 2 samples
IXb.2_Subtype_family <- IXb.2[IXb.2$GP60 == "IXb.2", "HI"]
IXb.3 <- Clades %>% filter(GP60 %in% c("IXb.3")) # 19 samples
IXb.3_Subtype_family <- IXb.3[IXb.3$GP60 == "IXb.3", "HI"]


t.test(IXb.2_Subtype_family, IXb.3_Subtype_family)
# 0.299

wilcox.test(IXb.2_Subtype_family, IXb.3_Subtype_family)
# p-value = 0.02236





IXb.1 <- Clades %>% filter(GP60 %in% c("IXb.3")) # 12 samples
IXb.1_Subtype_family <- IXb.3[IXb.3$GP60 == "IXb.3", "HI"]
IXb.23 <- Clades %>% filter(GP60 %in% c("IXb.3","IXb.2")) # 21 samples <- since IXb.1 is super tiny with 2 samples, we will join these two for analysis!
IXb.23_Subtype_family <- IXb.12[IXb.12$GP60 %in% c("IXb.1","IXb.2"), "HI"]


t.test(IXb.23_Subtype_family, IXb.1_Subtype_family)
# p-value = 0.8474

wilcox.test(IXb.23_Subtype_family, IXb.1_Subtype_family)
#p-value = 0.3365




################################################################################
################################################################################

#### COWP 
Clades <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Clade_Memberships.csv")

C1                 <- Clades %>% filter(COWP == "C1") # 22 samples
C1_Clade  <- C1[C1$COWP == "C1", "HI"]

C2                 <- Clades %>% filter(COWP == "C2") # 22 samples
C2_Clade  <- C2[C2$COWP == "C2", "HI"]



    t.test(C1_Clade, C2_Clade)
    # p-value < 2.2e-16
    
    wilcox.test(C1_Clade, C2_Clade)
    #p-value = 8.326e-10


HMHZ <- HMHZ %>% select(-GP60_Subtype_Protein)
    write.csv(HMHZ, "Analysis/HMHZ_Samples_Locations.csv")
    


################################################################################

a <- SOTA_Data_Product[SOTA_Data_Product$Mouse_ID %like% "AA_", ]

a <- a %>% filter(Ct_mean > 0)
b <- Crypto_Detection %>% filter(Ct_mean > 0)
c <- Crypto_Detection %>% filter(Ct_mean < 29, Ct_mean > 0)

  ggplot(b, aes(x = Ct_mean, y = Oocyst_Predict)) +
  geom_point() +
  geom_density_2d() +
  geom_label(data = c, aes(label = Mouse_ID), nudge_x = +1, nudge_y = +100000) +
    ggtitle("Distribution of Cryptosporidium Infection in the HMHZ", 
            subtitle = "Highest infections labeled with Mouse_ID")



################################################################################
  
Clades <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Clade_Memberships.csv")

MS1 <- Clades %>% filter(MSC6 == "MS1")
MS1 <- MS1[MS1$MSC6 == "MS1", "HI"]

MS2 <- Clades %>% filter(MSC6 == "MS2")
MS2 <- MS2[MS2$MSC6 == "MS2", "HI"]



# NOT SIGNIFICANT
t.test(MS1, MS2) #p-value = 0.09302
wilcox.test(MS1, MS2) # p-value = 0.146



################################################################################

Clades <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Clade_Memberships.csv")

CP1 <- Clades %>% filter(CP56 == "CP1")
CP1 <- CP1[CP1$CP56 == "CP1", "HI"]

CP2 <- Clades %>% filter(CP56 == "CP2")
CP2 <- CP2[CP2$CP56 == "CP2", "HI"]



t.test(CP1, CP2) # p-value = 0.02529
wilcox.test(CP1, CP2) #p-value = 0.3333



################################################################################

Clades <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Clade_Memberships.csv")

G1 <- Clades %>% filter(GST == "G1")
G1 <- G1[G1$GST == "G1", "HI"]

G2 <- Clades %>% filter(GST == "G2")
G2 <- G2[G2$GST == "G2", "HI"]


t.test(G1, G2) #p-value = 0.02503
wilcox.test(G1, G2) #p-value = 0.07652


################################################################################

Clades <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Clade_Memberships.csv")

S1 <- Clades %>% filter(SKSR == "S1")
S1 <- S1[S1$SKSR == "S1", "HI"]

S2 <- Clades %>% filter(SKSR == "S2")
S2 <- S2[S2$SKSR == "S2", "HI"]


t.test(S1, S2) # not enough observations
wilcox.test(S1, S2)  # not enough observations
################################################################################

Clades <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Clade_Memberships.csv")

ME1 <- Clades %>% filter(MEDLE == "ME1")
ME1 <- ME1[ME1$MEDLE == "ME1", "HI"]

ME2 <- Clades %>% filter(MEDLE == "ME2")
ME2 <- ME2[ME2$MEDLE == "ME2", "HI"]


t.test(ME1, ME2) # not enough observations
wilcox.test(ME1, ME2)  # not enough observations


################################################################################

Clades <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Clade_Memberships.csv")

A1 <- Clades %>% filter(Actin == "A1")
A1 <- A1[A1$Actin == "A1", "HI"]

A2 <- Clades %>% filter(Actin == "A2")
A2 <- A2[A2$Actin == "A2", "HI"]

A3 <- Clades %>% filter(Actin == "A3")
A3 <- A3[A3$Actin == "A3", "HI"]

A12 <- Clades %>% filter(Actin %in% c("A1","A2")) # 21 samples <- since IXb.1 is super tiny with 2 samples, we will join these two for analysis!
A12 <- A12[A12$Actin %in% c("A1","A2"), "HI"]



t.test(A1, A2) # p-value < 2.2e-16
wilcox.test(A1, A2) # p-value = 2.836e-05

t.test(A12, A3) # p-value = 0.0003592
wilcox.test(A12, A3) # p-value = 0.02538

################################################################################

Clades <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Clade_Memberships.csv")

T1 <- Clades %>% filter(TRAP_C1 == "T1")
T1 <- T1[T1$TRAP_C1 == "T1", "HI"]

T2 <- Clades %>% filter(TRAP_C1 == "T2")
T2 <- T2[T2$TRAP_C1 == "T2", "HI"]



t.test(T1, T2) #p-value = 1.822e-05
wilcox.test(T1, T2) # p-value = 5.055e-05

