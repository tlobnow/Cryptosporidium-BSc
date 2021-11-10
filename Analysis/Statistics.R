## STATISTICS MANN-WHITHNEY-U TEST



## GP60
HMHZ <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/HMHZ_Samples_Locations.csv", na.strings=c(""," ","NA")) %>% filter(!is.na(Longitude))
HMHZ <- HMHZ %>% 
  mutate(GP60_Ssp = case_when(Mouse_ID %in% c("CR_2085",
                                              "CR_2084",
                                              "CR_2090",
                                              "CR_2125",
                                              "CR_2126",
                                              "CR_2127",
                                              "CR_2128",
                                              "CR_2149",
                                              "CR_2208",
                                              "CR_2206") ~ "IXa",
                              Mouse_ID %in% c("AA_0523", "AA_0282") ~ "IXb.1",
                              Mouse_ID %in% c("G_2136",
                                              "CR_2163",
                                              "G_2099",
                                              "CR_2147",
                                              "G_2135",
                                              "G_2177",
                                              "G_2179",
                                              "G_2174",
                                              "AA_0534",
                                              "G_2103",
                                              "AA_0537",
                                              "AA_0585",
                                              "CR_4293",
                                              "G_2116",
                                              "G_2120",
                                              "G_2134",
                                              "G_2160",
                                              "G_2108",
                                              "G_2169",
                                              "G_3224") ~ "IXb.2",
                              Mouse_ID %in% c("AA_0209",
                                              "AA_0144",
                                              "AA_0325",
                                              "AA_0545",
                                              "AA_0553",
                                              "AA_0554",
                                              "AA_0580",
                                              "AA_0667",
                                              "AA_0689",
                                              "AA_0793",
                                              "AA_0805",
                                              "AA_0900",
                                              "NZ_1633",
                                              "NZ_1634",
                                              "NZ_1635",
                                              "NZ_1636",
                                              "NZ_1639",
                                              "NZ_1640",
                                              "NZ_1641",
                                              "NZ_1642",
                                              "NZ_1644",
                                              "Tyz-GP60") ~ "IXb.3"))


## for AA_0900 we will assume an HI of AA_0144, as they were sampled geographically very close to each other (HI =  0.869565217)

HMHZ$HI[ HMHZ$Mouse_ID %in% "AA_0900"] <- 0.869565217

## we will filter out 
## all Samples that didn't have their HI evaluated (e.g. New Zealand samples, we could also assume an HI of 0, but it simply wasn't measured)
## all Samples w/o geographic data (lat,lon missing)

## we will select necessary columns:
## Mouse_ID
## HI
## GP60_ssp

HMHZ <- HMHZ %>% 
  filter(!is.na(HI), !is.na(Latitude)) %>%
  select(Mouse_ID, HI, GP60_Ssp)

## Analysis for remaining 61/75 samples

IXa                 <- HMHZ %>% filter(GP60_Ssp == "IXa")# 10 samples
IXa_Subtype_family  <- IXa[IXa$GP60_Ssp == "IXa", "HI"]
IXb                 <- HMHZ %>% filter(GP60_Ssp %in% c("IXb.1", "IXb.2", "IXb.3")) # 33 samples
IXb_Subtype_family  <- IXb[IXb$GP60_Ssp %in% c("IXb.1", "IXb.2", "IXb.3"), "HI"]

t.test(IXa_Subtype_family, IXb_Subtype_family)
#data:  IXa_Subtype_family and IXb_Subtype_family
#t = 9.5924, df = 32.015, p-value = 6.183e-11
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.4927246 0.7583926
#sample estimates:
#  mean of x mean of y 
#0.9710000 0.3454414 

wilcox.test(IXa_Subtype_family, IXb_Subtype_family)
#data:  IXa_Subtype_family and IXb_Subtype_family
#W = 330, p-value = 1.868e-06
#alternative hypothesis: true location shift is not equal to 0




IXb.1 <- HMHZ %>% filter(GP60_Ssp %in% c("IXb.1")) # 2 samples
IXb.1_Subtype_family <- IXb.1[IXb.1$GP60_Ssp == "IXb.1", "HI"]
IXb.2 <- HMHZ %>% filter(GP60_Ssp %in% c("IXb.2")) # 19 samples
IXb.2_Subtype_family <- IXb.2[IXb.2$GP60_Ssp == "IXb.2", "HI"]


t.test(IXb.1_Subtype_family, IXb.2_Subtype_family)
#data:  IXb.1_Subtype_family and IXb.2_Subtype_family
#t = 1.9548, df = 1.0106, p-value = 0.299
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.308596  1.798596
#sample estimates:
#  mean of x mean of y 
#0.295     0.050 

wilcox.test(IXb.1_Subtype_family, IXb.2_Subtype_family)

#data:  IXb.1_Subtype_family and IXb.2_Subtype_family
#W = 38, p-value = 0.02236
#alternative hypothesis: true location shift is not equal to 0




IXb.3 <- HMHZ %>% filter(GP60_Ssp %in% c("IXb.3")) # 12 samples
IXb.3_Subtype_family <- IXb.3[IXb.3$GP60_Ssp == "IXb.3", "HI"]
IXb.12 <- HMHZ %>% filter(GP60_Ssp %in% c("IXb.1","IXb.2")) # 21 samples <- since IXb.1 is super tiny with 2 samples, we will join these two for analysis!
IXb.12_Subtype_family <- IXb.12[IXb.12$GP60_Ssp %in% c("IXb.1","IXb.2"), "HI"]


t.test(IXb.12_Subtype_family, IXb.3_Subtype_family)
#data:  IXb.12_Subtype_family and IXb.3_Subtype_family
#t = -27.359, df = 29.437, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.8042000 -0.6923942
#sample estimates:
#  mean of x  mean of y 
#0.07333333 0.82163043 

wilcox.test(IXb.12_Subtype_family, IXb.3_Subtype_family)
#data:  IXb.12_Subtype_family and IXb.3_Subtype_family
#W = 0, p-value = 2.214e-06
#alternative hypothesis: true location shift is not equal to 0




################################################################################
################################################################################

#### COWP 

HMHZ <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/HMHZ_Samples_Locations.csv", na.strings=c(""," ","NA")) %>% filter(!is.na(Longitude))
HMHZ <- HMHZ %>% 
  mutate(COWP_Ssp = case_when(Mouse_ID %in% c("CR_2085",
                                              "CR_2084",
                                              "AA_0578",
                                              "AA_0666",
                                              "AA_0580",
                                              "AA_0667",
                                              "AA_0679",
                                              "AA_0601",
                                              "AA_0553",
                                              "AA_0545",
                                              "AA_0546",
                                              "AA_0555",
                                              "AA_0557",
                                              "AA_0559",
                                              "CR_2090",
                                              "CR_2125",
                                              "CR_2126",
                                              "CR_2128",
                                              "CR_2149",
                                              "CR_2208",
                                              "CR_2206",
                                              "AA_0589") ~ "C1",
                              
                              Mouse_ID %in% c("NZ_1633",
                                              "NZ_1634",
                                              "NZ_1635",
                                              "NZ_1637",
                                              "NZ_1638",
                                              "NZ_1640",
                                              "NZ_1641",
                                              "NZ_1642",
                                              "AA_0534",
                                              "G_3224",
                                              "G_2177",
                                              "G_2179",
                                              "G_2174",
                                              "G_2120",
                                              "G_2134",
                                              "G_2135",
                                              "G_2136",
                                              "G_2160",
                                              "G_2108",
                                              "G_2169",
                                              "G_2110",
                                              "G_2103",
                                              "G_2099",
                                              "CR_4293",
                                              "CR_2163",
                                              "AA_0537",
                                              "AA_0571",
                                              "AA_0585",
                                              "AA_0660",
                                              "AA_0554") ~ "C2"))


## we will filter out 
## all Samples that didn't have their HI evaluated (e.g. New Zealand samples, we could also assume an HI of 0, but it simply wasn't measured)
## all Samples w/o geographic data (lat,lon missing)

## we will select necessary columns:
## Mouse_ID
## HI
## COWP_Ssp

HMHZ <- HMHZ %>% 
  filter(!is.na(HI), !is.na(Latitude)) %>%
  select(Mouse_ID, HI, COWP_Ssp)

## Analysis for remaining 61/75 samples

C1                 <- HMHZ %>% filter(COWP_Ssp == "C1") # 22 samples
C1_Clade  <- C1[C1$COWP_Ssp == "C1", "HI"]

C2                 <- HMHZ %>% filter(COWP_Ssp == "C2") # 22 samples
C2_Clade  <- C2[C2$COWP_Ssp == "C2", "HI"]



    t.test(C1_Clade, C2_Clade)
    #data:  C1_Clade and C2_Clade
    #t = 10.545, df = 41.142, p-value = 2.896e-13
    #alternative hypothesis: true difference in means is not equal to 0
    #95 percent confidence interval:
    #  0.5784457 0.8524633
    #sample estimates:
    # mean of x mean of y 
    #0.8413636 0.1259091 
    
    wilcox.test(C1_Clade, C2_Clade)
    #data:  C1_Clade and C2_Clade
    #W = 465.5, p-value = 1.373e-07
    #alternative hypothesis: true location shift is not equal to 0


