# Cryptosporidium-BSc

Column Documentation


## BASICS

Mouse_ID == Sample Number 

    AA_ Samples collected since 2016
    SK_ Samples collected prior to 2016
    ZZ_ other rodents          
           
Code == unique Abbreviation of Address

Transect	== Area where Samples were collected (rather roughly, doesn't match exactly)

    HZ_BR   Brandenburg
    HZ_BAV  Bavaria
    HZ_MV   Mecklenburg-Vorpommern
    HZ_CZ   Czech
    ALLO_PL   Poland
            
Sex ==  Female or Male

Longitude ==  11.214 - 21.750

Latitude  ==  47.788 - 57.150

Year  == 2010 - 2019 (SK_ and AA_)

HI_NLoci  == Hybrid Index is defined by 14 Loci that are used to differentiate Mmm and Mmd 
                		0 ... corresponds to a "pure Mmd",
                		14... corresponds to a "pure Mmm"
                		1:13 corresponds to Hybrids of different extend
            
HI  == depends on the Loci of HI_NLoci (although the definition and measurements differ across the
                		years, the HI is the most up-to-date variable to use for visualization, as Jarda tends 
                		to update and adjust measurements over the years)
                		takes values 0-1
            
State == PL, DE, or CZ


## qPCR DATA

qPCR_Date    == as far as known, date of qPCR performed, in case of multiple measurements it refers to the last qPCR performed

Tested_by   == Person that performed DNA extraction and qPCR

        YVT ... Yasmin, Victor, or Tabea performed DNA extraction and qPCR (Eppendorf only)
        Tes ... Tessa performed DNA extraction and qpCR
										
the column can show a variety of combinations, corresponding to multiple measurements:

        Tes
        Tes-Tes
        Tes-Tes-Tes
        YVT
        YVT-YVT
        YVT-YVT-YVT
        YVT-YVT-Tes
        YVT-Tes				

Machine	 == takes values Eppendorf, 7300_RealTimePCR, or Quantstudio1
                	
   BUT if a sample was run multiple times, the description is added
   
        Eppendorf
        Eppendorf + Eppendorf
        Eppendorf + Eppendorf + Eppendorf
        Eppendorf + Eppendorf + Quantstudio1
        Eppendorf + Quantstudio1
        Quantstudio1
        Quantstudio1 + Quantstudio1
        7300_RealTime_PCR + Quantstudio1
        7300_RealTime_PCR + Quantstudio1 + Quantstudio1
  
  
Measurements	 == number of qPCR measurements, 
				1 measurement in old samples stands for duplicates (Eppendorf, YVT)
				1 measurement in new Samples stands for triplicates (ABI machines, Tes)


CT_MEASUREMENTS:

    Ct_1_Ep + Ct_2_Ep = 1st Measurement
    Ct_3_Ep + Ct_4_Ep = 2nd Measurement    
    Ct_5_Ep + Ct_6_Ep = 3rd Measurement    
    Ct_mean_Ep = Mean of all available measurements
      
## FLAGS    
  qPCR quantification curves were double-checked, Flag == curve does not resemble normal quantification curve
  since new measurements since 2019 were performed as triplicates, and I have not performed the old qPCRs,
  the flagging system is only applied to older Samples (Eppendorf qPCRs performed by YVT)
              
		Flag_Ct_1_Ep      
		Flag_Ct_2_Ep  
		Flag_Ct_3_Ep  
		Flag_Ct_4_Ep  
		Flag_Ct_5_Ep  
		Flag_Ct_6_Ep  

Flags == sum of all Flags
              
Flag_perc == Percentage of Flags/ (Measurements * 2)
		This is calculated as Measurements * 2, since 1 measurement equals 2 individual Ct values, therefore also 2 potential flags
		By dividing by Measurements * 2, every Sample that was measured with the Eppendorf machine in previous years hereby receives a 'mistrust' value (percentage). 
		The higher the percentage, the less believable were the curves, thus increasing a need for re-running a qPCR
                 


#### since 2019, qPCRs were performed as Triplicates, not Duplicates
		Ct_1_ABI       
		Ct_2_ABI     
		Ct_3_ABI       
		Ct_mean_1_ABI  == mean of the 1st Measurement
		
		Ct_4_ABI       
		Ct_5_ABI       
		Ct_6_ABI       
		Ct_mean_2_ABI  == mean of the 2nd Measurement
		
		Ct_7_ABI       
		Ct_8_ABI
		Ct_9_ABI
		
		Ct_mean_3_ABI  == mean of the 3rd Measurement

		Ct_mean_ABI    == mean of all ABI measurements

		Ct_mean        == mean of ALL Ct measurements

## Oocyst Prediction
Oocyst Predictions were done based on Standard Curves, that were run in qPCRs with known Oocyst amounts 
(SC was done as 1:8 dilution steps, since this allows for a duplication of DNA amount every 3 quantification cycles == 2^3)
The Ct values should therefore differ in 3-ish steps, the best SC or the mean of all SCs was applied in a linear model:

		linear_model <- lm(log2(Amount_Oocysts) ~ Ct_mean, data = best_SC)
     
The Prediction is based on the linear model and the Ct values from the unkown Samples. From the Ct's we can now derive 
the amount of oocysts that was in the original Samples (to be exact, we get information on Oocyst equivalents, as Ileum Tissue
does not contain oocysts, but other parasite stages only).
              
		Oocyst_Predict <- 2^predict(linear_model, newdata = unknown_Samples)
              
Some Samples were positive, and for those the prediction model works just fine, but for a majority, there was no fluorescence signal
detected during qPCR, and since they are not NA (aka no qPCR performed), they received a Ct value of "0".
A Ct of Zero does however depict a wrong impression, as this would imply an incredibly high Oocyst amount.
To irradicate that misinterpretation, we now replace the Oocyst Prediction for negative Samples with "0".
A Ct value of 0 leads to an Oocyst_Predict of 4292821751815.77 Oocysts, so we set that number to 0 and turn it into a numeric type
              
		Crypto_Detection <- Crypto_Detection %>% mutate(Oocyst_Predict = replace(Oocyst_Predict, Oocyst_Predict == "4292821751815.77", "0"))
		Crypto_Detection$Oocyst_Predict <- as.numeric(Crypto_Detection$Oocyst_Predict)
              
              
		Oocyst_Predict  == takes values 0:974580 [Oocysts] in analyzed Samples so far 
                
#### DNA EXTRACTION DATA
Now this is just info for the storage of my Samples. Should I ever need more positive or negative Samples, I can simply look into the amount of DNA
I extracted from individual Samples and can re-dilute from the stock. This comes in handy when testing primers in PCRs.
      


