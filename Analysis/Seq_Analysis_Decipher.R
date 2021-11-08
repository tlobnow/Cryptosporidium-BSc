## ANALYZING SEQUENCES WITH DECIPHER AND OTHER TOOLS IN R

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")


library(tidyverse)
library(data.table)
library(BiocManager)
library(DECIPHER)



# WORKFLOW
# all paths are relative to the installed datasets
data_dir <- "/Users/FinnLo/Documents/Programming/R/HZ_SC_and_Raw_Data/Cryptosporidium-BSc/Analysis/LGC_Seq"


# Create a connection to an on-disk SQLite database
dbConn <- dbConnect(SQLite(), 
                    "./CRYPTO1.sqlite") # path to new database file

# Import sequences from a GenBank formatted file
Seqs2DB(paste(data_dir,
              "/MSC6-7_seqs.fasta", sep = ""),
        type="FASTA",
        dbFile=dbConn,
        identifier = "MSC6-7")

# View the database table that was constructed
BrowseDB(dbConn)

# Retrieve the imported sequences
dna <- SearchDB(dbConn)
dna

# Align the sequences based on their translations
DNA <- AlignTranslation(dna)
DNA

# Display the sequences in a web browser
BrowseSeqs(DNA)
# show differences with the first sequence
BrowseSeqs(DNA, highlight=1)
# show differences with the consensus sequence
BrowseSeqs(DNA, highlight=0)
# change the degree of consensus
BrowseSeqs(DNA, highlight=0, threshold=0.2)

# note the pattern common to most sequences
pattern <- DNAStringSet("TAGATTTAGCWATTTTTAGTTTACA")
BrowseSeqs(DNA,
           patterns=pattern)

# The protein sequences are very similar
AA <- AlignTranslation(dna, type="AAStringSet")
BrowseSeqs(AA, highlight=1)

# Choose a reference for frameshift correction
REF <- translate(dna[11]) # sequence #11

# Correct the frameshift in sequence #12
correct <- CorrectFrameshifts(myXStringSet=dna[12],
                              myAAStringSet=REF,
                              type="both")
correct
dna[12] <- correct$sequence

# Sequence #11 is now identical to #12
DNA <- AlignTranslation(dna)
BrowseSeqs(DNA, highlight=11)

# Identify clusters for primer design
d <- DistanceMatrix(DNA)
dim(d) # a symmetric matrix
c <- IdClusters(d,
                method="UPGMA",
                cutoff=0.05,
                show=TRUE)
head(c) # cluster numbers

# Identify sequences by cluster name in the database
Add2DB(data.frame(identifier=paste("cluster",
                                   c$cluster,
                                   sep="")),
       dbConn)
BrowseDB(dbConn)

# Design primers for next-generation sequencing
primers <- DesignSignatures(dbConn,
                            type="sequence",
                            resolution=5,
                            levels=5,
                            minProductSize=400,
                            maxProductSize=800,
                            annealingTemp=55,
                            maxPermutations=8)
primers[1,] # the top scoring primer set

# Highlight the primers' target sites
BrowseSeqs(DNA,
           patterns=c(DNAStringSet(primers[1, 1]),
                      reverseComplement(DNAStringSet(primers[1, 2]))))
