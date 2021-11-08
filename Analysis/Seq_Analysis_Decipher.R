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
                    "./CRYPTO_GP60.sqlite") # path to new database file

# Import sequences from a GenBank formatted file
Seqs2DB(paste(data_dir,
              "/GP60_Trim.fasta", sep = ""),
        type="FASTA",
        dbFile=dbConn,
        identifier = "GP60"
        )

# View the database table that was constructed
BrowseDB(dbConn)

# Retrieve the imported sequences
dna <- SearchDB(dbConn)
dna


# View the database table that was constructed
BrowseDB(dbConn)

# Retrieve the imported sequences
dna <- SearchDB(dbConn)
dna <- AlignSeqs(dna, iterations = 3)


# Display the sequences in a web browser
BrowseSeqs(dna)
# show differences with the first sequence
BrowseSeqs(dna, highlight=1)
# show differences with the consensus sequence
BrowseSeqs(dna, highlight=0)
# change the degree of consensus
BrowseSeqs(dna, highlight=0, threshold=0.2)

# note the pattern common to most sequences
pattern <- DNAStringSet("GCT")
BrowseSeqs(dna,
           patterns=pattern)








# Align the sequences based on their translations
#DNA <- AlignTranslation(dna)
#DNA

# Display the sequences in a web browser
#BrowseSeqs(DNA)
# show differences with the first sequence
#BrowseSeqs(DNA, highlight=1)
# show differences with the consensus sequence
#BrowseSeqs(DNA, highlight=0)
# change the degree of consensus
#BrowseSeqs(DNA, highlight=0, threshold=0.2)

# note the pattern common to most sequences
#pattern <- DNAStringSet("TAGATTTAGCWATTTTTAGTTTACA")
#BrowseSeqs(DNA,
#           patterns=pattern)

# The protein sequences are very similar
#AA <- AlignTranslation(dna, type="AAStringSet")
#BrowseSeqs(AA, highlight=1)

# Choose a reference for frameshift correction
#REF <- translate(dna[11]) # sequence #11

# Correct the frameshift in sequence #12
#correct <- CorrectFrameshifts(myXStringSet=dna[12],
#                              myAAStringSet=REF,
#                              type="both")
#correct
#dna[12] <- correct$sequence

# Sequence #11 is now identical to #12
DNA <- dna
BrowseSeqs(DNA)

# Identify clusters for primer design
d <- DistanceMatrix(DNA)
dim(d) # a symmetric matrix
c <- IdClusters(d,
                method="NJ",
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
                            minProductSize=250,
                            maxProductSize=300,
                            annealingTemp=55,
                            maxPermutations=8)
primers[1,] # the top scoring primer set

# Highlight the primers' target sites
BrowseSeqs(DNA,
           patterns=c(DNAStringSet(primers[1, 1]),
                      reverseComplement(DNAStringSet(primers[1, 2]))))

################################################################################
################################################################################
library(DECIPHER)

data_dir <- "/Users/FinnLo/Documents/Programming/R/HZ_SC_and_Raw_Data/Cryptosporidium-BSc/Analysis/LGC_Seq"


# Create a connection to an on-disk SQLite database
dbConn <- dbConnect(SQLite(), 
                    "./CRYPTO_GP60.sqlite") # path to new database file

# Import sequences from a GenBank formatted file
Seqs2DB(paste(data_dir,
              "/GP60_Trim.fasta", sep = ""),
        type="FASTA",
        dbFile=dbConn,
        identifier = "GP60")

# identify the sequences by their species
x <- dbGetQuery(dbConn,
                "select description from Seqs")$description

Add2DB(myData=data.frame(identifier=x,
                         stringsAsFactors=FALSE),
       dbConn)

dna <- SearchDB(dbConn)
dna <- AlignSeqs(dna, iterations = 3)

# form a consensus for each species
cons <- IdConsensus(dbConn,
                    threshold=0.3,
                    minInformation=0.1)

# calculate a maximum likelihood tree
d <- DistanceMatrix(dna, 
                    correction="Jukes-Cantor")

dend <- IdClusters(d,
                   method="NJ",
                   type="dendrogram",
                   myXStringSet = cons)


dend <- dendrapply(dend,
                   FUN=function(n) {
                     if(is.leaf(n)) 
                       attr(n, "label") <- 
                         as.expression(substitute(italic(leaf),
                                                  list(leaf=attr(n, "label"))))
                     n
                   })

# display the phylogenetic tree
p <- par(mar=c(1, 1, 1, 10),
         xpd=TRUE)
plot(dend,
     yaxt="n",
     horiz=TRUE)
arrows(-0.1, 6, -0.2, 6,
       angle=90,
       length=0.05,
       code=3)
text(-0.15, 6,
     "0.1",
     adj=c(0.5, -0.5))
par(p)

dbDisconnect(dbConn)

