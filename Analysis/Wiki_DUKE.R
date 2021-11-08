library(phangorn)
library(ape)


GP60  = read.phyDat("~/Documents/Programming/R/HZ_SC_and_Raw_Data/Cryptosporidium-BSc/Analysis/Trees 2.0/GP60-Tree.phy")
dm    = dist.dna(as.DNAbin(GP60))

treeUPGMA = upgma(dm)
treeNJ    = NJ(dm)

layout(matrix(c(1,2)), height=c(2,2.5))

plot(treeUPGMA, main = "UPGMA", cex = 0.2)
plot(treeNJ, "unrooted", main = "NJ", cex = 0.2)

parsimony(treeUPGMA, GP60)  # 317 <-- most parsimonious tree, has lowest score
parsimony(treeNJ, GP60)     # 322

optParsUPGMA = optim.parsimony(treeUPGMA, GP60)
optParsNJ = optim.parsimony(treeNJ, GP60)

plot(optParsUPGMA, main="UPGMA", cex = 0.8) # rooted tree on top
plot(optParsNJ, "unrooted", main="NJ", cex = 0.5) # unrooted tree on bottom

#### ML

fit_treeUPGMA = pml(unroot(treeUPGMA), data = GP60)