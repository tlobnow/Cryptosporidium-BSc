library(phangorn)
library(ape)
library(adegenet)
library(stats)
library(ade4)
library(DECIPHER)
library(tidyverse)
library(data.table)

##### MEDLE

MEDLE  = read.phyDat("~/Cryptosporidium-BSc/Analysis/Trees 2.0/MEDLE-Tree.phy")
dm    = dist.dna(as.DNAbin(MEDLE))

treeUPGMA = upgma(dm)
treeNJ    = NJ(dm)

layout(matrix(c(1,2)), height=c(2,2.5))

plot(treeUPGMA, main = "UPGMA", cex = 1)
plot(treeNJ, "unrooted", main = "NJ", cex = 1)

parsimony(treeUPGMA, MEDLE)  # 317 <-- most parsimonious tree, has lowest score
parsimony(treeNJ, MEDLE)     # 322

optParsUPGMA = optim.parsimony(treeUPGMA, MEDLE)
optParsNJ = optim.parsimony(treeNJ, MEDLE)

plot(optParsUPGMA, main="UPGMA", cex = 1) # rooted tree on top
plot(optParsNJ, "unrooted", main="NJ", cex = 1) # unrooted tree on bottom


fit_treeUPGMA = pml(unroot(treeUPGMA), data = MEDLE)

fit_treeUPGMA_opt1 = optim.pml(fit_treeUPGMA)
layout(matrix(c(1,2)), height=c(1,1))
par(mar = c(.1,.1,.1,.1))
plot(fit_treeUPGMA, main="default branches", cex = 0.8)   # top = default branch lengths
plot(fit_treeUPGMA_opt1, main="optimized branches", cex = 0.8)   # bottom = optimized branch lengths

# In the graphics window, the tree on top shows the default branches, while the tree on the bottom shows the branches after optimization. 
# Longer branches indicate that more molecular change has occurred along a given branch. 
# You will find that the likelihood of the data under optimized branch lengths has increased considerably - from 1573.4 to 1555. 
# We can compare these models using the 'AIC' function and we find unsurprisingly that the AIC is substantially lower for the tree with branches optimized.

AIC(fit_treeUPGMA, fit_treeUPGMA_opt1)

#We can also search for a better tree by re-arranging the branches, i.e., by setting the parameter 'optNni' to 'TRUE', which causes the function 'optim.pml' to optimize tree topology in addition to branch lengths:
  
  
fit_treeUPGMA_opt2 = optim.pml(fit_treeUPGMA, optNni=TRUE)

layout(matrix(c(1,2)), height=c(1,1))

plot(fit_treeUPGMA_opt1, cex = 0.8)  # top = original topology with optimized branch lengths
plot(fit_treeUPGMA_opt2, cex = 0.8)    # bottom = optimized topology AND branch lengths

AIC(fit_treeUPGMA_opt1, fit_treeUPGMA_opt2)


##### GST
GST  = read.phyDat("~/Cryptosporidium-BSc/Analysis/Trees 2.0/GST-Tree.phy")
dm    = dist.dna(as.DNAbin(GST))

treeUPGMA = upgma(dm)
treeNJ    = NJ(dm)

layout(matrix(c(1,2)), height=c(2,2.5))

plot(treeUPGMA, main = "UPGMA", cex = 1)
plot(treeNJ, "unrooted", main = "NJ", cex = 1)

parsimony(treeUPGMA, GST)  # 317 <-- most parsimonious tree, has lowest score
parsimony(treeNJ, GST)     # 322

optParsUPGMA = optim.parsimony(treeUPGMA, GST)
optParsNJ = optim.parsimony(treeNJ, GST)

plot(optParsUPGMA, main="UPGMA", cex = 1) # rooted tree on top
plot(optParsNJ, "unrooted", main="NJ", cex = 1) # unrooted tree on bottom


fit_treeUPGMA = pml(unroot(treeUPGMA), data = GST)

fit_treeUPGMA_opt1 = optim.pml(fit_treeUPGMA)
layout(matrix(c(1,2)), height=c(1,1))
par(mar = c(.1,.1,.1,.1))
plot(fit_treeUPGMA, main="default branches", cex = 0.8)   # top = default branch lengths
plot(fit_treeUPGMA_opt1, main="optimized branches", cex = 0.8)   # bottom = optimized branch lengths

# In the graphics window, the tree on top shows the default branches, while the tree on the bottom shows the branches after optimization. 
# Longer branches indicate that more molecular change has occurred along a given branch. 
# You will find that the likelihood of the data under optimized branch lengths has increased considerably - from 1573.4 to 1555. 
# We can compare these models using the 'AIC' function and we find unsurprisingly that the AIC is substantially lower for the tree with branches optimized.

AIC(fit_treeUPGMA, fit_treeUPGMA_opt1)

#We can also search for a better tree by re-arranging the branches, i.e., by setting the parameter 'optNni' to 'TRUE', which causes the function 'optim.pml' to optimize tree topology in addition to branch lengths:


fit_treeUPGMA_opt2 = optim.pml(fit_treeUPGMA, optNni=TRUE)

layout(matrix(c(1,2)), height=c(1,1))

plot(fit_treeUPGMA_opt1, cex = 0.8)  # top = original topology with optimized branch lengths
plot(fit_treeUPGMA_opt2, cex = 0.8)    # bottom = optimized topology AND branch lengths

AIC(fit_treeUPGMA_opt1, fit_treeUPGMA_opt2)


##### GP60
GP60  = read.phyDat("~/Cryptosporidium-BSc/Analysis/Trees 2.0/GP60-Tree.phy")
dm    = dist.dna(as.DNAbin(GP60))

treeUPGMA = upgma(dm)
treeNJ    = NJ(dm)

layout(matrix(c(1,2)), height=c(2,2.5))

plot(treeUPGMA, main = "UPGMA", cex = 1)
plot(treeNJ, "unrooted", main = "NJ", cex = 1)

parsimony(treeUPGMA, GP60)  # 317 <-- most parsimonious tree, has lowest score
parsimony(treeNJ, GP60)     # 322

optParsUPGMA = optim.parsimony(treeUPGMA, GP60)
optParsNJ = optim.parsimony(treeNJ, GP60)

plot(optParsUPGMA, main="UPGMA", cex = 1) # rooted tree on top
plot(optParsNJ, "unrooted", main="NJ", cex = 1) # unrooted tree on bottom


fit_treeUPGMA = pml(unroot(treeUPGMA), data = GP60)

fit_treeUPGMA_opt1 = optim.pml(fit_treeUPGMA)
layout(matrix(c(1,2)), height=c(1,1))
par(mar = c(.1,.1,.1,.1))
plot(fit_treeUPGMA, main="default branches", cex = 0.8)   # top = default branch lengths
plot(fit_treeUPGMA_opt1, main="optimized branches", cex = 0.8)   # bottom = optimized branch lengths

# In the graphics window, the tree on top shows the default branches, while the tree on the bottom shows the branches after optimization. 
# Longer branches indicate that more molecular change has occurred along a given branch. 
# You will find that the likelihood of the data under optimized branch lengths has increased considerably - from 1573.4 to 1555. 
# We can compare these models using the 'AIC' function and we find unsurprisingly that the AIC is substantially lower for the tree with branches optimized.

AIC(fit_treeUPGMA, fit_treeUPGMA_opt1)

#We can also search for a better tree by re-arranging the branches, i.e., by setting the parameter 'optNni' to 'TRUE', which causes the function 'optim.pml' to optimize tree topology in addition to branch lengths:


fit_treeUPGMA_opt2 = optim.pml(fit_treeUPGMA, optNni=TRUE)

layout(matrix(c(1,2)), height=c(1,1))

plot(fit_treeUPGMA_opt1, cex = 0.8)  # top = original topology with optimized branch lengths
plot(fit_treeUPGMA_opt2, cex = 0.8)    # bottom = optimized topology AND branch lengths

AIC(fit_treeUPGMA_opt1, fit_treeUPGMA_opt2)

getwd()

##### MSC6-7

