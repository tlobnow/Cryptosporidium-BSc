#install CRAN Task View for phylogenetics install.packages('ctv')
library('ctv') 
#install.views('Phylogenetics')
#update.views('Phylogenetics')


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

library(ctv)
library(ape)
library(ggtree)
library(phytools)
library(ggtree)
library(tidyverse)



###############################################################################
###############################################################################
# 1. COWP

tree <- read.nexus('~/Trees 2.0/COWP-Tree.nex')


p <- ggtree(tree, layout = "rectangular") + 
  geom_tiplab() +
  xlim(NA, 0.1) +
  geom_treescale()

p +
  geom_hilight(node=74
               , fill="lightgreen", alpha=.6, type = "roundrect") +
  geom_strip(29, 51, barsize=2, color='red', 
             label="C 1", offset.text=.001, offset = -.03, extend = 0.85) +
  geom_strip('G_2136', 'NZ_1642', barsize=2, color='blue', 
             label="C 2", offset.text=.001, offset = -.03) +
  geom_strip('Mel-COWP', 'Par-COWP', barsize=2, color='darkgreen', 
             label=
  "C.meleagridis
C.hominis
C.parvum", offset.text=.001, offset = .01, extend = 1) +
  geom_strip('AA_0559', 'NZ_1642', barsize=2, color='purple', 
             label="C.tyzzeri", offset.text=.001, offset = .01) +
  ggtitle("COWP NJ Tree")

p

### PHY-ML
tree <- read.tree("~/PhyML-3.1/COWP-Tree.phy_phyml_tree.txt")
p <- ggtree(tree) +
  #geom_text(aes(label=node)) +
  geom_tiplab() +
  xlim(NA, 0.1) +
  geom_treescale()
p


###############################################################################
###############################################################################
# 2. TRAP-C1


tree <- read.nexus('~/Trees 2.0/TRAP-C1-Tree.nex')

p <-ggtree(tree) + 
  #geom_text(aes(label=node), nudge_x = -0.002) +
  geom_tiplab(hjust = -0.1) +
  xlim(NA, 0.1) +
  geom_treescale()

 p%>% flip( 16, 17)  %>% flip(42, 56) %>% flip(40,41) + 
  geom_hilight(node=54, fill="lightgreen", alpha=.6, type = "roundrect") +
  geom_strip('NZ_1640', 'AA_0571', barsize=2, color='red', 
             label="T 1", offset.text=.001, offset = -0.055, extend = 0.2) +
  geom_strip('CR_2084', 'AA_0578', barsize=2, color='blue', 
             label="T 2", offset.text=.001, offset = -0.055) +
  geom_strip('Hom-TRAP-C1', 'Par-TRAP-C1', barsize=2, color='darkgreen', 
             label=
               "C.hominis
C.meleagridis
C.parvum", offset.text=.001, offset = .01, extend = 1) +
  geom_strip('AA_0559', 'NZ_1642', barsize=2, color='purple', 
             label="C.tyzzeri", offset.text=.001, offset = .01)




###############################################################################
###############################################################################
# 3. GP60


tree <- read.nexus('~/Trees 2.0/GP60-Tree.nex')


p <- ggtree(tree) +
  #geom_text(aes(label=node)) +
  geom_tiplab() +
  xlim(NA, 1) +
  geom_treescale()
p

p %>% flip(67,69) %>% flip(44,66) %>% flip(64,65) +
  geom_hilight(node=54, fill="lightgreen", alpha=.6, type = "roundrect") +
  geom_strip('NZ_1640', 'AA_0571', barsize=2, color='red', 
             label="T 1", offset.text=.001, offset = -0.055, extend = 0.2) +
  geom_strip('CR_2084', 'AA_0578', barsize=2, color='blue', 
             label="T 2", offset.text=.001, offset = -0.055) +
  geom_strip('Hom-TRAP-C1', 'Par-TRAP-C1', barsize=2, color='darkgreen', 
             label=
               "C.hominis
C.meleagridis
C.parvum", offset.text=.001, offset = .01, extend = 1) +
  geom_strip('AA_0559', 'NZ_1642', barsize=2, color='purple', 
             label="C.tyzzeri", offset.text=.001, offset = .01)



### PHY-ML
tree <- read.tree("~/PhyML-3.1/GP60-Tree.phy_phyml_tree.txt")
p <- ggtree(tree) +
  geom_text(aes(label=node)) +
  geom_tiplab() +
  xlim(NA, 1) +
  geom_treescale()
p


###############################################################################
###############################################################################
# 4. MEDLE


tree <- read.nexus('~/Trees 2.0/MEDLE-Tree.nex')

p <- ggtree(tree) + 
  geom_tiplab()
p


### PHY-ML
tree <- read.tree("~/PhyML-3.1/MEDLE-Tree.phy_phyml_tree.txt")
p <- ggtree(tree) +
  #geom_text(aes(label=node)) +
  geom_tiplab() +
  xlim(NA, 0.3) +
  geom_treescale()
p

###############################################################################
###############################################################################
# 5. MSC6-7


tree <- read.nexus('~/Trees 2.0/MSC6-7-Tree.nex')

p <- ggtree(tree) + 
  geom_tiplab()
p

###############################################################################
###############################################################################
# 7. SKSR


tree <- read.nexus('~/Trees 2.0/SKSR-Tree.nex')

p <- ggtree(tree) + 
  geom_tiplab()
p

###############################################################################
###############################################################################
# 8. GST


tree <- read.nexus('~/Trees 2.0/GST-Tree.nex')

p <- ggtree(tree) + 
  geom_tiplab()
p

tree <- read.tree("~/PhyML-3.1/GST-Tree.phy_phyml_tree.txt")
p <- ggtree(tree) +
 # geom_text(aes(label=node)) +
  geom_tiplab() +
  xlim(NA, 0.1) +
  geom_treescale()
p

