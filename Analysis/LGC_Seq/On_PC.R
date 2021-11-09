#install CRAN Task View for phylogenetics install.packages('ctv')
#install.views('Phylogenetics')
#update.views('Phylogenetics')


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")

library(ctv)
library(ape)
library(ggtree)
library(phytools)
library(ggtree)
library(tidyverse)
library(tibble)


###############################################################################
###############################################################################
# 1. COWP


tree <- read.tree("~/Cryptosporidium-BSc/Analysis/Trees 2.0/COWP-Tree.phy_phyml_tree.txt")


bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 0.01, "", tree$node.label))

ggtree(tree) %<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=2, nudge_y = 1, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab() +
  xlim(NA, 0.1) +
  geom_treescale() +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)


p %>% flip(55,56)  +
  geom_strip("NZ_1642", "AA_0585", barsize=2, color='blue', 
             label="C 2", offset.text=.001, offset = 0.01, extend = 1) +
  geom_strip('CR_2084', 'CR_2128', barsize=2, color='red', 
             label="C 1", offset.text=.001, offset = 0.01, extend = 0.5) +
  geom_strip('Par-COWP', 'Mel-COWP', barsize=2, color='darkgreen', 
             label=
               "C.parvum
C.hominis
C.meleagridis", offset.text=.001, offset = .02, fontsize = 3) +
  geom_strip('Tyz-COWP', 'CR_2128', barsize=2, color='purple', 
             label="C.tyzzeri", offset.text=.001, offset = .02, extend = 0.5) +
  ggtitle("COWP Neighbor-Joining Tree", subtitle = "Genetic Distance Tamura-Nei, C.parvum Outgroup") +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        plot.subtitle = element_text(hjust = 0.5, size = 20))
  




###############################################################################
###############################################################################
# 3. GP60


tree <- read.tree("~/Cryptosporidium-BSc/Analysis/Trees 2.0/GP60-Tree.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 0.9, "", tree$node.label))

ggtree(tree)%<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=1, nudge_y = 0.5, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab() +
  xlim(NA,0.59) +
  geom_treescale()  +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)

p %>% flip(81,84) +
  #geom_text(aes(label=node), nudge_x = -0.002) +
  geom_strip('AA_0523', 'Tyz-GP60', barsize=2, color='blue', 
             label="IXb", offset.text=.01, offset = 0.07, extend = 0.5) +
  geom_strip('CR_2085', 'CR_2206', barsize=2, color='red', 
             label="IXa", offset.text=.01, offset = 0.07) +
  geom_strip('AA_0523', 'CR_2206', barsize=2, color='purple', 
             label="C.tyzzeri", offset.text=.01, offset = .1, extend = 0.6) +
  geom_strip('AA_0523', 'G_2136', barsize = 2, color = "blue",
             label = 'IXb.1', offset.text = .006, offset = 0.03, extend = 0.5) +
  geom_strip('CR_2163', 'G_3224', barsize = 2, color = "blue",
             label = 'IXb.2', offset.text = .006, offset = 0.03, extend = 0.4) +
  geom_strip('AA_0144', 'NZ_1644', barsize = 2, color = "blue",
             label = 'IXb.3', offset.text = .006, offset = 0.03, extend = 0.8) +
  geom_strip('Hom-GP60', 'Par-GP60', barsize=2, color='darkgreen', 
             label=
"C.meleagridis
C.hominis
C.parvum", offset.text=.01, offset = .1, extend = 1, fontsize = 2.5) +
  ggtitle("GP60 Neighbor-Joining Tree", subtitle = "Genetic Distance Tamura-Nei, C.parvum Outgroup") +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        plot.subtitle = element_text(hjust = 0.5, size = 20))
  




###############################################################################
###############################################################################
# 4. MEDLE


tree <- read.tree("~/Cryptosporidium-BSc/Analysis/Trees 2.0/MEDLE-Tree.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 0.9, "", tree$node.label))

ggtree(tree)%<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=1, nudge_y = 0.5, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab() +
  xlim(NA,0.8) +
  geom_treescale()  +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)

p +
  #geom_text(aes(label=node), nudge_x = -0.002) +
  geom_strip('Tyz-MEDLE', 'AA_0209', barsize=2, color='blue', 
             label="IXb", offset.text=.01, offset = -0.07, extend = 0.1) +
  geom_strip('Tyz-MEDLE', 'AA_0209', barsize=2, color='purple', 
             label="C.tyzzeri", offset.text=.01, offset = .1, extend = 0.1) +
  geom_strip('Hom-MEDLE', 'Par-MEDLE', barsize=2, color='darkgreen', 
             label=
               "C.hominis
C.meleagridis
C.parvum", offset.text=.01, offset = .1, extend = 0.2, fontsize = 4) +
  ggtitle("MEDLE Neighbor-Joining Tree", subtitle = "Genetic Distance Tamura-Nei, C.parvum Outgroup") +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        plot.subtitle = element_text(hjust = 0.5, size = 20))


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

