library(ctv)
library(ape)
library(ggtree)
library(phytools)
library(ggtree)
library(tidyverse)
library(tibble)
library(cowplot)


## GST #########################################################################
tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%203.0/GST.phy_phyml_tree.txt")


bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 50, "", tree$node.label))


ggtree(tree) %<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=2, nudge_y = 1, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab(fontsize = 2) +
  xlim(NA, 0.2) +
  geom_treescale() +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)
p

## COWP #########################################################################
tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%203.0/COWP.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 50, "", tree$node.label))


ggtree(tree) %<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=2, nudge_y = 1, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab(fontsize = 2) +
  xlim(NA, 0.05) +
  geom_treescale() +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)
p



## MEDLE #######################################################################
tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%203.0/MEDLE.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 50, "", tree$node.label))


ggtree(tree) %<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=2, nudge_y = 1, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab(fontsize = 2) +
  xlim(NA, 0.5) +
  geom_treescale() +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)
p

## MEDLE_2 #######################################################################
tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%203.0/MEDLE_2.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 50, "", tree$node.label))


ggtree(tree) %<+% bs_tibble +
  geom_text(aes(label=node), nudge_x = -0.002) +
  geom_text(aes(label=parent), nudge_x =  -0.2) +
  geom_text(aes(label=bootstrap), hjust=2, nudge_y = 1, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  #geom_text(aes(label=parent), nudge_x =  -0.2) +
  geom_tiplab(fontsize = 2) +
  xlim(NA, 0.5) +
  geom_treescale() +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)
p %>% 
  #flip(1,18) %>% 
  flip(19, 20) %>% 
  flip(3,4)

## MEDLE_3 #######################################################################
tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%203.0/MEDLE_3.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 50, "", tree$node.label))


ggtree(tree) %<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=2, nudge_y = 0, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab(fontsize = 2) +
  xlim(NA, 0.5) +
  geom_treescale() +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0, size = 3)
p


## CP56 #########################################################################

tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%203.0/CP56.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 50, "", tree$node.label))


ggtree(tree) %<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=2, nudge_y = 1, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab(fontsize = 2) +
  xlim(NA, 2) +
  geom_treescale() +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)
p

## CP56-2 #########################################################################

tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%203.0/CP56_2.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 50, "", tree$node.label))


ggtree(tree) %<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=2, nudge_y = 1, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab(fontsize = 2) +
  xlim(NA, 2) +
  geom_treescale() +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)
p

## CP56-3 #########################################################################

tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%203.0/CP56_3.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 100, "", tree$node.label))


ggtree(tree) %<+% bs_tibble +
  geom_text(aes(label=node), nudge_x = -0.002) +
  geom_text(aes(label=parent), nudge_x =  -0.2) +
  geom_text(aes(label=bootstrap), hjust=2, nudge_y = 1, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab(fontsize = 2) +
  xlim(NA, 2) +
  geom_treescale() +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)
p

## GP60 #########################################################################
tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%203.0/GP60.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 50, "", tree$node.label))


ggtree(tree) %<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=2, nudge_y = 1, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab(fontsize = 2) +
  xlim(NA, 1) +
  geom_treescale() +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)
p


## MSC6-7 #########################################################################
tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%203.0/MSC6-7.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 50, "", tree$node.label))


ggtree(tree) %<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=1.1, nudge_y = 0, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab(fontsize = 2) +
  xlim(NA, 2) +
  geom_treescale() +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)
p

## SKSR #########################################################################
tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%203.0/SKSR.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 50, "", tree$node.label))


ggtree(tree) %<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=1.1, nudge_y = 0, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab(fontsize = 2) +
  xlim(NA, 2) +
  geom_treescale() +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)
p


## Concatenated #########################################################################
tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%203.0/cc.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 1, "", tree$node.label))


ggtree(tree) %<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=1.1, nudge_y = 0, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab(fontsize = 2) +
  xlim(NA, 0.9) +
  geom_treescale() +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0, size = 3)
p %>% flip(5,15)




