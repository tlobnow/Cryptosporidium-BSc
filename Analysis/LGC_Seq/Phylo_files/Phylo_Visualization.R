library(ggtree)
library(ape)
library(ctv)
library(ggimage)
library(shadowtext)
#install.views('Phylogenetics')
#update.views('Phylogenetics')

library(phylobase)
library(phylotools)
library("treeio")
library("ggtree")

xfun::session_info(
  c('ape', 'aplot', 'ggplot2', 'ggtree', 'ggtreeExtra', 'tidytree', 'treeio'),
  dependencies = FALSE)

data_dir <- "/Users/FinnLo/Documents/Programming/R/HZ_SC_and_Raw_Data/Cryptosporidium-BSc/Analysis/LGC_Seq/Phylo_files"
setwd("/Users/FinnLo/Documents/Programming/R/HZ_SC_and_Raw_Data/Cryptosporidium-BSc/Analysis/LGC_Seq/Phylo_files/")

tree <- read.nexus("MSC6-7 Tree.nex")

p <- ggtree(tree, layout = "rectangular") + 
  geom_tiplab(size = 3, align = F, hjust = -0.1) +
  xlim(NA, 0.4) 

p + 
  geom_strip('AA_0900', 'AA_0144',  barsize=2, color='red', 
             label="Clade 1", offset = -0.1, offset.text = 0.005, extend = 0.3) + 
  geom_strip('AA_0523', 'AA_0282', barsize=2, color='blue', 
             label = "Clade 2", offset = -0.1, offset.text = 0.005, extend = 0.3) +
  geom_strip('AA_0900', 'AA_0282', barsize=2, color='black', 
             label = "C.tyzzeri", offset= 0.06, offset.text = 0.005) +
  geom_strip('Hom-MSC6-7', 'Hom-MSC6-7', barsize=2, color='black', 
             label = "C.hominis", offset= 0.06, offset.text = 0.005, extend = 0.5) +
  geom_strip('Par-MSC6-7', 'AA_0325', barsize=2, color='black', 
             label = "C.parvum", offset= 0.06, offset.text = 0.005, extend = 0.2) +
  geom_treescale() +
  ggtitle("MSC6-7 Tree (NJ-Tree)") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

################################################################################

## GP60

tree <- read.nexus("GP60 Tree.nex")


p <- ggtree(tree, layout = "rectangular") + 
  geom_tiplab(size = 3, align = T) +
  xlim(NA, 20)


p + 
  geom_strip(taxa1 = 'GP60-AA_523',  taxa2 = 'GP60-AA_0900', barsize=2, color='blue', 
             label = "Clade 1", offset = 0, offset.text = 0, extend = 0.1) +
  geom_treescale() +
  ggtitle(" Tree (NJ-Tree)") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  


