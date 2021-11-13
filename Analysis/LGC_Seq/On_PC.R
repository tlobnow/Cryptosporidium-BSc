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
library(cowplot)


###############################################################################
###############################################################################
# 1. COWP


tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%202.0/COWP-Tree.phy_phyml_tree.txt")


bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 50, "", tree$node.label))


COWP_HI <- HMHZ %>% filter(!is.na(COWP_Ssp)) %>% select(Mouse_ID, HI)

ggtree(tree) %<+% bs_tibble %<+% COWP_HI +
  geom_text(aes(label=bootstrap), hjust=2, nudge_y = 1, size = 3) +
  geom_tiplab(aes(label=label)) +
  geom_text(aes(label = HI))


p <-ggtree(tree) %<+% bs_tibble %<+% COWP_HI +
  geom_tiplab(fontsize = 2) +
  xlim(NA, 0.07) +
  geom_treescale() +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)
#  geom_text(aes(label = HI), hjust=-5, size = 3)

p

p <- p %>% flip(55,56)  +
  #geom_text(aes(label=node), nudge_x = -0.002) +
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
  ggtitle("COWP Maximum-Likelihood Tree") +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        plot.subtitle = element_text(hjust = 0.5, size = 20))
  
p

#save_plot(p, filename = "COWP-ML-Tree.jpg", base_height = 20, base_width = 30)




###############################################################################
###############################################################################
# 3. GP60


tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%202.0/GP60-Tree.phy_phyml_tree.txt")


bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 50, "", tree$node.label))

ggtree(tree)%<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=1, nudge_y = 0.5, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab() +
  xlim(NA,0.59) +
  geom_treescale()  +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)

p <- p %>% flip(81,84) +
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
  ggtitle("GP60 Maximum-Likelihood Tree", subtitle = "Genetic Distance Tamura-Nei, C.parvum Outgroup") +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        plot.subtitle = element_text(hjust = 0.5, size = 20))
  

p

save_plot(p, filename = "GP60-ML-Tree.jpg", base_height = 10, base_width = 15)



###############################################################################
###############################################################################
# 4. MEDLE


tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%202.0/MEDLE-Tree.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 0.9, "", tree$node.label))

ggtree(tree)%<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=1, nudge_y = 0.5, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab() +
  xlim(NA,0.4) +
  geom_treescale()  +
  geom_text(aes(label=bootstrap), hjust=1.3, nudge_y = 0.5, size = 3)

p <- p +
  #geom_text(aes(label=node), nudge_x = -0.002) +
  geom_strip('Tyz-MEDLE', 'AA_0209', barsize=2, color='blue', 
             label="IXb", offset.text=.01, offset = -0.1, extend = 0.1) +
  geom_strip('Tyz-MEDLE', 'AA_0209', barsize=2, color='purple', 
             label="C.tyzzeri", offset.text=.01, offset = .06, extend = 0.1) +
  geom_strip('Hom-MEDLE', 'Par-MEDLE', barsize=2, color='darkgreen', 
             label=
               "C.hominis
C.meleagridis
C.parvum", offset.text=.01, offset = .06, extend = 0.2, fontsize = 4) +
  ggtitle("MEDLE Maximum-Likelihood Tree", subtitle = "Genetic Distance Tamura-Nei, C.parvum Outgroup") +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        plot.subtitle = element_text(hjust = 0.5, size = 20))

p

#save_plot(p, filename = "MEDLE-ML-Tree.jpg", base_height = 20, base_width = 30)


###############################################################################
###############################################################################
# 4. GST


tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%202.0/GST-Tree.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 0.9, "", tree$node.label))

ggtree(tree)%<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=1, nudge_y = 0.1, nudge_x = -0.001, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab() +
  xlim(NA,0.2) +
  geom_treescale()  +
  geom_text(aes(label=bootstrap), hjust=1.1, nudge_y = 0.11, size = 3)

p <- p %>% flip(8,9)+
  #geom_text(aes(label=node), nudge_x = -0.002) +
  geom_strip('Tyz-GST', 'AA_0282', barsize=2, color='blue', 
             label="IXb.1", offset.text=.003, offset = -0.05, extend = 0.2) +
  geom_strip('AA_0209', 'AA_0144', barsize=2, color='blue', 
             label="IXb.2", offset.text=.003, offset = -0.05, extend = 0.2) +
  geom_strip('Tyz-GST', 'AA_0144', barsize=2, color='purple', 
             label="C.tyzzeri", offset.text=.003, offset = .02, extend = 0.2) +
  geom_strip('Hom-GST', 'Mel-GST', barsize=2, color='darkgreen', 
             label=
               "C.hominis
C.parvum
C.meleagridis", offset.text=.003, offset = .02, extend = 0.3, fontsize = 4) +
  ggtitle("GST Maximum-Likelihood Tree", subtitle = "Genetic Distance Tamura-Nei, C.parvum Outgroup") +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        plot.subtitle = element_text(hjust = 0.5, size = 20))



p

#save_plot(p, filename = "GST-ML-Tree.jpg", base_height = 20, base_width = 30)








###############################################################################
###############################################################################
# 4. CP56


tree <- read.tree("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%202.0/CP56-Tree.phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 0.9, "", tree$node.label))

ggtree(tree)%<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=1, nudge_y = 0.1, nudge_x = -0.001, size = 3) +
  geom_tiplab(aes(label=label))


p <-ggtree(tree) %<+% bs_tibble +
  geom_tiplab() +
  xlim(NA,1.2) +
  geom_treescale()  +
  geom_text(aes(label=bootstrap), hjust=1.1, nudge_y = 0.11, size = 3)
p

p <- p +
  #geom_text(aes(label=node), nudge_x = -0.002) +
  geom_strip('Tyz-CP56', 'AA_0282', barsize=2, color='blue', 
             label="CP1", offset.text=.01, offset = 0.08, extend = 0.3) +
  geom_strip('AA_0209', 'AA_0144', barsize=2, color='blue', 
             label="CP2", offset.text=.01, offset = 0.08, extend = 0.3) +
  geom_strip('Tyz-CP56', 'AA_0144', barsize=2, color='purple', 
             label="C.tyzzeri", offset.text=.01, offset = .15, extend = 0.3) +
  geom_strip('Hom-CP56', 'Mel-CP56', barsize=2, color='darkgreen', 
             label=
               "C.hominis
C.parvum
C.meleagridis", offset.text=.01, offset = .15, extend = 0.3, fontsize = 4) +
  ggtitle("CP56 Maximum-Likelihood Tree", subtitle = "Genetic Distance Tamura-Nei, C.parvum Outgroup") +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        plot.subtitle = element_text(hjust = 0.5, size = 20))



p

save_plot(p, filename = "CP56-ML-Tree.jpg", base_height = 20, base_width = 30)

























###############################################################################
###############################################################################
# 2. TRAP-C1


tree <- read.nexus('https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%202.0/TRAP-C1-Tree.nex')

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
###############################################################################
# 5. MSC6-7


tree <- read.nexus('https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%202.0/MSC6-7-Tree.nex')

p <- ggtree(tree) + 
  geom_tiplab()
p

###############################################################################
###############################################################################
# 7. SKSR


tree <- read.nexus('https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Trees%202.0/SKSR-Tree.nex')

p <- ggtree(tree) + 
  geom_tiplab()
p

