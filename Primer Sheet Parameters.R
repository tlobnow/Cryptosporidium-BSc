library(ggplot2)
library(dplyr)



Spreadsheet <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/master/Primer%20Spreadsheet%20-%20Primer%20Sheet%20no%20merged%20cells.csv")
Spreadsheet


Primer_Pair           <- Spreadsheet$Primer.Pair
Locus                 <- Spreadsheet$Locus
Locus_Direction       <- Spreadsheet$Locus.Direction
Gene_ID               <- Spreadsheet$Gene.ID
CDS_ID                <- Spreadsheet$CDS.ID
Chr                   <- Spreadsheet$Chr.
Name_Primer           <- Spreadsheet$Name
Size                  <- Spreadsheet$Size.including.Primers
Primer_Length         <- Spreadsheet$Primer.Length
Min_Pos               <- Spreadsheet$Min.binding.position.in.genome
Max_Pos               <- Spreadsheet$Max.binding.position.in.genome
Gen_Seq               <- Spreadsheet$Sequence
GC_Content            <- Spreadsheet$X.GC
Tm_Geneious           <- Spreadsheet$Tm
Tm_basic              <- Spreadsheet$Basic.Tm..C.
Tm_Salt_adjust        <- Spreadsheet$Salt.Adjusted.Tm..C.
Tm_Nearest_Neighbor   <- Spreadsheet$Nearest.Neighbor.Tm..C.
Hairpin_Geneious      <- Spreadsheet$Hairpin.Tm
Self_Dimer_Geneious   <- Spreadsheet$Self.Dimer.Tm
Pair_Dimer_Geneious   <- Spreadsheet$Pair.Dimer.Tm
Binds_tyzparhom       <- Spreadsheet$Pair.binds.to..1...tyz...2...tyz.par...3...tyz.par.hom...4...tyz.hom
Exon                  <- Spreadsheet$Exon
Target_Seq            <- Spreadsheet$Target..nt.
par_Identical_nt      <- Spreadsheet$Identical.Sites..nt.
Pairwise_nt           <- Spreadsheet$Pairwise.Identity..nt...
OutORF_SNPs           <- Spreadsheet$SNPs.Outside.ORF
InORF_SNPs            <- Spreadsheet$SNPs.in.ORF..nt.
Syn_SNPs              <- Spreadsheet$Synonymous.SNPs..nt.
NonSyn_SNPs           <- Spreadsheet$Non..Synonymous.SNPs..nt.
Indel_InORF           <- Spreadsheet$Indels.in.ORF.Sum..nt.
Indel_OutORF          <- Spreadsheet$Indels.outside.ORF.sum..nt.
Insertion_inORF       <- Spreadsheet$Insertions.in.ORF..nt.
Insertion_OutORF      <- Spreadsheet$Insertions.outside.ORF..nt.
Deletion_inORF        <- Spreadsheet$Deletions.in.ORF..nt.
Deletion_OutORF       <- Spreadsheet$Deletions.outside.ORF..nt.
NNN_unknown_nt        <- Spreadsheet$NNNNN..unknown.nt.
amplified_aa          <- Spreadsheet$X..aa.amplified
par_Identical_aa      <- Spreadsheet$Identical.Sites.with.parvum..aa.
Pairwise_aa           <- Spreadsheet$Pairwise.Identity..nt....1
aa_inserted           <- Spreadsheet$aa.inserted.sum
aa_deleted            <- Spreadsheet$aa.deleted.sum
aa_changes            <- Spreadsheet$aa.changes.sum
aa_X_exchange         <- Spreadsheet$X...aa.exchanges
aa_B_exchange         <- Spreadsheet$B...D.N.exchange
aa_J_exchange         <- Spreadsheet$J...I.L.exchange
aa_Z_exchange         <- Spreadsheet$Z...E.Q.exchange
XXX_unknown_aa        <- Spreadsheet$X......unknown.aa.
Single_Base_Pass      <- Spreadsheet$Single.Base.Runs.Pass
Dinucleotide_Pass     <- Spreadsheet$Dinucleotide.Base.Runs.Pass
Length_Pass           <- Spreadsheet$Length..Pass
Tm_Neighbor_Pass      <- Spreadsheet$Tm.Nearest.Neighbor.Pass
GC_Pass               <- Spreadsheet$X.GC..Pass
GC_Clamp_Pass         <- Spreadsheet$GC...Clamp.Pass
Self_Annealing_Pass   <- Spreadsheet$Self.Annealing.Pass
Hairpin_Pass          <- Spreadsheet$Hairpin.formation..Pass





all_parameters <- data.frame(Primer_Pair ,
                       Locus                ,
                       Locus_Direction      ,
                       Gene_ID             ,
                       CDS_ID                ,
                       Chr                   ,
                       Name_Primer          ,
                       Size                  ,
                       Primer_Length        ,
                       Min_Pos              ,
                       Max_Pos              ,
                       Gen_Seq               ,
                       GC_Content           ,
                       Tm_Geneious          ,
                       Tm_basic             ,
                       Tm_Salt_adjust        ,
                       Tm_Nearest_Neighbor   ,
                       Hairpin_Geneious      ,
                       Self_Dimer_Geneious   ,
                       Pair_Dimer_Geneious  ,
                       Binds_tyzparhom      ,
                       Exon                  ,
                       Target_Seq            ,
                       par_Identical_nt     ,
                       Pairwise_nt          ,
                       OutORF_SNPs           ,
                       InORF_SNPs           ,
                       Syn_SNPs             ,
                       NonSyn_SNPs           ,
                       Indel_OutORF          ,
                       Indel_InORF,
                       Insertion_inORF ,
                       Insertion_OutORF ,
                       Deletion_inORF    ,
                       Deletion_OutORF ,
                       NNN_unknown_nt        ,
                       amplified_aa         ,
                       par_Identical_aa     ,
                       Pairwise_aa          ,
                       aa_inserted           ,
                       aa_deleted            ,
                       aa_changes            ,
                       aa_X_exchange         ,
                       aa_B_exchange         ,
                       aa_J_exchange        ,
                       aa_Z_exchange         ,
                       XXX_unknown_aa        ,
                       Single_Base_Pass     ,
                       Dinucleotide_Pass    ,
                       Length_Pass           ,
                       Tm_Neighbor_Pass      ,
                       GC_Pass               ,
                       GC_Clamp_Pass         ,
                       Self_Annealing_Pass  ,
                       Hairpin_Pass          
                       )
all_parameters

sub_parameters <- filter(all_parameters, Size <= 300, Size >= 250)
sub_parameters

subset_data <- sub_parameters %>% 
                group_by(Size) %>%
                filter(Single_Base_Pass == T) %>%
                filter(Dinucleotide_Pass == T) %>%
                filter(Length_Pass == T) %>%
                filter(GC_Pass == T) %>%
                filter(GC_Clamp_Pass == T) %>%
                filter(Self_Annealing_Pass == T) %>%
                filter(Hairpin_Pass == T) %>%
                filter(Hairpin_Geneious <= 55 | Hairpin_Geneious == "None")
subset_data


plot(GC_Content, Tm_Geneious)
GCModel <- lm(GC_Content ~ Tm_Geneious, data = subset_data)   # that's how you plot a linear regression model, has to be supplied with a formula, allows to predict defense 
GCModel


table(Chr)
table(OutORF_SNPs)
table(InORF_SNPs)
table(InORF_Indels)
table(Locus_Direction)
table(Binds_tyzparhom)



  

ggplot(all_parameters, aes(Chr, Locus)) +
  geom_bar(stat = "summary", fun = "median", col = "black", position = "dodge") 
  
ggplot(all_parameters, aes(Gene_ID)) +
  geom_bar(col = "black", position = "dodge")

ggplot(all_parameters, aes(Chr, Binds_tyzparhom)) +
  geom_bar(stat = "identity", fill = Binds_tyzparhom, position = "dodge") 



table_Locus_Chr <- table(Locus, Chr)
table_Locus_Chr


              
# GP40 ----            
  GP40_Marker <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                     Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                     Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                     NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                     aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                     Hairpin_Pass) %>% 
                  filter(Locus == "gp40")                  
                  GP40_Marker

  
# Actin ----
  Actin_Marker <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                       Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                       Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                       NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                       aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                       Hairpin_Pass) %>% 
                  filter(Locus == "Actin 1" | Locus == "Actin 2" | Locus == "Actin 3"| Locus == "Actin 4"| Locus == "Actin 5"| Locus == "Actin 6")                  
                  Actin_Marker
  
  
# COWP ----
  COWP_Marker <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                       Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                       Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                       NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                       aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                       Hairpin_Pass) %>% 
                  filter(Locus == "COWP 3" | Locus == "COWP 4" | Locus == "COWP 6"| Locus == "COWP 7"| Locus == "COWP 8"| Locus == "COWP 9")                  
                  COWP_Marker
  
# CP47 ----
  CP47_Marker <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                       Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                       Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                       NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                       aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                       Hairpin_Pass) %>% 
                  filter(Locus == "CP47")                  
                  CP47_Marker

  
# CP56 ----
  CP56_Marker <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                       Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                       Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                       NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                       aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                       Hairpin_Pass) %>% 
                  filter(Locus == "CP56")                  
                  CP56_Marker
  

# GST ----
  GST_Marker <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                        Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                        Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                        NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                        aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                         Hairpin_Pass) %>% 
                  filter(Locus == "GST 1" | Locus == "GST 2" | Locus == "GST 3")                  
                  GP40_Marker


                  
# MEDLE ----
  GST_Marker <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                                       Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                                       Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                                       NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                                       aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                                       Hairpin_Pass) %>% 
                  
                  
# Mucin-1 ----
  Mucin1_table <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                       Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                       Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                       NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                       aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                       Hairpin_Pass) %>% 
                    filter(Locus == "Mucin-1")                  
                    Mucin1_table
  
  
# SKSR ----
  
# TRAP-C1 ----
  TRAPC1_table <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                       Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                       Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                       NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                       aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                       Hairpin_Pass) %>% 
                    filter(Locus == "TRAP-C1")                  
                    TRAPC1_table
  
  
  
  
  
  
  
  
  







  
  
  
  
  
  
  
  
  