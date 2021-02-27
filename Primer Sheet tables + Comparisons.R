
display.brewer.all()

table(Chr)
table(OutORF_SNPs)
table(InORF_SNPs)
table(InORF_Indels)
table(Locus_Direction)
table(Binds_tyzparhom)
table_Locus_Chr <- table(Locus, Chr)
  table_Locus_Chr

# Loci on different Chr
ggplot(all_parameters, aes(Chr, Locus)) +
  geom_bar(stat = "summary", fun = "median", col = "black", position = "dodge") 

# diversity trends per Chr
ggplot(all_parameters, aes(Chr, Binds_tyzparhom)) +
  geom_bar(stat = "identity", fill = Binds_tyzparhom, position = "dodge") 


# comparing number of InORF_SNPs per Locus and colored by Chromosome
  # point plot
ggplot(all_parameters, aes(InORF_SNPs, Locus, color = Chr)) + 
  geom_point()
  # bar plot
ggplot(all_parameters, aes(InORF_SNPs, Locus, fill = Chr)) +
  geom_bar(stat = "identity", position = "dodge")  +
  scale_color_brewer(palette = "Set2")


plot()





# GP40 ----            
GP40_Marker <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                      Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                      Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                      NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                      aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                      Hairpin_Pass) %>% 
  filter(Locus == "gp40")                  
GP40_Marker

GP40_parameter1 <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, 
                          Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs) %>%
  filter(Locus == "gp40")   

GP40_parameter1
  plot(GP40_parameter1)






GP40_Tibble <- tibble(Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                      Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                      Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                      NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                      aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                      Hairpin_Pass) %>% 
  filter(Locus == "gp40")
GP40_Tibble



# Actin ----
Actin_Marker <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                       Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                       Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                       NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                       aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                       Hairpin_Pass) %>% 
  filter(Locus == "Actin 1" | Locus == "Actin 2" | Locus == "Actin 3"| Locus == "Actin 4"| Locus == "Actin 5"| Locus == "Actin 6")                  
Actin_Marker


Actin_Tibble <- tibble(Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                       Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                       Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                       NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                       aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                       Hairpin_Pass) %>% 
  filter(Locus == "Actin 1" | Locus == "Actin 2" | Locus == "Actin 3"| Locus == "Actin 4"| Locus == "Actin 5"| Locus == "Actin 6")                  
Actin_Tibble


# COWP ----
COWP_Marker <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                      Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                      Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                      NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                      aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                      Hairpin_Pass) %>% 
  filter(Locus == "COWP 3" | Locus == "COWP 4" | Locus == "COWP 6"| Locus == "COWP 7"| Locus == "COWP 8"| Locus == "COWP 9")                  
COWP_Marker
plot(COWP_Marker)



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
GST_Marker



# MEDLE ----
MEDLE_Marker <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                       Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                       Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                       NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                       aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                       Hairpin_Pass) %>% 
  filter(Locus == "MEDLE 1" | Locus == "MEDLE 2")                 
MEDLE_Marker

# Mucin-1 ----
Mucin1_Marker <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                       Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                       Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                       NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                       aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                       Hairpin_Pass) %>% 
  filter(Locus == "Mucin-1")                  
Mucin1_Marker


# SKSR ----
SKSR_Marker <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                     Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                     Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                     NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                     aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                     Hairpin_Pass) %>% 
  filter(Locus == "SKSR 1 + SKSR 2" | Locus == "SKSR 3" | Locus == "SKSR 4" | Locus == "SKSR 5" | Locus == "SKSR 6" | Locus == "SKSR 7" | Locus == "SKSR 8" | Locus == "SKSR 9" | Locus == "SKSR 10")                  
SKSR_Marker      


# TRAP-C1 ----
TRAPC1_Marker <- select(all_parameters, Primer_Pair, Locus, Locus_Direction, Gene_ID, CDS_ID, Chr, Name_Primer, Size, Primer_Length, Min_Pos, Max_Pos, Gen_Seq, GC_Content, Tm_Geneious, Tm_basic, Tm_Salt_adjust,
                       Tm_Nearest_Neighbor, Hairpin_Geneious, Self_Dimer_Geneious, Pair_Dimer_Geneious, Binds_tyzparhom, Exon, Target_Seq, par_Identical_nt, Pairwise_nt, OutORF_SNPs, InORF_SNPs,
                       Syn_SNPs, NonSyn_SNPs, Indel_OutORF, Indel_InORF, Insertion_inORF, Insertion_OutORF, Deletion_inORF, Deletion_OutORF,
                       NNN_unknown_nt, amplified_aa, par_Identical_aa, Pairwise_aa, aa_inserted, aa_deleted, aa_changes, aa_X_exchange, aa_B_exchange,
                       aa_J_exchange, aa_Z_exchange, XXX_unknown_aa, Single_Base_Pass, Dinucleotide_Pass, Length_Pass, Tm_Neighbor_Pass, GC_Pass, GC_Clamp_Pass, Self_Annealing_Pass,
                       Hairpin_Pass) %>% 
  filter(Locus == "TRAP-C1")                  
TRAPC1_Marker


