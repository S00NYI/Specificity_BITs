################################################################################
## eCLIP Analysis Plotting V4: 
## Written by Soon Yi
## Created: 2023-12-01
## Last Edited: 2024-02-06
## Figure 3
################################################################################

library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(corrplot)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(GenomicRanges)
library(IRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(ggsignif)

## Custom Functions:
################################################################################
# Min-Max normalization:
min_max_norm = function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Specificity index calculation:
SI_Calculation = function(x) {
  x / (median(x, na.rm = FALSE))
}

Motif_Variants = function(motif, nucleotides) {
  variants = c()
  # Loop over each position in the motif
  for (i in 1:nchar(motif)) {
    # Loop over each nucleotide
    for (nuc in nucleotides) {
      if (substr(motif, i, i) != nuc) {
        # Construct the new variant
        variant = paste0(substr(motif, 1, i - 1), nuc, substr(motif, i + 1, nchar(motif)))
        variants = c(variants, variant)
      }
    }
  }
  return(variants)
}

################################################################################

## Set up basic parameters:
################################################################################
baseDir = '~/Desktop/Genomics/Specificity/Data_Final/'
sampleTable_eCLIP = read_csv(paste0(baseDir, 'sample_table_eCLIP.csv'), col_names = T, show_col_types = F)
sampleTable_eCLIP = data.frame(sampleTable_eCLIP)

cell_lines = c('K562', 'HepG2')

K562_RBPs = as.character(unique((sampleTable_eCLIP %>% filter(Cell == 'K562'))[, 'RBP']))
HepG2_RBPs = as.character(unique((sampleTable_eCLIP %>% filter(Cell == 'HepG2'))[, 'RBP']))

K = 5
NTs = c("A", "C", "G", "T")
KMer_combinations = expand.grid(replicate(K, NTs, simplify = FALSE))
KMer_combinations = apply(KMer_combinations, 1, paste, collapse = "")
KMer_DSS = DNAStringSet(KMer_combinations)
KMer_PDict = PDict(KMer_DSS)

DATE = '20240217'
################################################################################


## RBNS  Stats
################################################################################
RBNS_Stats =  read_csv(paste0(baseDir, '/RBNS/5mer_metrics.csv'), col_names = T, show_col_types = F)

RBNS_5mer = read_csv(paste0(baseDir, 'RBNS/RBNS_normalized_', K, 'mer.csv'), col_names = T, show_col_types = F)
colnames(RBNS_5mer) = c('Motif', colnames(RBNS_5mer)[2:ncol(RBNS_5mer)])
################################################################################

# Compile all eCLIP Data:
################################################################################
K562_RBPs_noStructure = as.character(unique((sampleTable_eCLIP %>% filter(structure_context == 'L') %>% filter(Cell == 'K562'))[, 'RBP']))
HepG2_RBPs_noStructure = as.character(unique((sampleTable_eCLIP %>% filter(structure_context == 'L') %>% filter(Cell == 'HepG2'))[, 'RBP']))
RBPs_noStructure = unique(c(K562_RBPs_noStructure, HepG2_RBPs_noStructure))

cell_line = 'HepG2'
RBP = 'PCBP2'
PCBP2 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
PCBP2 = PCBP2[, -1]
colnames(PCBP2) = c('Motif', 'Score')
PCBP2$Motif = gsub("T", "U", PCBP2$Motif)
PCBP2 = PCBP2 %>% arrange(Motif)

RBP = 'FUBP3'
FUBP3 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
FUBP3 = FUBP3[, -1]
colnames(FUBP3) = c('Motif', 'Score')
FUBP3$Motif = gsub("T", "U", FUBP3$Motif)
FUBP3 = FUBP3 %>% arrange(Motif)

cell_line = 'K562'
RBP = 'IGF2BP2'
IGF2BP2 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
IGF2BP2 = IGF2BP2[, -1]
colnames(IGF2BP2) = c('Motif', 'Score')
IGF2BP2$Motif = gsub("T", "U", IGF2BP2$Motif)
IGF2BP2 = IGF2BP2 %>% arrange(Motif)

eCLIP_noStrucutures = data.frame(MOTIF = PCBP2$Motif,
                                    HNRNPC = NA,
                                    HNRNPL = NA,
                                    PCBP1 = NA,
                                    TIA1 = NA,
                                    KHSRP = NA,
                                    FUS = NA,
                                    IGF2BP2 = IGF2BP2$Score,
                                    FUBP3 = FUBP3$Score,
                                    PCBP2 = PCBP2$Score)


for (RBP in c('HNRNPC', 'HNRNPL', 'PCBP1', 'TIA1', 'KHSRP', 'FUS')) {
  K562 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/K562/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
  K562 = K562[, -1]
  colnames(K562) = c('Motif', 'Score')
  K562$Motif = gsub("T", "U", K562$Motif)
  K562 = K562 %>% arrange(Motif)
  
  HepG2 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/HepG2/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
  HepG2 = HepG2[, -1]
  colnames(HepG2) = c('Motif', 'Score')
  HepG2$Motif = gsub("T", "U", HepG2$Motif)
  HepG2 = HepG2 %>% arrange(Motif)
  
  eCLIP_noStrucutures[, RBP] = min_max_norm((HepG2$Score + K562$Score)/2)
}

K562_RBPs_Structure = as.character(unique((sampleTable_eCLIP %>% filter(structure_context == 'H') %>% filter(Cell == 'K562'))[, 'RBP']))
HepG2_RBPs_Structure = as.character(unique((sampleTable_eCLIP %>% filter(structure_context == 'H') %>% filter(Cell == 'HepG2'))[, 'RBP']))
RBPs_Structure = unique(c(K562_RBPs_Structure, HepG2_RBPs_Structure))

cell_line = 'K562'
RBP = 'HNRNPK_1'
HNRNPK_1 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
HNRNPK_1 = HNRNPK_1[, -1]
colnames(HNRNPK_1) = c('Motif', 'Score')
HNRNPK_1$Motif = gsub("T", "U", HNRNPK_1$Motif)
HNRNPK_1 = HNRNPK_1 %>% arrange(Motif)

RBP = 'HNRNPK_2'
HNRNPK_2 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
HNRNPK_2 = HNRNPK_2[, -1]
colnames(HNRNPK_2) = c('Motif', 'Score')
HNRNPK_2$Motif = gsub("T", "U", HNRNPK_2$Motif)
HNRNPK_2 = HNRNPK_2 %>% arrange(Motif)

HNRNPK = HNRNPK_1
HNRNPK$Score = min_max_norm((HNRNPK_1$Score + HNRNPK_2$Score)/2)

RBP = 'PUM1'
PUM1 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
PUM1 = PUM1[, -1]
colnames(PUM1) = c('Motif', 'Score')
PUM1$Motif = gsub("T", "U", PUM1$Motif)
PUM1 = PUM1 %>% arrange(Motif)

RBP = 'EWSR1'
EWSR1 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
EWSR1 = EWSR1[, -1]
colnames(EWSR1) = c('Motif', 'Score')
EWSR1$Motif = gsub("T", "U", EWSR1$Motif)
EWSR1 = EWSR1 %>% arrange(Motif)

RBP = 'EIF4G2'
EIF4G2 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
EIF4G2 = EIF4G2[, -1]
colnames(EIF4G2) = c('Motif', 'Score')
EIF4G2$Motif = gsub("T", "U", EIF4G2$Motif)
EIF4G2 = EIF4G2 %>% arrange(Motif)

cell_line = 'HepG2'
RBP = 'SRSF9'
SRSF9 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
SRSF9 = SRSF9[, -1]
colnames(SRSF9) = c('Motif', 'Score')
SRSF9$Motif = gsub("T", "U", SRSF9$Motif)
SRSF9 = SRSF9 %>% arrange(Motif)

eCLIP_Strucutures = data.frame(MOTIF = PUM1$Motif,
                                  RBFOX2 = NA,
                                  RBM22 = NA,
                                  TRA2A = NA,
                                  TAF15 = NA,
                                  HNRNPK = HNRNPK$Score,
                                  PUM1 = PUM1$Score,
                                  EWSR1 = EWSR1$Score,
                                  EIF4G2 = EIF4G2$Score,
                                  SRSF9 = SRSF9$Score)


for (RBP in c('RBFOX2', 'RBM22', 'TRA2A', 'TAF15')) {
  K562 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/K562/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
  K562 = K562[, -1]
  colnames(K562) = c('Motif', 'Score')
  K562$Motif = gsub("T", "U", K562$Motif)
  K562 = K562 %>% arrange(Motif)
  
  HepG2 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/HepG2/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
  HepG2 = HepG2[, -1]
  colnames(HepG2) = c('Motif', 'Score')
  HepG2$Motif = gsub("T", "U", HepG2$Motif)
  HepG2 = HepG2 %>% arrange(Motif)
  
  eCLIP_Strucutures[, RBP] = min_max_norm((HepG2$Score + K562$Score)/2)
}


K562_RBPs_unKnown = as.character(unique((sampleTable_eCLIP %>% filter(structure_context == '.') %>% filter(Cell == 'K562'))[, 'RBP']))
HepG2_RBPs_unKnown = as.character(unique((sampleTable_eCLIP %>% filter(structure_context == '.') %>% filter(Cell == 'HepG2'))[, 'RBP']))
RBPs_unKnown = unique(c(K562_RBPs_unKnown, HepG2_RBPs_unKnown))

cell_line = 'K562'
RBP = 'SAFB2'
SAFB2 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
SAFB2 = SAFB2[, -1]
colnames(SAFB2) = c('Motif', 'Score')
SAFB2$Motif = gsub("T", "U", SAFB2$Motif)
SAFB2 = SAFB2 %>% arrange(Motif)

RBP = 'ZRANB2'
ZRANB2 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
ZRANB2 = ZRANB2[, -1]
colnames(ZRANB2) = c('Motif', 'Score')
ZRANB2$Motif = gsub("T", "U", ZRANB2$Motif)
ZRANB2 = ZRANB2 %>% arrange(Motif)

RBP = 'XRCC6'
XRCC6 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
XRCC6 = XRCC6[, -1]
colnames(XRCC6) = c('Motif', 'Score')
XRCC6$Motif = gsub("T", "U", XRCC6$Motif)
XRCC6 = XRCC6 %>% arrange(Motif)

RBP = 'AKAP8L'
AKAP8L = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
AKAP8L = AKAP8L[, -1]
colnames(AKAP8L) = c('Motif', 'Score')
AKAP8L$Motif = gsub("T", "U", AKAP8L$Motif)
AKAP8L = AKAP8L %>% arrange(Motif)

RBP = 'APOBEC3C'
APOBEC3C = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
APOBEC3C = APOBEC3C[, -1]
colnames(APOBEC3C) = c('Motif', 'Score')
APOBEC3C$Motif = gsub("T", "U", APOBEC3C$Motif)
APOBEC3C = APOBEC3C %>% arrange(Motif)

cell_line = 'HepG2'
RBP = 'RBM5'
RBM5 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
RBM5 = RBM5[, -1]
colnames(RBM5) = c('Motif', 'Score')
RBM5$Motif = gsub("T", "U", RBM5$Motif)
RBM5 = RBM5 %>% arrange(Motif)

RBP = 'IGF2BP3'
IGF2BP3 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
IGF2BP3 = IGF2BP3[, -1]
colnames(IGF2BP3) = c('Motif', 'Score')
IGF2BP3$Motif = gsub("T", "U", IGF2BP3$Motif)
IGF2BP3 = IGF2BP3 %>% arrange(Motif)

RBP = 'EIF3D'
EIF3D = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
EIF3D = EIF3D[, -1]
colnames(EIF3D) = c('Motif', 'Score')
EIF3D$Motif = gsub("T", "U", EIF3D$Motif)
EIF3D = EIF3D %>% arrange(Motif)

eCLIP_unKnown = data.frame(MOTIF = SAFB2$Motif,
                              LIN28B = NA,
                              TROVE2 = NA,
                              XRCC6 = XRCC6$Score,
                              SAFB2 = SAFB2$Score,
                              ZRANB2 = ZRANB2$Score,
                              AKAP8L = AKAP8L$Score,
                              APOBEC3C = APOBEC3C$Score,
                              RBM5 = RBM5$Score,
                              IGF2BP3 = IGF2BP3$Score,
                              EIF3D = EIF3D$Score)


for (RBP in c('LIN28B', 'TROVE2')) {
  K562 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/K562/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
  K562 = K562[, -1]
  colnames(K562) = c('Motif', 'Score')
  K562$Motif = gsub("T", "U", K562$Motif)
  K562 = K562 %>% arrange(Motif)
  
  HepG2 = read_csv(paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/HepG2/2_Scaled/', RBP, '_', K, 'mer.csv'), col_names = T,  show_col_types = F)
  HepG2 = HepG2[, -1]
  colnames(HepG2) = c('Motif', 'Score')
  HepG2$Motif = gsub("T", "U", HepG2$Motif)
  HepG2 = HepG2 %>% arrange(Motif)
  
  eCLIP_unKnown[, RBP] = min_max_norm((HepG2$Score + K562$Score)/2)
}

eCLIP_all = cbind(eCLIP_noStrucutures, eCLIP_Strucutures[, -1], eCLIP_unKnown[, -1])
################################################################################

# Make a summary table for RBPs:
################################################################################
RBPs = data.frame(RBP = c(colnames(eCLIP_noStrucutures)[-1],
                          colnames(eCLIP_Strucutures)[-1],
                          colnames(eCLIP_unKnown)[-1]),
                  Structure = c(rep('N', length(colnames(eCLIP_noStrucutures)[-1])),
                                rep('Y', length(colnames(eCLIP_Strucutures)[-1])),
                                rep('?', length(colnames(eCLIP_unKnown)[-1]))))
################################################################################

# Make table containing all necessary analysis:
################################################################################
eCLIP_Analysis = data.frame(RBP = RBPs$RBP,
                            Structure = RBPs$Structure,
                            RBNS_Motif = NA,
                            RBNS_Specificity = NA,
                            RBNS_Sensitivity = NA,
                            eCLIP_Motif = NA,
                            eCLIP_Specificity = NA,
                            eCLIP_Sensitivity = NA,
                            RBNS_Motif_eCLIP_Score = NA,
                            eCLIP_Motif_RBNS_Score = NA,
                            RBNS_Motif_eCLIP_Specificity = NA,
                            RBNS_Motif_eCLIP_Sensitivity = NA,
                            eCLIP_Motif_RBNS_Specificity = NA,
                            eCLIP_Motif_RBNS_Sensitivity = NA,
                            RBNS_v_eCLIP_Pearson = NA,
                            RBNS_v_eCLIP_Spearman = NA)

for (RBP in RBPs$RBP) {
  RBNS = RBNS_5mer[, c('Motif', RBP)]
  colnames(RBNS) = c('Motif', 'Score')
  RBNS$Motif = gsub("T", "U", RBNS$Motif)
  RBNS = RBNS %>% arrange(Motif)
  RBNS$Rank = rank(-RBNS$Score)
  
  eCLIP = eCLIP_all[, c('MOTIF', RBP)]
  colnames(eCLIP) = c('Motif', 'Score')
  eCLIP$Motif = gsub("T", "U", eCLIP$Motif)
  eCLIP = eCLIP %>% arrange(Motif)
  eCLIP$Rank = rank(-eCLIP$Score)
  
  IS = data.frame(Motif = eCLIP$Motif,
                  Inherent_Specificity = SI_Calculation(eCLIP$Score))
  IS$Motif = gsub("T", "U", IS$Motif)
  IS = IS %>% arrange(-Inherent_Specificity)
  
  MS = data.frame(Motif = eCLIP$Motif,
                  Mutational_Sensitivity = NA)
  
  for (Motif in MS$Motif) {
    Motif_Vars = gsub("T", "U", Motif_Variants(Motif, NTs))
    MSs = c()
    for (idx in c(1:length(Motif_Vars))) {
      Motif_Var = Motif_Vars[idx]
      if (Motif_Var != Motif) {
        MSs[idx] = 1 - eCLIP[which(eCLIP$Motif == Motif_Var), ]$Score
      }
    }
    MS$Mutational_Sensitivity[which(MS$Motif == Motif)] = mean(MSs, na.rm = T)
  }
  
  eCLIP_Analysis[which(eCLIP_Analysis$RBP == RBP), 'RBNS_Motif'] = RBNS$Motif[which(RBNS$Rank == 1)]
  eCLIP_Analysis[which(eCLIP_Analysis$RBP == RBP), 'RBNS_Specificity'] = max(RBNS$Score)/median(RBNS$Score)
  eCLIP_Analysis[which(eCLIP_Analysis$RBP == RBP), 'RBNS_Sensitivity'] = RBNS_Stats$MS[which(RBNS_Stats$RBP == RBP)]
  
  eCLIP_Analysis[which(eCLIP_Analysis$RBP == RBP), 'eCLIP_Motif'] = eCLIP$Motif[which(eCLIP$Rank == 1)]
  eCLIP_Analysis[which(eCLIP_Analysis$RBP == RBP), 'eCLIP_Specificity'] = IS$Inherent_Specificity[which(IS$Motif == eCLIP$Motif[which(eCLIP$Rank == 1)])]
  eCLIP_Analysis[which(eCLIP_Analysis$RBP == RBP), 'eCLIP_Sensitivity'] = MS$Mutational_Sensitivity[which(MS$Motif == eCLIP$Motif[which(eCLIP$Rank == 1)])]

  eCLIP_Analysis[which(eCLIP_Analysis$RBP == RBP), 'RBNS_Motif_eCLIP_Score'] = eCLIP$Score[which(eCLIP$Motif == RBNS$Motif[which(RBNS$Rank == 1)])]
  eCLIP_Analysis[which(eCLIP_Analysis$RBP == RBP), 'eCLIP_Motif_RBNS_Score'] = RBNS$Score[which(RBNS$Motif == eCLIP$Motif[which(eCLIP$Rank == 1)])]

  eCLIP_Analysis[which(eCLIP_Analysis$RBP == RBP), 'RBNS_Motif_eCLIP_Specificity'] = IS$Inherent_Specificity[which(IS$Motif == RBNS$Motif[which(RBNS$Rank == 1)])]
  eCLIP_Analysis[which(eCLIP_Analysis$RBP == RBP), 'RBNS_Motif_eCLIP_Sensitivity'] = MS$Mutational_Sensitivity[which(MS$Motif == RBNS$Motif[which(RBNS$Rank == 1)])]
  eCLIP_Analysis[which(eCLIP_Analysis$RBP == RBP), 'eCLIP_Motif_RBNS_Specificity'] = RBNS$Score[which(RBNS$Motif == eCLIP$Motif[which(eCLIP$Rank == 1)])]/median(RBNS$Score)
  eCLIP_Analysis[which(eCLIP_Analysis$RBP == RBP), 'eCLIP_Motif_RBNS_Sensitivity'] = mean(MSs)

  eCLIP_Analysis[which(eCLIP_Analysis$RBP == RBP), 'RBNS_v_eCLIP_Pearson'] = cor(RBNS$Score, eCLIP$Score, method = c('pearson'))
  eCLIP_Analysis[which(eCLIP_Analysis$RBP == RBP), 'RBNS_v_eCLIP_Spearman'] = cor(RBNS$Score, eCLIP$Score, method = c('spearman'))
}


# write_csv(eCLIP_Analysis, paste0(baseDir, 'Data/', DATE, '_', K, '_mer_RBNS_eCLIP_Analysis.csv'), quote = 'none', col_names = T)
################################################################################

## hnRNPC: eCLIP vs RBNS
################################################################################
# Read RBNS file:
hnRNPC_RBNS = RBNS_5mer[, c('Motif', 'HNRNPC')]
colnames(hnRNPC_RBNS) = c('Motif', 'Score')
hnRNPC_RBNS = hnRNPC_RBNS %>% arrange(Motif)

hnRNPC_eCLIP = eCLIP_all[, c('MOTIF', 'HNRNPC')]
colnames(hnRNPC_eCLIP) = c('Motif', 'Score')
hnRNPC_eCLIP$Motif = gsub("T", "U", hnRNPC_eCLIP$Motif)
hnRNPC_eCLIP = hnRNPC_eCLIP %>% arrange(Motif)

hnRNPC_plot = data.frame(MOTIF = hnRNPC_eCLIP$Motif,
                         ECLIP = hnRNPC_eCLIP$Score,
                         RBNS = hnRNPC_RBNS$Score)

temp = hnRNPC_plot %>% mutate(selection = pmax(ECLIP, RBNS))
temp = temp %>% mutate(index = row_number())
temp = temp %>% arrange(desc(selection))

Corrs = c()
Corrs[1] = cor(hnRNPC_plot$ECLIP, hnRNPC_plot$RBNS, method = c('pearson'))

for (idx in c(1:9)) {
  Corrs[idx+1] = cor(hnRNPC_plot[-temp$index[1:idx], 'ECLIP'], hnRNPC_plot[-temp$index[1:idx], 'RBNS'], method = c('pearson'))
}

ggplot() +
  geom_point(data = hnRNPC_plot, aes(x = (RBNS), y = (ECLIP)), alpha = 0.7, shape = 19, size = 2) +
  geom_smooth(data = hnRNPC_plot, aes(x = (RBNS), y = (ECLIP)), method = "lm", linetype = 'solid', color = 'salmon', se = F, linewidth = 0.5, alpha = 0.1) +
  geom_smooth(data = hnRNPC_plot[-temp$index[1:1], ], aes(x = (RBNS), y = (ECLIP)), method = "lm", linetype = 'solid', color = 'skyblue', se = F, linewidth = 0.5, alpha = 0.1) +
  theme_minimal() +  
  labs(x = 'RBNS', 
       y = 'eCLIP', 
       title = paste0('hnRNPC RBNS vs eCLIP')) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  annotate("text", x = Inf, y = 0.25, label = sprintf("%.2f", Corrs[1]), hjust = 1.1, vjust = 0, size = 5, color = 'salmon') + theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) +
  annotate("text", x = Inf, y = 0.2, label = sprintf("%.2f", Corrs[7]), hjust = 1.1, vjust = 0, size = 5, color = 'skyblue') + theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))

################################################################################

## PCBP2: eCLIP vs RBNS
################################################################################
PCBP2_RBNS = RBNS_5mer[, c('Motif', 'PCBP2')]
colnames(PCBP2_RBNS) = c('Motif', 'Score')
PCBP2_RBNS = PCBP2_RBNS %>% arrange(Motif)

PCBP2_eCLIP = eCLIP_all[, c('MOTIF', 'PCBP2')]
colnames(PCBP2_eCLIP) = c('Motif', 'Score')
PCBP2_eCLIP$Motif = gsub("T", "U", PCBP2_eCLIP$Motif)
PCBP2_eCLIP = PCBP2_eCLIP %>% arrange(Motif)

PCBP2_plot = data.frame(MOTIF = PCBP2_eCLIP$Motif,
                        ECLIP = PCBP2_eCLIP$Score,
                        RBNS = PCBP2_RBNS$Score)

temp = PCBP2_plot %>% mutate(selection = pmax(ECLIP, RBNS))
temp = temp %>% mutate(index = row_number())
temp = temp %>% arrange(desc(selection))

Corrs = c()
Corrs[1] = cor(PCBP2_plot$ECLIP, PCBP2_plot$RBNS, method = c('pearson'))

for (idx in c(1:9)) {
  Corrs[idx+1] = cor(PCBP2_plot[-temp$index[1:idx], 'ECLIP'], PCBP2_plot[-temp$index[1:idx], 'RBNS'], method = c('pearson'))
}

ggplot() +
  geom_point(data = PCBP2_plot, aes(x = (RBNS), y = (ECLIP)), alpha = 0.7, shape = 19, size = 2) +
  geom_smooth(data = PCBP2_plot, aes(x = (RBNS), y = (ECLIP)), method = "lm", linetype = 'solid', color = 'salmon', se = F, linewidth = 0.5, alpha = 0.1) +
  geom_smooth(data = PCBP2_plot[-temp$index[1:9], ], aes(x = (RBNS), y = (ECLIP)), method = "lm", linetype = 'solid', color = 'skyblue', se = F, linewidth = 0.5, alpha = 0.1) +
  theme_minimal() +  
  labs(x = 'RBNS', 
       y = 'eCLIP', 
       title = paste0('PCBP2 RBNS vs eCLIP')) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  annotate("text", x = Inf, y = 0.25, label = sprintf("%.2f", Corrs[4]), hjust = 1.1, vjust = 0, size = 5, color = 'salmon') + theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) +
  annotate("text", x = Inf, y = 0.2, label = sprintf("%.2f", Corrs[10]), hjust = 1.1, vjust = 0, size = 5, color = 'skyblue') + theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))


################################################################################

## RBFOX2: eCLIP vs RBNS
################################################################################
RBFOX2_RBNS = RBNS_5mer[, c('Motif', 'RBFOX2')]
colnames(RBFOX2_RBNS) = c('Motif', 'Score')
RBFOX2_RBNS = RBFOX2_RBNS %>% arrange(Motif)

RBFOX2_eCLIP = eCLIP_all[, c('MOTIF', 'RBFOX2')]
colnames(RBFOX2_eCLIP) = c('Motif', 'Score')
RBFOX2_eCLIP$Motif = gsub("T", "U", RBFOX2_eCLIP$Motif)
RBFOX2_eCLIP = RBFOX2_eCLIP %>% arrange(Motif)

RBFOX2_plot = data.frame(MOTIF = RBFOX2_eCLIP$Motif,
                         ECLIP = RBFOX2_eCLIP$Score,
                         RBNS = RBFOX2_RBNS$Score)

temp = RBFOX2_plot %>% mutate(selection = pmax(ECLIP, RBNS))
temp = temp %>% mutate(index = row_number())
temp = temp %>% arrange(desc(selection))

Corrs = c()
Corrs[1] = cor(RBFOX2_plot$ECLIP, RBFOX2_plot$RBNS, method = c('pearson'))

for (idx in c(1:9)) {
  Corrs[idx+1] = cor(RBFOX2_plot[-temp$index[1:idx], 'ECLIP'], RBFOX2_plot[-temp$index[1:idx], 'RBNS'], method = c('pearson'))
}

ggplot() +
  geom_point(data = RBFOX2_plot, aes(x = (RBNS), y = (ECLIP)), alpha = 0.7, shape = 19, size = 2) +
  geom_smooth(data = RBFOX2_plot, aes(x = (RBNS), y = (ECLIP)), method = "lm", linetype = 'solid', color = 'salmon', se = F, linewidth = 0.5, alpha = 1) +
  geom_smooth(data = RBFOX2_plot[-temp$index[1:2], ], aes(x = (RBNS), y = (ECLIP)), method = "lm", linetype = 'solid', color = 'skyblue', se = F, linewidth = 0.5, alpha = 1) +
  theme_minimal() +  
  labs(x = 'RBNS', 
       y = 'eCLIP', 
       title = paste0('RBFOX2 RBNS vs eCLIP')) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  annotate("text", x = Inf, y = 0.2, label = sprintf("%.2f", Corrs[1]), hjust = 1.1, vjust = 0, size = 5, color = 'salmon') + theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) +
  annotate("text", x = Inf, y = 0.25, label = sprintf("%.2f", Corrs[6]), hjust = 1.1, vjust = 0, size = 5, color = 'skyblue') + theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))

################################################################################

## EIF4G2: eCLIP vs RBNS
################################################################################
EIF4G2_RBNS = RBNS_5mer[, c('Motif', 'EIF4G2')]
colnames(EIF4G2_RBNS) = c('Motif', 'Score')
EIF4G2_RBNS = EIF4G2_RBNS %>% arrange(Motif)

EIF4G2_eCLIP = eCLIP_all[, c('MOTIF', 'EIF4G2')]
colnames(EIF4G2_eCLIP) = c('Motif', 'Score')
EIF4G2_eCLIP$Motif = gsub("T", "U", EIF4G2_eCLIP$Motif)
EIF4G2_eCLIP = EIF4G2_eCLIP %>% arrange(Motif)

EIF4G2_plot = data.frame(MOTIF = EIF4G2_eCLIP$Motif,
                         ECLIP = EIF4G2_eCLIP$Score,
                         RBNS = EIF4G2_RBNS$Score)

temp = EIF4G2_plot %>% mutate(selection = pmax(ECLIP, RBNS))
temp = temp %>% mutate(index = row_number())
temp = temp %>% arrange(desc(selection))

Corrs = c()
Corrs[1] = cor(EIF4G2_plot$ECLIP, EIF4G2_plot$RBNS, method = c('pearson'))

for (idx in c(1:9)) {
  Corrs[idx+1] = cor(EIF4G2_plot[-temp$index[1:idx], 'ECLIP'], EIF4G2_plot[-temp$index[1:idx], 'RBNS'], method = c('pearson'))
}

ggplot() +
  geom_point(data = EIF4G2_plot, aes(x = (RBNS), y = (ECLIP)), alpha = 0.7, shape = 19, size = 2) +
  geom_smooth(data = EIF4G2_plot, aes(x = (RBNS), y = (ECLIP)), method = "lm", linetype = 'solid', color = 'salmon', se = F, linewidth = 0.5, alpha = 1) +
  geom_smooth(data = EIF4G2_plot[-temp$index[1:2], ], aes(x = (RBNS), y = (ECLIP)), method = "lm", linetype = 'solid', color = 'skyblue', se = F, linewidth = 0.5, alpha = 1) +
  theme_minimal() +  
  labs(x = 'RBNS', 
       y = 'eCLIP', 
       title = paste0('EIF4G2 RBNS vs eCLIP')) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  annotate("text", x = Inf, y = 1, label = sprintf("%.2f", Corrs[10]), hjust = 1.1, vjust = 0, size = 5, color = 'salmon') + theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) +
  annotate("text", x = Inf, y = 0.9, label = sprintf("%.2f", Corrs[1]), hjust = 1.1, vjust = 0, size = 5, color = 'skyblue') + theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))

################################################################################



## Create dataframe for plot:
################################################################################
temp_plotData = data.frame(RBP = eCLIP_Analysis$RBP,
                           IS = eCLIP_Analysis$RBNS_Specificity,
                           IT = eCLIP_Analysis$RBNS_Sensitivity,
                           AS = eCLIP_Analysis$eCLIP_Specificity,
                           AT = eCLIP_Analysis$eCLIP_Sensitivity,
                           Structure = eCLIP_Analysis$Structure,
                           Pearson = eCLIP_Analysis$RBNS_v_eCLIP_Pearson,
                           Spearman = eCLIP_Analysis$RBNS_v_eCLIP_Spearman)
################################################################################

## Specificity vs Sensitivity for Apparent and Inherent
################################################################################
Corr_Y = cor(temp_plotData$AS, temp_plotData$AT, method = c('spearman'))
Corr_N = cor(temp_plotData$IS, temp_plotData$IT, method = c('spearman'))

ggplot() +
  geom_point(data = temp_plotData, aes(x = (AS), y = AT), fill = 'cornflowerblue', alpha = 1, shape = 21, size = 4) +
  geom_point(data = temp_plotData, aes(x = (IS), y = IT), fill = 'darkgoldenrod1', alpha = 1, shape = 21, size = 4) +
  geom_smooth(data = temp_plotData, aes(x = (AS), y = AT), method = 'lm', linetype = 'solid', color = 'cornflowerblue', se = F, linewidth = 0.5, alpha = 0.1) +
  geom_smooth(data = temp_plotData, aes(x = (IS), y = IT), method = 'lm', linetype = 'solid', color = 'darkgoldenrod1', se = F, linewidth = 0.5, alpha = 0.1) +
  theme_minimal() +
  labs(x = 'Specificity',
       y = 'Sensitivity',
       title = paste0('Specificity vs Sensitivity Correlation')) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(trans = 'log2', limits = c(1, 64)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14)) +
  annotate("text", x = Inf, y = 0.2, label = sprintf("%.2f", Corr_Y), hjust = 1.1, vjust = 0, size = 5, color = 'cornflowerblue') +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) +
  annotate("text", x = Inf, y = 0.15, label = sprintf("%.2f", Corr_N), hjust = 1.1, vjust = 0, size = 5, color = 'darkgoldenrod1') +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))

################################################################################

## Apparent vs Inherent Specificity
################################################################################
Corr_Y = cor((temp_plotData %>% filter(Structure == 'Y'))$AS, (temp_plotData %>% filter(Structure == 'Y'))$IS, method = c('spearman'))
Corr_N = cor((temp_plotData %>% filter(Structure == 'N'))$AS, (temp_plotData %>% filter(Structure == 'N'))$IS, method = c('spearman'))

ggplot() +
  geom_point(data = temp_plotData %>% filter(Structure == 'Y'), aes(x = IS, y = AS), fill = 'antiquewhite2', alpha = 1, shape = 21, size = 4) +
  geom_smooth(data = temp_plotData %>% filter(Structure == 'Y'), aes(x = IS, y = AS), method = "lm", linetype = 'solid', color = 'antiquewhite2', se = F, linewidth = 0.5, alpha = 0.1) +
  geom_point(data = temp_plotData %>% filter(Structure == 'N'), aes(x = IS, y = AS), fill = 'darkseagreen4', alpha = 1, shape = 21, size = 4) +
  geom_smooth(data = temp_plotData %>% filter(Structure == 'N'), aes(x = IS, y = AS), method = "lm", linetype = 'solid', color = 'darkseagreen4', se = F, linewidth = 0.5, alpha = 0.1) +
  theme_minimal() +  
  labs(x = 'Inherent Specificity', 
       y = 'Apparent Specificity', 
       title = paste0('Apparent vs Inherent Specificity')) +
  scale_x_continuous(trans = 'log2', limits = c(1, 64)) +
  scale_y_continuous(trans = 'log2', limits = c(1, 64)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  annotate("text", x = Inf, y = 2, label = sprintf("%.2f", Corr_Y), hjust = 1.1, vjust = 0, size = 5, color = 'antiquewhite2') +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) +
  annotate("text", x = Inf, y = 1.5, label = sprintf("%.2f", Corr_N), hjust = 1.1, vjust = 0, size = 5, color = 'darkseagreen4') +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
################################################################################

## Apparent vs Inherent Sensitivity
################################################################################
Corr_Y = cor((temp_plotData %>% filter(Structure == 'Y'))$AT, (temp_plotData %>% filter(Structure == 'Y'))$IT, method = c('pearson'))
Corr_N = cor((temp_plotData %>% filter(Structure == 'N'))$AT, (temp_plotData %>% filter(Structure == 'N'))$IT, method = c('pearson'))

ggplot() +
  geom_point(data = temp_plotData %>% filter(Structure == 'Y'), aes(x = IT, y = AT), fill = 'antiquewhite2', alpha = 1, shape = 21, size = 4) +
  geom_smooth(data = temp_plotData %>% filter(Structure == 'Y'), aes(x = IT, y = AT), method = "lm", linetype = 'solid', color = 'antiquewhite2', se = F, linewidth = 0.5, alpha = 0.1) +
  geom_point(data = temp_plotData %>% filter(Structure == 'N'), aes(x = IT, y = AT), fill = 'darkseagreen4', alpha = 1, shape = 21, size = 4) +
  geom_smooth(data = temp_plotData %>% filter(Structure == 'N'), aes(x = IT, y = AT), method = "lm", linetype = 'solid', color = 'darkseagreen4', se = F, linewidth = 0.5, alpha = 0.1) +
  theme_minimal() +  
  labs(x = 'Inherent Mutational Sensitivity', 
       y = 'Apparent Mutational Sensitivity', 
       title = paste0('Apparent vs Inherent Sensitivity')) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  annotate("text", x = Inf, y = 0.1, label = sprintf("%.2f", Corr_Y), hjust = 1.1, vjust = 0, size = 5, color = 'antiquewhite2') +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) +
  annotate("text", x = Inf, y = 0.2, label = sprintf("%.2f", Corr_N), hjust = 1.1, vjust = 0, size = 5, color = 'darkseagreen4') +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
################################################################################

## Inherent Specificity Boxplot by Pearson Correlation ##
################################################################################
temp_plotData = data.frame(RBP = eCLIP_Analysis$RBP,
                           IS = eCLIP_Analysis$RBNS_Specificity,
                           IT = eCLIP_Analysis$RBNS_Sensitivity,
                           AS = eCLIP_Analysis$eCLIP_Specificity,
                           AT = eCLIP_Analysis$eCLIP_Sensitivity,
                           Structure = eCLIP_Analysis$Structure,
                           Pearson = eCLIP_Analysis$RBNS_v_eCLIP_Pearson,
                           Spearman = eCLIP_Analysis$RBNS_v_eCLIP_Spearman)

temp_plotData = temp_plotData %>% filter(Structure != '?')

ggplot(temp_plotData, aes(x = Structure, y = (Spearman))) +
  geom_boxplot(notch = F, outlier.shape = NA) +
  geom_jitter(aes(shape = Structure), width = 0.1, size = 3, alpha = 1) +
  theme_minimal() +
  labs(title = "Spearman Correlation by Structure",
       x = "Structure",
       y = "Spearman") +
  scale_y_continuous(limits = c(-0.5, 1), breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))+
  geom_signif(comparisons = list(c("N", "Y")),
              test = c('ks.test'),
              map_signif_level = F,
              textsize = 5, y_position = 0.8, tip_length = 0.01)

ggplot(temp_plotData, aes(x = Structure, y = (Pearson))) +
  geom_boxplot(notch = F, outlier.shape = NA) +
  geom_jitter(aes(shape = Structure), width = 0.1, size = 3, alpha = 1) +
  theme_minimal() +
  labs(title = "Pearson Correlation by Structure",
       x = "Structure",
       y = "Pearson") +
  scale_y_continuous(limits = c(-0.5, 1), breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))+
  geom_signif(comparisons = list(c("N", "Y")),
              test = c('ks.test'),
              map_signif_level = F,
              textsize = 5, y_position = 0.8, tip_length = 0.01)

################################################################################
