################################################################################
## RBNS Metric calculation
## Written by Soon Yi
## Created: 2023-12-09
## Last Edited: 2024-02-16
################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(tibble)

# Set some basic parameters and custom functions:
################################################################################
baseDir = '~/Desktop/Genomics/Specificity/Data_Final/'

NTs = c("A", "C", "G", "U")

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

# Read RBNS file:
################################################################################
RBNS_4mer = read_csv(paste0(baseDir, 'RBNS/RBNS_normalized_4mer.csv'), col_names = T, show_col_types = F)
RBNS_5mer = read_csv(paste0(baseDir, 'RBNS/RBNS_normalized_5mer.csv'), col_names = T, show_col_types = F)
RBNS_6mer = read_csv(paste0(baseDir, 'RBNS/RBNS_normalized_6mer.csv'), col_names = T, show_col_types = F)
RBNS_7mer = read_csv(paste0(baseDir, 'RBNS/RBNS_normalized_7mer.csv'), col_names = T, show_col_types = F)

RBPs_4mer = colnames(RBNS_4mer)[2:ncol(RBNS_4mer)]
RBPs_5mer = colnames(RBNS_5mer)[2:ncol(RBNS_5mer)]
RBPs_6mer = colnames(RBNS_6mer)[2:ncol(RBNS_6mer)]
RBPs_7mer = colnames(RBNS_7mer)[2:ncol(RBNS_7mer)]
################################################################################

# Calculate SI and MS for 4-mer:
################################################################################
RBPs = RBPs_4mer
RBNS = RBNS_4mer

Metrics_Kmer = as.data.frame(matrix(NA, nrow = length(RBPs), ncol = 3))
colnames(Metrics_Kmer) = c('RBP', 'SI', 'MS')
Metrics_Kmer$RBP = RBPs

for (RBP in RBPs) {
  Metrics_Kmer[which(Metrics_Kmer$RBP == RBP), 'SI'] = max(unlist(RBNS[, RBP])) / median(unlist(RBNS[, RBP]))
  
  Motif_Top = RBNS$Motif[which(RBNS[, RBP] == max(RBNS[, RBP]))]
  Motif_Vars = Motif_Variants(Motif_Top, NTs)
  MSs = c()
  for (idx in c(1:length(Motif_Vars))) {
    Motif_Var = Motif_Vars[idx]
    MSs[idx] = 1 - RBNS[which(RBNS$Motif == Motif_Var), ][[RBP]][1]
  }
  Metrics_Kmer[which(Metrics_Kmer$RBP == RBP), 'MS'] = mean(MSs)

}

Metrics_Kmer = Metrics_Kmer %>% arrange(RBP)
Metrics_4mer = Metrics_Kmer
################################################################################

# Calculate SI and MS for 5-mer:
################################################################################
RBPs = RBPs_5mer
RBNS = RBNS_5mer

Metrics_Kmer = as.data.frame(matrix(NA, nrow = length(RBPs), ncol = 3))
colnames(Metrics_Kmer) = c('RBP', 'SI', 'MS')
Metrics_Kmer$RBP = RBPs

for (RBP in RBPs) {
  Metrics_Kmer[which(Metrics_Kmer$RBP == RBP), 'SI'] = max(unlist(RBNS[, RBP])) / median(unlist(RBNS[, RBP]))
  
  Motif_Top = RBNS$Motif[which(RBNS[, RBP] == max(RBNS[, RBP]))]
  Motif_Vars = Motif_Variants(Motif_Top, NTs)
  MSs = c()
  for (idx in c(1:length(Motif_Vars))) {
    Motif_Var = Motif_Vars[idx]
    MSs[idx] = 1 - RBNS[which(RBNS$Motif == Motif_Var), ][[RBP]][1]
  }
  Metrics_Kmer[which(Metrics_Kmer$RBP == RBP), 'MS'] = mean(MSs)
  
}

Metrics_Kmer = Metrics_Kmer %>% arrange(RBP)
Metrics_5mer = Metrics_Kmer
################################################################################

# Calculate SI and MS for 6-mer:
################################################################################
RBPs = RBPs_6mer
RBNS = RBNS_6mer

Metrics_Kmer = as.data.frame(matrix(NA, nrow = length(RBPs), ncol = 3))
colnames(Metrics_Kmer) = c('RBP', 'SI', 'MS')
Metrics_Kmer$RBP = RBPs

for (RBP in RBPs) {
  Metrics_Kmer[which(Metrics_Kmer$RBP == RBP), 'SI'] = max(unlist(RBNS[, RBP])) / median(unlist(RBNS[, RBP]))
  
  Motif_Top = RBNS$Motif[which(RBNS[, RBP] == max(RBNS[, RBP]))]
  Motif_Vars = Motif_Variants(Motif_Top, NTs)
  MSs = c()
  for (idx in c(1:length(Motif_Vars))) {
    Motif_Var = Motif_Vars[idx]
    MSs[idx] = 1 - RBNS[which(RBNS$Motif == Motif_Var), ][[RBP]][1]
  }
  Metrics_Kmer[which(Metrics_Kmer$RBP == RBP), 'MS'] = mean(MSs)
  
}

Metrics_Kmer = Metrics_Kmer %>% arrange(RBP)
Metrics_6mer = Metrics_Kmer
################################################################################

# Calculate SI and MS for 7-mer:
################################################################################
RBPs = RBPs_7mer
RBNS = RBNS_7mer

Metrics_Kmer = as.data.frame(matrix(NA, nrow = length(RBPs), ncol = 3))
colnames(Metrics_Kmer) = c('RBP', 'SI', 'MS')
Metrics_Kmer$RBP = RBPs

for (RBP in RBPs) {
  Metrics_Kmer[which(Metrics_Kmer$RBP == RBP), 'SI'] = max(unlist(RBNS[, RBP])) / median(unlist(RBNS[, RBP]))
  
  Motif_Top = RBNS$Motif[which(RBNS[, RBP] == max(RBNS[, RBP]))]
  Motif_Vars = Motif_Variants(Motif_Top, NTs)
  MSs = c()
  for (idx in c(1:length(Motif_Vars))) {
    Motif_Var = Motif_Vars[idx]
    MSs[idx] = 1 - RBNS[which(RBNS$Motif == Motif_Var), ][[RBP]][1]
  }
  Metrics_Kmer[which(Metrics_Kmer$RBP == RBP), 'MS'] = mean(MSs)
  
}

Metrics_Kmer = Metrics_Kmer %>% arrange(RBP)
Metrics_7mer = Metrics_Kmer
################################################################################

# Read RBNS file:
################################################################################
write.csv(Metrics_4mer, paste0(baseDir, 'RBNS/4mer_metrics.csv'), quote = F, row.names = F)
write.csv(Metrics_5mer, paste0(baseDir, 'RBNS/5mer_metrics.csv'), quote = F, row.names = F)
write.csv(Metrics_6mer, paste0(baseDir, 'RBNS/6mer_metrics.csv'), quote = F, row.names = F)
write.csv(Metrics_7mer, paste0(baseDir, 'RBNS/7mer_metrics.csv'), quote = F, row.names = F)
################################################################################
