################################################################################
## Generate Distribution from RBNS
## Written by Soon Yi
## Created: 2023-11-07
## Last Edited: 2024-02-11
## Figure 1
################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

## Set some basic parameters:
################################################################################
baseDir = '~/Desktop/Genomics/Specificity/Data_Final/'

return_MS = function(data) {
  colnames(data) = c('Motif', 'Score')
  MS = data.frame(Nucleotide = c('A', 'C', 'G', 'U'),
                  Pos1 = NA,
                  Pos2 = NA,
                  Pos3 = NA,
                  Pos4 = NA,
                  Pos5 = NA)
  
  Motif_Top = data$Motif[which(data[, 'Score'] == max(data[, 'Score']))]
  
  for (pos in c(1, 2, 3, 4, 5)) {
    for (nt in c('A', 'C', 'G', 'U')) {
      Motif_Var = Motif_Top
      substring(Motif_Var, pos, pos) = nt
      MS[which(MS$Nucleotide == nt), paste0('Pos', pos)] = 1 - data[which(data$Motif == Motif_Var), 'Score']
    }
  }
  
  return(MS)
}

plot_MS = function(data) {
  MS_long = reshape2::melt(data, id.vars = 'Nucleotide', measure.vars = c('Pos1', 'Pos2', 'Pos3', 'Pos4', 'Pos5'))
  colnames(MS_long) = c('Nucleotide', 'Position', 'Value')
  MS_long$Position = gsub("Pos", "", MS_long$Position)
  MS_long$Position = as.numeric(as.character(MS_long$Position))
  
  plot = ggplot(MS_long, aes(x = Position, y = Value, fill = Nucleotide)) +
    geom_point(pch = 21, size = 10, stroke = 1) +
    scale_fill_manual(values = c("G" = "#F5C714", "A" = "#70BF52", "C" = "#3D94D1", "U" = "#E0546C")) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = 1:5) +
    scale_y_continuous(limits = c(1.05, -0.05), breaks = seq(0, 1, by = 0.2), trans = 'reverse') +
    coord_cartesian(xlim = c(0.9, 5.1)) + 
    labs(x = "Position") + 
    theme_minimal() +
    theme_bw() + 
    theme(axis.text = element_text(size=14), 
          axis.title = element_text(size=14, face = 'bold'), 
          legend.text = element_text(size=14))
  
  return(plot)
  
}

K = 5
RBNS = read_csv(paste0(baseDir, 'RBNS/RBNS_normalized_', K, 'mer.csv'), col_names = T, show_col_types = F)
colnames(RBNS) = c('Motif', colnames(RBNS)[2:ncol(RBNS)])
RBPs = colnames(RBNS)[2:ncol(RBNS)]
################################################################################

# Read RBNS file:
hnRNPC = 'hnRNPC'
hnRNPC = RBNS[, c('Motif', 'HNRNPC')]
colnames(hnRNPC) = c('motif', 'enrichment')

EIF4G2 = 'EIF4G2'
EIF4G2 = RBNS[, c('Motif', 'EIF4G2')]
colnames(EIF4G2) = c('motif', 'enrichment')


ggplot() +
  geom_histogram(data = hnRNPC, aes(x = enrichment), binwidth = 0.001, fill = "black", alpha = 1.0) +
  labs(title = "hnRNPC", x = "RBNS Score", y = "Frequency") +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 25)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot() +
  geom_histogram(data = EIF4G2, aes(x = enrichment), binwidth = 0.001, fill = "black", alpha = 1.0) +
  labs(title = "EIF4G2", x = "RBNS Score", y = "Frequency") +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 25)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


hnRNPC_MS = plot_MS(return_MS(hnRNPC))

EIF4G2_MS = plot_MS(return_MS(EIF4G2))

