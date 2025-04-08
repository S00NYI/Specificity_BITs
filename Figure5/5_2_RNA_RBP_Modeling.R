################################################################################
## RNA-RBP Modeling 
## Written by Soon Yi
## Created: 2023-12-11
## Last Edited: 2023-12-11
## Figure 5
################################################################################

## Load libraries
################################################################################
library(readr)
library(data.table)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(seqLogo)
library(Biostrings)
################################################################################

## Set basic parameters:
################################################################################
baseDir = './Specificity/'

featureScale = function(X, MAX, MIN) {
  MIN + ((X - min(X, na.rm = TRUE))*(MAX - MIN)) / (max(X, na.rm = TRUE) - min(X, na.rm = TRUE))
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

return_probMatrix = function(data, N) {
  nts = c('A', 'C', 'G', 'U')
  
  colnames(data) = c('Motif', 'Score')
  data = data %>% arrange(-Score)
  data_pseudoCounts = data[c(1:N), ]
  data_pseudoCounts$Score = round(data_pseudoCounts$Score*100, 0)
  
  data_countMatrix = matrix(0, nrow = 4, ncol = 5)
  rownames(data_countMatrix) = nts
  
  for(seqID in 1:nrow(data_pseudoCounts)) {
    seq = as.character(data_pseudoCounts$Motif[seqID])
    for(pos in 1:5) {
      count2Add = as.integer(data_pseudoCounts[seqID, 2]) # Assuming the counts are in the first column
      data_countMatrix[substr(seq, pos, pos), pos] = data_countMatrix[substr(seq, pos, pos), pos] + count2Add
    }
  }
  data_probMatrix = sweep(data_countMatrix, 2, colSums(data_countMatrix), FUN = "/")
  
  return(data_probMatrix)
}


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

return_IS = function(x) {
  x / (median(x, na.rm = FALSE))
}
################################################################################

## Set up query RNA:
################################################################################
# K:              K-mer
# query:          RNA sequence to test.
# L:              length of the query sequence.
# W:              number of possible K-mers in the query sequence.
# query_W:        all possible K-mers within the query sequence.
# unique_query_W: all unique K-mers within the query sequence.

K = 5
query = 'GUUUUUGUGGUGCAGAAGAUAGAGAGAGCGGGUUGCGAGUCUAGAAUGCGGGUAAGGAGUGCAGAGAGAUAGAGUUUAGUAUUAUGUUUUUUUUUUGUAU'
L = nchar(query)
W = L - K + 1

query_W = sapply(1:W, function(i) substr(query, i, i + K - 1))

unique_query_W = data.frame(motif = unique(query_W),
                            count = NA)

for (idx in c(1:nrow(unique_query_W))) {
  unique_query_W$count[idx] = sum(query_W == unique_query_W$motif[idx])
}
################################################################################

## Set up model RBP(s):
################################################################################
# Read in the model RBPs:
Model_HH = fread(paste0(baseDir, 'Data_Final/Model_RBP/Model_HH.csv'))
Model_HL = fread(paste0(baseDir, 'Data_Final/Model_RBP/Model_HL.csv'))
Model_LH = fread(paste0(baseDir, 'Data_Final/Model_RBP/Model_LH.csv'))
Model_LL = fread(paste0(baseDir, 'Data_Final/Model_RBP/Model_LL.csv'))

Model_HH = Model_HH %>% arrange(Motif)
Model_HL = Model_HL %>% arrange(Motif)
Model_LH = Model_LH %>% arrange(Motif)
Model_LL = Model_LL %>% arrange(Motif)
################################################################################

## Histograms for model RBP comparison
################################################################################
y_lim = 200
y_break = 50

ggplot() +
  geom_histogram(data = Model_HH, aes(x = featureScale(HH, MIN = 0, MAX = 1)), binwidth = 0.0025, fill = "black", alpha = 1.0) +
  labs(title = "RBP HH: High Inherent Specificity, High Mutational Sensitivity", x = "Normalized Affinity", y = "Frequency") +
  scale_y_continuous(limits = c(0, y_lim), breaks = seq(0, y_lim, by = y_break)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot() +
  geom_histogram(data = Model_HL, aes(x = featureScale(HL, MIN = 0, MAX = 1)), binwidth = 0.0025, fill = "black", alpha = 1.0) +
  labs(title = "RBP HL: High Inherent Specificity, Low Mutational Sensitivity", x = "Normalized Affinity", y = "Frequency") +
  scale_y_continuous(limits = c(0, y_lim), breaks = seq(0, y_lim, by = y_break)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

y_lim = 100
y_break = 25

ggplot() +
  geom_histogram(data = Model_LH, aes(x = featureScale(LH, MIN = 0, MAX = 1)), binwidth = 0.0025, fill = "black", alpha = 1.0) +
  labs(title = "RBP LH: Low Inherent Specificity, High Mutational Sensitivity", x = "Normalized Affinity", y = "Frequency") +
  scale_y_continuous(limits = c(0, y_lim), breaks = seq(0, y_lim, by = y_break)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot() +
  geom_histogram(data = Model_LL, aes(x = featureScale(LL, MIN = 0, MAX = 1)), binwidth = 0.0025, fill = "black", alpha = 1.0) +
  labs(title = "RBP LL: Low Inherent Specificity, Low Mutational Sensitivity", x = "Normalized Affinity", y = "Frequency") +
  scale_y_continuous(limits = c(0, y_lim), breaks = seq(0, y_lim, by = y_break)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
################################################################################

## Sequence Logos for model RBP comparison
################################################################################
numMotif = 10

## HH:
data = data.frame(Model_HH)
data = data[order(-data[, 'HH']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$Motif[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## HL:
data = data.frame(Model_HL)
data = data[order(-data[, 'HL']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$Motif[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## LH:
data = data.frame(Model_LH)
data = data[order(-data[, 'LH']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$Motif[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## LL:
data = data.frame(Model_LL)
data = data[order(-data[, 'LL']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$Motif[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

################################################################################

## Make plot for MS
################################################################################
HH_Plot = plot_MS(return_MS(Model_HH))

HL_Plot = plot_MS(return_MS(Model_HL))

LH_Plot = plot_MS(return_MS(Model_LH))

LL_Plot = plot_MS(return_MS(Model_LL))

HH_Plot
HL_Plot
LH_Plot
LL_Plot
################################################################################


## Calculate IS and MS
################################################################################
# nts = c('A', 'C', 'G', 'U')
# K = 5
Models = c('HH', 'HL', 'LH', 'LL')

stats = data.frame(RBP = Models,
                   IS = NA,
                   MS = NA)

HH_IS = return_IS(Model_HH$HH)
HL_IS = return_IS(Model_HL$HL)
LH_IS = return_IS(Model_LH$LH)
LL_IS = return_IS(Model_LL$LL)

stats$IS = c(max(HH_IS), max(HL_IS), max(LH_IS), max(LL_IS))

HH_MS = return_MS(Model_HH)[, c(2:6)]
HH_MS[HH_MS == 0] = NA
HH_MS = mean(c(as.matrix(HH_MS)), na.rm = TRUE)

HL_MS = return_MS(Model_HL)[, c(2:6)]
HL_MS[HL_MS == 0] = NA
HL_MS = mean(c(as.matrix(HL_MS)), na.rm = TRUE)

LH_MS = return_MS(Model_LH)[, c(2:6)]
LH_MS[LH_MS == 0] = NA
LH_MS = mean(c(as.matrix(LH_MS)), na.rm = TRUE)

LL_MS = return_MS(Model_LL)[, c(2:6)]
LL_MS[LL_MS == 0] = NA
LL_MS = mean(c(as.matrix(LL_MS)), na.rm = TRUE)

stats$MS = c(HH_MS, HL_MS, LH_MS, LL_MS)

ggplot() +
  geom_point(data = stats, aes(x = IS, y = MS)) +
  labs(title = "Model RBP", x = "Inherent Specificity", y = "Mutational Sensitivity") +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10))+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


################################################################################


# Define your RNA sequence
################################################################################
RNA = 'UAGCGGUGCGAUUGGCCCGUGGACCGCGUUUUUGCACUCAUCGUUUCGCACUAAGUACAUAUAGUUGCGACAAAGCCGCUUAUGAGUUGGGGGUAUAUUC'

# Create a dataframe
nucleotides_df = data.frame(
  position = 1:nchar(RNA),
  nucleotide = strsplit(RNA, '')[[1]],
  color = factor(
    strsplit(RNA, '')[[1]],
    levels = c('G', 'C', 'A', 'U'),
    labels = c('#CCBB44', '#66CCEE', '#228833', '#AA3377')
  )
)

ggplot(nucleotides_df, aes(x = position, y = 1, fill = color)) +
  geom_tile(width = 1, height = 1, color = "black") +
  scale_fill_identity() +
  theme_void() +
  labs(title = "RNA Sequence Plot") +
  scale_x_continuous(
    position = "top",
    breaks = 1:nchar(RNA),
    labels = strsplit(RNA, '')[[1]]
  ) +
  scale_y_continuous(
    breaks = 1,
    labels = NULL  # Remove y-axis labels
  ) +
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank()  # Remove y-axis ticks
  ) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=10), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

################################################################################


## Handle simulation data
################################################################################
resultsDir = './Data/Model_RBP/sim_results'
R0s = c(10, 25, 50, 100, 250, 500, 1000)
Us = c(29, 30, 31, 32, 33)
Gs = c(89, 90, 91, 92, 93)

simResult_summary = data.frame(
  R0 = R0s,
  UUUUU_HH = NA,
  UUUUU_HL = NA,
  UUUUU_LH = NA,
  UUUUU_LL = NA,
  GGGGG_HH = NA,
  GGGGG_HL = NA,
  GGGGG_LH = NA,
  GGGGG_LL = NA,
  UUUUU_HH_n = NA,
  UUUUU_HL_n = NA,
  UUUUU_LH_n = NA,
  UUUUU_LL_n = NA,
  GGGGG_HH_n = NA,
  GGGGG_HL_n = NA,
  GGGGG_LH_n = NA,
  GGGGG_LL_n = NA
)

for (R0 in R0s) {
  # simResult = fread(paste0(resultsDir, '/PerPos_Equimolar_RBPs_', R0, '_RNA.csv'))
  simResult = fread(paste0(resultsDir, '/PerPos_Lower_High_IS_RBPs_', R0, '_RNA.csv'))
  simResult$HH_norm = simResult$`HH/R0` / mean(simResult$`HH/R0`)
  simResult$HL_norm = simResult$`HL/R0` / mean(simResult$`HL/R0`)
  simResult$LH_norm = simResult$`LH/R0` / mean(simResult$`LH/R0`)
  simResult$LL_norm = simResult$`LL/R0` / mean(simResult$`LL/R0`)
  
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_HH = mean(simResult[Us, ]$`HH/R0`)
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_HL = mean(simResult[Us, ]$`HL/R0`)
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_LH = mean(simResult[Us, ]$`LH/R0`)
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_LL = mean(simResult[Us, ]$`LL/R0`)
  
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_HH = mean(simResult[Gs, ]$`HH/R0`)
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_HL = mean(simResult[Gs, ]$`HL/R0`)
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_LH = mean(simResult[Gs, ]$`LH/R0`)
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_LL = mean(simResult[Gs, ]$`LL/R0`)
  
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_HH_n = mean(simResult[Us, ]$HH_norm)
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_HL_n = mean(simResult[Us, ]$HL_norm)
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_LH_n = mean(simResult[Us, ]$LH_norm)
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_LL_n = mean(simResult[Us, ]$LL_norm)
  
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_HH_n = mean(simResult[Gs, ]$HH_norm)
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_HL_n = mean(simResult[Gs, ]$HL_norm)
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_LH_n = mean(simResult[Gs, ]$LH_norm)
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_LL_n = mean(simResult[Gs, ]$LL_norm)
  
}

ggplot(simResult_summary, aes(x = R0)) +
  geom_line(aes(y = UUUUU_HH_n, color = "HH", linetype = "HH"), linewidth = 1) +
  geom_line(aes(y = UUUUU_HL_n, color = "HL", linetype = "HL"), linewidth = 1) +
  geom_line(aes(y = UUUUU_LH_n, color = "LH", linetype = "LH"), linewidth = 1) +
  geom_line(aes(y = UUUUU_LL_n, color = "LL", linetype = "LL"), linewidth = 1) +
  labs(title = "Line Graph of HH/HL/LH/LL Values with R0 as X-axis") +
  xlab("Initial RNA Concentration") +
  ylab("Relative Binding Probability") +
  # scale_y_continuous(limits = c(0, .5), breaks = seq(0, 2.5, by = 0.1)) +
  scale_x_log10(breaks = R0s,
                labels = R0s) +
  scale_color_manual(values = c("HH" = "#4478ab", "HL" = "#4478ab", "LH" = "#ed6677", "LL" = "#ed6677")) +
  scale_linetype_manual(values = c("HH" = "solid", "HL" = "dashed", "LH" = "dashed", "LL" = "solid")) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=10), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(simResult_summary, aes(x = R0)) +
  geom_line(aes(y = GGGGG_HH_n, color = "HH", linetype = "HH"), linewidth = 1) +
  geom_line(aes(y = GGGGG_HL_n, color = "HL", linetype = "HL"), linewidth = 1) +
  geom_line(aes(y = GGGGG_LH_n, color = "LH", linetype = "LH"), linewidth = 1) +
  geom_line(aes(y = GGGGG_LL_n, color = "LL", linetype = "LL"), linewidth = 1) +
  labs(title = "Line Graph of HH/HL/LH/LL Values with R0 as X-axis") +
  xlab("Initial RNA Concentration") +
  ylab("Relative Binding Probability") +
  # scale_y_continuous(limits = c(0, .5), breaks = seq(0, 2.5, by = 0.1)) +
  scale_x_log10(breaks = R0s,
                labels = R0s) +
  scale_color_manual(values = c("HH" = "#4478ab", "HL" = "#4478ab", "LH" = "#ed6677", "LL" = "#ed6677")) +
  scale_linetype_manual(values = c("HH" = "solid", "HL" = "dashed", "LH" = "dashed", "LL" = "solid")) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=10), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

################################################################################

## Handle simulation data
################################################################################
resultsDir = './Data/Model_RBP/sim_results'
R0s = c(10, 25, 50, 100, 250, 500, 1000)
Us = c(29, 30, 31, 32, 33)
Gs = c(89, 90, 91, 92, 93)

simResult_summary = data.frame(
  R0 = R0s,
  UUUUU_HH = NA,
  UUUUU_HL = NA,
  UUUUU_LH = NA,
  UUUUU_LL = NA,
  GGGGG_HH = NA,
  GGGGG_HL = NA,
  GGGGG_LH = NA,
  GGGGG_LL = NA,
  UUUUU_HH_n = NA,
  UUUUU_HL_n = NA,
  UUUUU_LH_n = NA,
  UUUUU_LL_n = NA,
  GGGGG_HH_n = NA,
  GGGGG_HL_n = NA,
  GGGGG_LH_n = NA,
  GGGGG_LL_n = NA
)

for (R0 in R0s) {
  # simResult = fread(paste0(resultsDir, '/PerPos_Equimolar_RBPs_', R0, '_RNA.csv'))
  simResult = fread(paste0(resultsDir, '/PerPos_Lower_High_IS_RBPs_', R0, '_RNA.csv'))
  simResult$HH_norm = simResult$`HH/R0` / mean(simResult$`HH/R0`)
  simResult$HL_norm = simResult$`HL/R0` / mean(simResult$`HL/R0`)
  simResult$LH_norm = simResult$`LH/R0` / mean(simResult$`LH/R0`)
  simResult$LL_norm = simResult$`LL/R0` / mean(simResult$`LL/R0`)
  
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_HH = mean(simResult[Us, ]$`HH/R0`)
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_HL = mean(simResult[Us, ]$`HL/R0`)
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_LH = mean(simResult[Us, ]$`LH/R0`)
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_LL = mean(simResult[Us, ]$`LL/R0`)
  
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_HH = mean(simResult[Gs, ]$`HH/R0`)
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_HL = mean(simResult[Gs, ]$`HL/R0`)
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_LH = mean(simResult[Gs, ]$`LH/R0`)
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_LL = mean(simResult[Gs, ]$`LL/R0`)
  
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_HH_n = mean(simResult[Us, ]$HH_norm)
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_HL_n = mean(simResult[Us, ]$HL_norm)
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_LH_n = mean(simResult[Us, ]$LH_norm)
  simResult_summary[which(simResult_summary$R0 == R0), ]$UUUUU_LL_n = mean(simResult[Us, ]$LL_norm)
  
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_HH_n = mean(simResult[Gs, ]$HH_norm)
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_HL_n = mean(simResult[Gs, ]$HL_norm)
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_LH_n = mean(simResult[Gs, ]$LH_norm)
  simResult_summary[which(simResult_summary$R0 == R0), ]$GGGGG_LL_n = mean(simResult[Gs, ]$LL_norm)
  
}

ggplot(simResult_summary, aes(x = R0)) +
  geom_line(aes(y = UUUUU_HH, color = "HH", linetype = "HH"), linewidth = 1) +
  geom_line(aes(y = UUUUU_HL, color = "HL", linetype = "HL"), linewidth = 1) +
  geom_line(aes(y = UUUUU_LH, color = "LH", linetype = "LH"), linewidth = 1) +
  geom_line(aes(y = UUUUU_LL, color = "LL", linetype = "LL"), linewidth = 1) +
  labs(title = "Line Graph of HH/HL/LH/LL Values with R0 as X-axis") +
  xlab("Initial RNA Concentration") +
  ylab("Relative Binding Probability") +
  # scale_y_continuous(limits = c(0, .5), breaks = seq(0, 2.5, by = 0.1)) +
  scale_x_log10(breaks = R0s,
                labels = R0s) +
  scale_color_manual(values = c("HH" = "#4478ab", "HL" = "#4478ab", "LH" = "#ed6677", "LL" = "#ed6677")) +
  scale_linetype_manual(values = c("HH" = "solid", "HL" = "dashed", "LH" = "dashed", "LL" = "solid")) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=10), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(simResult_summary, aes(x = R0)) +
  geom_line(aes(y = GGGGG_HH, color = "HH", linetype = "HH"), linewidth = 1) +
  geom_line(aes(y = GGGGG_HL, color = "HL", linetype = "HL"), linewidth = 1) +
  geom_line(aes(y = GGGGG_LH, color = "LH", linetype = "LH"), linewidth = 1) +
  geom_line(aes(y = GGGGG_LL, color = "LL", linetype = "LL"), linewidth = 1) +
  labs(title = "Line Graph of HH/HL/LH/LL Values with R0 as X-axis") +
  xlab("Initial RNA Concentration") +
  ylab("Relative Binding Probability") +
  # scale_y_continuous(limits = c(0, .5), breaks = seq(0, 2.5, by = 0.1)) +
  scale_x_log10(breaks = R0s,
                labels = R0s) +
  scale_color_manual(values = c("HH" = "#4478ab", "HL" = "#4478ab", "LH" = "#ed6677", "LL" = "#ed6677")) +
  scale_linetype_manual(values = c("HH" = "solid", "HL" = "dashed", "LH" = "dashed", "LL" = "solid")) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=10), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

################################################################################


## Per Pos Binding Probaability
################################################################################
resultsDir = './Data/Model_RBP/sim_results'

simResult_EquiMolar = fread(paste0(resultsDir, '/PerPos_Equimolar_RBPs_10_RNA.csv'))
simResult_EquiMolar = simResult_EquiMolar[, c('pos', 'nt', "HH/R0", "HL/R0", "LH/R0", "LL/R0")]

ggplot(simResult_EquiMolar, aes(x = pos)) +
  geom_line(aes(y = `HH/R0`, color = "HH", linetype = "HH"), linewidth = 0.5) +
  geom_line(aes(y = `HL/R0`, color = "HL", linetype = "HL"), linewidth = 0.5) +
  geom_line(aes(y = `LH/R0`, color = "LH", linetype = "LH"), linewidth = 0.5) +
  geom_line(aes(y = `LL/R0`, color = "LL", linetype = "LL"), linewidth = 0.5) +
  labs(title = "Line Graph of HH/HL/LH/LL Values with R0 as X-axis") +
  xlab("Initial RNA Concentration") +
  ylab("Relative Binding Probability") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(
    position = "bottom",
    breaks = 1:nchar(RNA),
    labels = strsplit(RNA, '')[[1]]
  ) +
  scale_color_manual(values = c("HH" = "#4478ab", "HL" = "#4478ab", "LH" = "#ed6677", "LL" = "#ed6677")) +
  scale_linetype_manual(values = c("HH" = "solid", "HL" = "dashed", "LH" = "dashed", "LL" = "solid")) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=10), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

simResult_DiffMolar = fread(paste0(resultsDir, '/PerPos_Lower_High_IS_RBPs_250_RNA.csv'))
simResult_DiffMolar = simResult_DiffMolar[, c('pos', 'nt', "HH/R0", "HL/R0", "LH/R0", "LL/R0")]

ggplot(simResult_DiffMolar, aes(x = pos)) +
  geom_line(aes(y = `HH/R0`, color = "HH", linetype = "HH"), linewidth = 0.5) +
  geom_line(aes(y = `HL/R0`, color = "HL", linetype = "HL"), linewidth = 0.5) +
  geom_line(aes(y = `LH/R0`, color = "LH", linetype = "LH"), linewidth = 0.5) +
  geom_line(aes(y = `LL/R0`, color = "LL", linetype = "LL"), linewidth = 0.5) +
  labs(title = "Line Graph of HH/HL/LH/LL Values with R0 as X-axis") +
  xlab("Initial RNA Concentration") +
  ylab("Relative Binding Probability") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(
    position = "bottom",
    breaks = 1:nchar(RNA),
    labels = strsplit(RNA, '')[[1]]
  ) +
  scale_color_manual(values = c("HH" = "#4478ab", "HL" = "#4478ab", "LH" = "#ed6677", "LL" = "#ed6677")) +
  scale_linetype_manual(values = c("HH" = "solid", "HL" = "dashed", "LH" = "dashed", "LL" = "solid")) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=10), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


################################################################################

## Per Pos Binding Probaability
################################################################################
library(pracma)

simResult_EquiMolar = fread('./Data/plot_data.csv')
simResult_EquiMolar = simResult_EquiMolar[, c('pos', 'nt', "HHReq", "HLReq", "LHReq", "LLReq")]

HH0 = 100
HL0 = 100
LH0 = 100
LL0 = 100

RNA_len = nchar(RNA)

simResult_EquiMolar$HHReq = (simResult_EquiMolar$HHReq / trapz(simResult_EquiMolar$pos, simResult_EquiMolar$HHReq)) / (HH0/RNA_len) * 100
simResult_EquiMolar$HLReq = (simResult_EquiMolar$HLReq / trapz(simResult_EquiMolar$pos, simResult_EquiMolar$HLReq)) / (HL0/RNA_len) * 100
simResult_EquiMolar$LHReq = (simResult_EquiMolar$LHReq / trapz(simResult_EquiMolar$pos, simResult_EquiMolar$LHReq)) / (LH0/RNA_len) * 100
simResult_EquiMolar$LLReq = (simResult_EquiMolar$LLReq / trapz(simResult_EquiMolar$pos, simResult_EquiMolar$LLReq)) / (LL0/RNA_len) * 100


ggplot(simResult_EquiMolar, aes(x = pos)) +
  geom_line(aes(y = `HHReq`, color = "HH", linetype = "HH"), linewidth = 0.5) +
  geom_line(aes(y = `HLReq`, color = "HL", linetype = "HL"), linewidth = 0.5) +
  geom_line(aes(y = `LHReq`, color = "LH", linetype = "LH"), linewidth = 0.5) +
  geom_line(aes(y = `LLReq`, color = "LL", linetype = "LL"), linewidth = 0.5) +
  labs(title = "Line Graph of HH/HL/LH/LL Values with R0 as X-axis") +
  xlab("Initial RNA Concentration") +
  ylab("Relative Binding Probability") +
  # scale_y_continuous(trans = log2_trans(), limits = c(0.1, 4),
  #                    breaks = c(0.1, 0.5, 1, 2, 4),
  #                    labels = c("0.1", '0.5', "1", "2", "4")) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 0.5)) +
  scale_x_continuous(
    position = "bottom",
    breaks = 1:nchar(RNA),
    labels = strsplit(RNA, '')[[1]]
  ) +
  scale_color_manual(values = c("HH" = "#4478ab", "HL" = "#4478ab", "LH" = "#ed6677", "LL" = "#ed6677")) +
  scale_linetype_manual(values = c("HH" = "solid", "HL" = "dashed", "LH" = "dashed", "LL" = "solid")) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=10), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))




################################################################################




