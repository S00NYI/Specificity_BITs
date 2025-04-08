################################################################################
## Compare RBNS to PWM
## Written by Soon Yi
## Created: 2023-11-08
## Last Edited: 2024-02-11
## Figure 1
################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Set some basic parameters:
################################################################################
baseDir = '~/Desktop/Genomics/Specificity/Data_Final/'

# Read RBNS file:
K = 5

nts_id = c('A', 'C', 'G', 'U')
Kmers = expand.grid(rep(list(nts_id), K))
colnames(Kmers) = c(1, 2, 3, 4, 5)

RBNS = read_csv(paste0(baseDir, 'RBNS/RBNS_normalized_', K, 'mer.csv'), col_names = T, show_col_types = F)
colnames(RBNS) = c('Motif', colnames(RBNS)[2:ncol(RBNS)])
RBPs = colnames(RBNS)[2:ncol(RBNS)]
################################################################################

################################################################################
numMotif = 20

RBP = 'HNRNPC'

data = data.frame(RBNS[, c('Motif', RBP)])
data = data[order(-data[, RBP]), ]
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

# Container for scores:
Scores = matrix(NA, nrow = 1024, ncol = 5)
colnames(Scores) = c(1, 2, 3, 4, 5)

# Calculate scores per motif:
for (id in 1:nrow(Kmers)) {
  Kmer = Kmers[id, ]
  for (pos in 1:ncol(Kmer)) {
    Scores[id, pos] = PWM@pwm[Kmer[[pos]], pos]
  }
}

PWM_Motif_Scores = data.frame(Motif = apply(Kmers, 1, paste, collapse = ""),
                              Score = rowSums(Scores))

PWM_Motif_Scores_Normalized = PWM_Motif_Scores
PWM_Motif_Scores_Normalized$Score = (PWM_Motif_Scores_Normalized$Score - min(PWM_Motif_Scores_Normalized$Score)) / (max(PWM_Motif_Scores_Normalized$Score) - min(PWM_Motif_Scores_Normalized$Score))

data = data %>% arrange(Motif)
PWM_Motif_Scores_Normalized = PWM_Motif_Scores_Normalized %>% arrange(Motif)

Scores_combined = cbind(data, PWM_Motif_Scores_Normalized$Score)
colnames(Scores_combined) = c('Motif', 'RBNS', 'PWM')

ggplot(Scores_combined, aes(x = RBNS, y = PWM)) +
  geom_point() +  
  theme_minimal() + 
  labs(x = "RBNS Score", y = "PWM Score", title = "Scatter Plot of RBNS vs. PWM") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

################################################################################


################################################################################
numMotif = 10

RBP = 'EIF4G2'

data = data.frame(RBNS[, c('Motif', RBP)])
data = data[order(-data[, RBP]), ]
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

# Container for scores:
Scores = matrix(NA, nrow = 1024, ncol = 5)
colnames(Scores) = c(1, 2, 3, 4, 5)

# Calculate scores per motif:
for (id in 1:nrow(Kmers)) {
  Kmer = Kmers[id, ]
  for (pos in 1:ncol(Kmer)) {
    Scores[id, pos] = PWM@pwm[Kmer[[pos]], pos]
  }
}

PWM_Motif_Scores = data.frame(Motif = apply(Kmers, 1, paste, collapse = ""),
                              Score = rowSums(Scores))

PWM_Motif_Scores_Normalized = PWM_Motif_Scores
PWM_Motif_Scores_Normalized$Score = (PWM_Motif_Scores_Normalized$Score - min(PWM_Motif_Scores_Normalized$Score)) / (max(PWM_Motif_Scores_Normalized$Score) - min(PWM_Motif_Scores_Normalized$Score))

data = data %>% arrange(Motif)
PWM_Motif_Scores_Normalized = PWM_Motif_Scores_Normalized %>% arrange(Motif)

Scores_combined = cbind(data, PWM_Motif_Scores_Normalized$Score)
colnames(Scores_combined) = c('Motif', 'RBNS', 'PWM')

ggplot(Scores_combined, aes(x = RBNS, y = PWM)) +
  geom_point() +  
  theme_minimal() + 
  labs(x = "RBNS Score", y = "PWM Score", title = "Scatter Plot of RBNS vs. PWM") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

################################################################################

################################################################################
numMotif = 20

RBP = 'RBM25'

data = data.frame(RBNS[, c('Motif', RBP)])
data = data[order(-data[, RBP]), ]
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

# Container for scores:
Scores = matrix(NA, nrow = 1024, ncol = 5)
colnames(Scores) = c(1, 2, 3, 4, 5)

# Calculate scores per motif:
for (id in 1:nrow(Kmers)) {
  Kmer = Kmers[id, ]
  for (pos in 1:ncol(Kmer)) {
    Scores[id, pos] = PWM@pwm[Kmer[[pos]], pos]
  }
}

PWM_Motif_Scores = data.frame(Motif = apply(Kmers, 1, paste, collapse = ""),
                              Score = rowSums(Scores))

PWM_Motif_Scores_Normalized = PWM_Motif_Scores
PWM_Motif_Scores_Normalized$Score = (PWM_Motif_Scores_Normalized$Score - min(PWM_Motif_Scores_Normalized$Score)) / (max(PWM_Motif_Scores_Normalized$Score) - min(PWM_Motif_Scores_Normalized$Score))

data = data %>% arrange(Motif)
PWM_Motif_Scores_Normalized = PWM_Motif_Scores_Normalized %>% arrange(Motif)

Scores_combined = cbind(data, PWM_Motif_Scores_Normalized$Score)
colnames(Scores_combined) = c('Motif', 'RBNS', 'PWM')

ggplot(Scores_combined, aes(x = RBNS, y = PWM)) +
  geom_point() +  
  theme_minimal() + 
  labs(x = "RBNS Score", y = "PWM Score", title = "Scatter Plot of RBNS vs. PWM") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

################################################################################