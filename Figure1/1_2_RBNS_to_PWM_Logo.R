################################################################################
## Generate PWM/Logo from RBNS 
## Written by Soon Yi
## Created: 2023-10-20
## Last Edited: 2024-01-31
## Figure 1
################################################################################

library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(Biostrings)
library(seqLogo)


## Set some basic parameters:
################################################################################
baseDir = '~/Desktop/Genomics/Specificity/Data_Final/'

K = 5
RBNS = read_csv(paste0(baseDir, 'RBNS/RBNS_normalized_', K, 'mer.csv'), col_names = T, show_col_types = F)
colnames(RBNS) = c('Motif', colnames(RBNS)[2:ncol(RBNS)])
RBPs = colnames(RBNS)[2:ncol(RBNS)]
################################################################################

## For hnRNPC and EIF4G2:
################################################################################
numMotif = 10

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
seqLogo(PWM, ic.scale = F)


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
seqLogo(PWM, ic.scale = F)
################################################################################





