################################################################################
## RBNS Data Parsing
## Written by Soon Yi
## Last Edited: 2024-02-16
################################################################################

library(readr)
library(dplyr)
library(tidyr)

## Define some functions:
################################################################################
# Min-Max Normalization:
min_max_norm = function(x, a = 0, b = 1) {
  a + (x - min(x, na.rm = TRUE)) * (b - a) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# RBNS Processing:
process_RBNS = function(K, input_df, special_org, special_map) {
  RBPs = unique(input_df$RBP)
  
  Kmer = expand.grid(replicate(K, NTs, simplify = FALSE))
  Kmer = apply(Kmer, 1, paste, collapse = "")
  
  output_df = data.frame(Motif = Kmer)
  output_df = output_df %>% arrange(Motif)
  
  for (RBP in RBPs) {
    # Read in the RBNS data for each RBP
    RBP_ID = input_df$RBNS[which(input_df$RBP == RBP)]
    RBP_Data = read_tsv(paste0(baseDir, 'RBNS/raw/', RBP, '/', RBP_ID, '.tsv'), col_names = T, show_col_types = F)
    columns = colnames(RBP_Data)
    RBP_Data = data.frame(RBP_Data)
    colnames(RBP_Data) = columns
    
    # Check if RBP is a special case
    RBP_org = RBP
    if (RBP %in% special_org) {
      RBP = special_map[match(RBP, special_org)]
    }
    
    # Col 1, Row 1 value contains the RBP name, use it to check if we are loading right RBP.
    assigned_name = gsub("[", "", gsub("]", "", colnames(RBP_Data)[1], fixed = T), fixed = T)
    assigned_name = toupper(assigned_name)
    assigned_name = strsplit(assigned_name, '_')[[1]][1]
    
    if (assigned_name == RBP) {
      print('RBP name in the sample list matches assigned name in the data.')
      print(paste0(RBP_org, ' vs ', assigned_name))
    } else {
      print('RBP name in the sample list does not match assigned name in the data.')
      print(paste0(RBP_org, ' vs ', assigned_name))
      break
    }
    
    RBP = RBP_org
    
    # Check number of rows to see if we are loading right K-mer data.
    if (nrow(RBP_Data) != nrow(output_df)) {
      print(paste0('Number of motif does not match for ', RBP))
      print(paste0('We need: ', nrow(output_df)))
      print(paste0('We have: ', nrow(RBP_Data)))
      if (nrow(RBP_Data) %in% c(4^4, 4^5, 4^6, 4^7)) {
        print('Wrong K-mer data has been imported, check the accession number.')
        break
      } else {
        print('We have missing entries: ')
        missing_Kmers = output_df$Motif[!(output_df$Motif %in% RBP_Data[, 1])]
        print(missing_Kmers)
        print('Filling them with 0s.')
        filler = data.frame(matrix(data = 0, nrow = length(missing_Kmers), ncol = length(columns)))
        colnames(filler) = columns
        filler[, 1] = missing_Kmers
        RBP_Data = rbind(RBP_Data, filler)
      }
    }
    
    colOriginal = colnames(RBP_Data)[-1]
    colnames(RBP_Data) = c('Motif', c(1:(length(colnames(RBP_Data))-1)))
    
    RBP_Data$Motif = sapply(RBP_Data$Motif, function(x) gsub("[^ACGT]", "", x, perl=TRUE))
    RBP_Data = RBP_Data %>% arrange(Motif)
    
    colMax = apply(RBP_Data[, -1], MARGIN = 2, max, na.rm = TRUE)
    colChoice = which(colMax == max(colMax))
    
    print(paste0('Column with maximum R: ', colOriginal[colChoice]))
    RBP_Data = RBP_Data[, c('Motif', colChoice)]
    colnames(RBP_Data) = c('Motif', RBP)
    RBP_Data[, RBP] = data.frame(log(min_max_norm(RBP_Data[, RBP], a = 1, b = exp(1))))
    
    output_df = cbind(output_df, setNames(RBP_Data[, RBP, drop = FALSE], RBP))
  }
  output_df$Motif = gsub("T", "U", output_df$Motif)
  return(output_df)
}

################################################################################

## Set some basic parameters:
################################################################################
baseDir = '~/Desktop/Genomics/Specificity/Data_Final/'
sampleTable_RBNS = read_csv(paste0(baseDir, 'sample_table_RBNS.csv'), col_names = T, show_col_types = F)

RBNS_4mer = sampleTable_RBNS %>% filter(`K-mer` == '4-mer')
RBNS_5mer = sampleTable_RBNS %>% filter(`K-mer` == '5-mer')
RBNS_6mer = sampleTable_RBNS %>% filter(`K-mer` == '6-mer')
RBNS_7mer = sampleTable_RBNS %>% filter(`K-mer` == '7-mer')

NTs = c("A", "C", "G", "T")

special_case_org = c('MBNL1', 'RBFOX2', 'TRNAU1AP', 'ZFP36L1', 'ZFP36L2')
special_case_map = c('MBNL', 'FOX', 'SECP43', 'TIS11B', 'TIS11D')
################################################################################

## Parse 4-mer data:
################################################################################
Data_4mer = process_RBNS(4, RBNS_4mer, special_case_org, special_case_map)

write.csv(Data_4mer, paste0(baseDir, 'RBNS/RBNS_normalized_4mer.csv'), quote = F, row.names = F)

RBNS_4mer_singleRBD = RBNS_4mer %>% filter(singleRBD == 'Y')
Data_4mer_singleRBD = Data_4mer[, c('Motif', RBNS_4mer_singleRBD$RBP)]
write.csv(Data_4mer_singleRBD, paste0(baseDir, 'RBNS/RBNS_singleRBD_normalized_4mer.csv'), quote = F, row.names = F)

################################################################################

## Parse 5-mer data:
################################################################################
Data_5mer = process_RBNS(5, RBNS_5mer, special_case_org, special_case_map)

write.csv(Data_5mer, paste0(baseDir, 'RBNS/RBNS_normalized_5mer.csv'), quote = F, row.names = F)

RBNS_5mer_singleRBD = RBNS_5mer %>% filter(singleRBD == 'Y')
Data_5mer_singleRBD = Data_5mer[, c('Motif', RBNS_5mer_singleRBD$RBP)]
write.csv(Data_5mer_singleRBD, paste0(baseDir, 'RBNS/RBNS_singleRBD_normalized_5mer.csv'), quote = F, row.names = F)

################################################################################

## Parse 6-mer data:
################################################################################
Data_6mer = process_RBNS(6, RBNS_6mer, special_case_org, special_case_map)

write.csv(Data_6mer, paste0(baseDir, 'RBNS/RBNS_normalized_6mer.csv'), quote = F, row.names = F)

RBNS_6mer_singleRBD = RBNS_6mer %>% filter(singleRBD == 'Y')
Data_6mer_singleRBD = Data_6mer[, c('Motif', RBNS_6mer_singleRBD$RBP)]
write.csv(Data_6mer_singleRBD, paste0(baseDir, 'RBNS/RBNS_singleRBD_normalized_6mer.csv'), quote = F, row.names = F)

################################################################################

## Parse 7-mer data:
################################################################################
Data_7mer = process_RBNS(7, RBNS_7mer, special_case_org, special_case_map)

write.csv(Data_7mer, paste0(baseDir, 'RBNS/RBNS_normalized_7mer.csv'), quote = F, row.names = F)

RBNS_7mer_singleRBD = RBNS_7mer %>% filter(singleRBD == 'Y')
Data_7mer_singleRBD = Data_7mer[, c('Motif', RBNS_7mer_singleRBD$RBP)]
write.csv(Data_7mer_singleRBD, paste0(baseDir, 'RBNS/RBNS_singleRBD_normalized_7mer.csv'), quote = F, row.names = F)

################################################################################