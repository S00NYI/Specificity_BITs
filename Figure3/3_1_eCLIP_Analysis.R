################################################################################
## eCLIP Analysis V5: 
## Written by Soon Yi
## Created: 2023-12-08
## Last Edited: 2024-02-05
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

## Custom Functions:
################################################################################
# Min-Max normalization:
min_max_norm = function(x, a = 0, b = 1) {
  a + (x - min(x, na.rm = TRUE)) * (b - a) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Specificity index calculation:
SI_Calculation = function(x) {
  x / (median(x, na.rm = FALSE))
}

# Scramble given sequence
ScrambleDNA = function(seq) {
  scrambled_seq = sample(unlist(strsplit(as.character(seq), "")))
  return(DNAString(paste(scrambled_seq, collapse = "")))
}

# Generate N random distances with a minimum of X nt away
generate_random_distances = function(N, X) {
  distances_negative = runif(N / 2, min = -1000, max = -X)
  distances_positive = runif(N / 2, min = X, max = 1000)
  distances = sample(c(distances_negative, distances_positive))
  return(distances)
}

################################################################################

## Set up basic parameters:
################################################################################
# baseDir = '~/Desktop/Genomics/Specificity/'

baseDir = '~/Desktop/Genomics/Specificity/Data_Final/'
sampleTable_eCLIP = read_csv(paste0(baseDir, 'sample_table_eCLIP.csv'), col_names = T, show_col_types = F)
sampleTable_eCLIP = data.frame(sampleTable_eCLIP)

cell_lines = c('K562', 'HepG2')

K562_RBPs = as.character(unique((sampleTable_eCLIP %>% filter(Cell == 'K562'))[, 'RBP']))
HepG2_RBPs = as.character(unique((sampleTable_eCLIP %>% filter(Cell == 'HepG2'))[, 'RBP']))

DATE = '20240217'
extension = 25

K = 5
NTs = c("A", "C", "G", "T")
KMer_combinations = expand.grid(replicate(K, NTs, simplify = FALSE))
KMer_combinations = apply(KMer_combinations, 1, paste, collapse = "")
KMer_DSS = DNAStringSet(KMer_combinations)
KMer_PDict = PDict(KMer_DSS)
################################################################################

################################################################################
for (cell_line in cell_lines) {
  print(paste0('Working on ', cell_line, '....'))
  if (cell_line == 'K562') {RBPs = K562_RBPs} else {RBPs = HepG2_RBPs}
  for (RBP in RBPs) {
    # Read peaks and make GR object
    print(paste0('Working on ', RBP, '....'))
    print(paste0('     Reading and processing peak data....'))
    file_ID = sampleTable_eCLIP[sampleTable_eCLIP$RBP == RBP & sampleTable_eCLIP$Cell == cell_line, ]$eCLIP
    Peak = fread(paste0(baseDir, 'eCLIP_peak/raw/', file_ID,  '.bed'))
    Peak = Peak %>% mutate(V1 = case_when(V1 == 'chrMT' ~ 'chrM', TRUE ~ as.character(V1)))
    Peak = Peak %>% filter(V1 %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 
                                     'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
    Peak_GR = GRanges(seqnames = Rle(Peak$V1),
                      ranges = IRanges(start = Peak$V2 - extension, end = Peak$V3),
                      strand = Rle(Peak$V6))
    
    # Get sequence of peaks and count K-mers
    Peak_Seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, Peak_GR, as.character = TRUE)
    Peak_Seqs = DNAStringSet(Peak_Seqs)
    Peak_KMers = data.frame(MOTIF = KMer_combinations, 
                            COUNT = vcountPDict(KMer_PDict, Peak_Seqs, collapse = 1))
    
    # Make background set of peaks that are random distance away from the peaks, but with the same size.
    # We will do random selection 100 times, count K-mers, and take average
    print(paste0('     Making and processing background data....'))
    iterations = 100
    CNT = data.frame(MOTIF = KMer_combinations)
    AVG = data.frame(MOTIF = KMer_combinations)
    
    pb = txtProgressBar(min = 0, max = iterations, style = 3)
    
    for (idx in c(1:iterations)) {
      
      Rand_Dist = generate_random_distances(length(Peak_GR), 500)
      Rand_GR = shift(Peak_GR, round(Rand_Dist, 0))
      overlaps = findOverlaps(Rand_GR, Peak_GR)
      Rand_GR = Rand_GR[-queryHits(overlaps)]
      
      Rand_Seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, Rand_GR, as.character = TRUE)
      Rand_Seqs = DNAStringSet(lapply(Rand_Seqs, ScrambleDNA))
      Bkg_KMer = data.frame(MOTIF = KMer_combinations,
                            COUNT = vcountPDict(KMer_PDict, Rand_Seqs, collapse = 1))
      
      CNT = cbind(CNT, (Bkg_KMer$COUNT))
      
      if (idx == 1) {
        AVG$AVG = CNT[, -1]
      } else {
        AVG = cbind(AVG, data.frame(rowMeans(CNT[, -1])))
      }
      setTxtProgressBar(pb, idx)
      
    }
    print(paste0(""))
    colnames(AVG) = c('MOTIF', paste0('Iter_', c(1:iterations)))
    # write.csv(AVG, 
    #           paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/0_Bkg/', RBP, '_', K, 'mer.csv'), quote = F)
    
    # Set the 100th average as the background K-mer counts:
    Bkg_KMer$COUNT = AVG[, length(colnames(AVG))]
    
    # First, calculate the background corrected K-mer counts:
    Peak_KMers_C = Peak_KMers
    # Peak_KMers_C$COUNT = Peak_KMers_C$COUNT - Bkg_KMer$COUNT
    Peak_KMers_C$COUNT = Peak_KMers_C$COUNT / Bkg_KMer$COUNT
    
    # Second, scale the normalized counts to scale of 0-to-1:
    KMer_Matches_Scaled = Peak_KMers_C
    # KMer_Matches_Scaled = KMer_Matches_Scaled %>% mutate_if(is.numeric, min_max_norm) %>% mutate_if(is.numeric, round, 10)
    KMer_Matches_Scaled = KMer_Matches_Scaled %>% mutate_if(is.numeric, min_max_norm, a = 1, b = exp(1)) %>% mutate_if(is.numeric, round, 10)
    KMer_Matches_Scaled$COUNT = log(KMer_Matches_Scaled$COUNT)

    # Third, calculate the specificity index for each K-mer:
    KMer_SpecificityIndex = KMer_Matches_Scaled
    KMer_SpecificityIndex = KMer_Matches_Scaled %>% mutate_if(is.numeric, SI_Calculation) %>% mutate_if(is.numeric, round, 10)
    
    print(paste0('     Writing all Kmer information....'))
    write.csv(Peak_KMers_C, 
              paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/1_Raw/', RBP, '_', K, 'mer.csv'), quote = F)
    write.csv(KMer_Matches_Scaled, 
              paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/2_Scaled/', RBP, '_', K, 'mer.csv'), quote = F)
    write.csv(KMer_SpecificityIndex, 
              paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', cell_line, '/3_SpecificityIndex/', RBP, '_', K, 'mer.csv'), quote = F)
  }
}
################################################################################
















