################################################################################
## Motif Analysis: 
## Written by Soon Yi
## Created: 2025-01-24
## Last Edited: 2025-03-01
## Figure 4
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

return_MS = function(data, top_motif = NULL) {
  colnames(data) = c('Motif', 'Score')
  MS = data.frame(Nucleotide = c('A', 'C', 'G', 'U'),
                  Pos1 = NA,
                  Pos2 = NA,
                  Pos3 = NA,
                  Pos4 = NA,
                  Pos5 = NA)
  
  if (!is.null(top_motif)) {
    Motif_Top = top_motif
  } else {
    Motif_Top = data$Motif[which(data[, 'Score'] == max(data[, 'Score']))]
  }
  
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
################################################################################

## Set up basic parameters:
################################################################################
baseDir = '~/Desktop/Genomics/Specificity/Sequencing/Processed/IndividualPeaks'

# previous figures used 25
# may need to optimize extension 10x, 25x, 50x
# extension = 25
extensions = c(0, 5, 10, 25, 50, 100)

K = 5
NTs = c("A", "C", "G", "T")
KMer_combinations = expand.grid(replicate(K, NTs, simplify = FALSE))
KMer_combinations = apply(KMer_combinations, 1, paste, collapse = "")
KMer_DSS = DNAStringSet(KMer_combinations)
KMer_PDict = PDict(KMer_DSS)
################################################################################

## Filter peaks and generate isolated files for filtered peaks:
################################################################################
## Filters
# HepG2: HepG2_BC == 2, HepG2_avg_RPM > 0.3
# K562: K562_BC == 2, K562_avg_RPM > 0.3


# eCLIP: Total_BC == 4, HepG2_avg_RPM > 0.2, K562_avg_RPM > 0.2
# Hela: Hela >= 2, hela_avg_RPM > 2.0
# Huh7: 
# hnRNPC_WT_BC >= 3, hnRNPC_WT_avg_RPM > 0.2
# hnRNPC_Mut_BC >= 3, hnRNPC_Mut_avg_RPM > 0.2
# RBM25_WT_BC >= 2, RBM25_WT_avg_RPM > 0.1
# RBM25_Mut_BC >= 2, RBM25_Mut_avg_RPM > 0.1

## eCLIP
Peak_eCLIP = fread(paste0(baseDir, '/eCLIP_peaks_filtered.txt'))
Peak_eCLIP = Peak_eCLIP %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_eCLIP = Peak_eCLIP %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                  'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))


## K562
Peak_K562 = fread(paste0(baseDir, '/K562_peaks_filtered.txt'))
Peak_K562 = Peak_K562 %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_K562 = Peak_K562 %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                  'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))

## HepG2
Peak_HepG2 = fread(paste0(baseDir, '/HepG2_peaks_filtered.txt'))
Peak_HepG2 = Peak_HepG2 %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_HepG2 = Peak_HepG2 %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                  'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))

## Hela
Peak_HeLa = fread(paste0(baseDir, '/HeLa_peaks_filtered.txt'))
Peak_HeLa = Peak_HeLa %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_HeLa = Peak_HeLa %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                  'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))


## Our data
Peak_hnRNPC_WT = fread(paste0(baseDir, '/Huh7_WT_HNRNPC_peaks_filtered.txt'))
Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                  'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))

Peak_hnRNPC_Mut = fread(paste0(baseDir, '/Huh7_Mut_HNRNPC_peaks_filtered.txt'))
Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                                      'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))

Peak_RBM25_WT = fread(paste0(baseDir, '/Huh7_WT_RBM25_peaks_filtered.txt'))
Peak_RBM25_WT = Peak_RBM25_WT %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_RBM25_WT = Peak_RBM25_WT %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                                      'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))

Peak_RBM25_Mut = fread(paste0(baseDir, '/Huh7_Mut_RBM25_peaks_filtered.txt'))
Peak_RBM25_Mut = Peak_RBM25_Mut %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_RBM25_Mut = Peak_RBM25_Mut %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                                      'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
################################################################################

## Motif Analysis
################################################################################
samples = c('eCLIP_hnRNPC_WT', 
            'K562_hnRNP_WT', 'HepG2_hnRNPC_WT', 
            'HeLa_hnRNPC_WT', 
            'Huh75_hnRNPC_WT', 'Huh75_hnRNPC_Mut', 
            'Huh75_RBM25_WT', 'Huh75_RBM25_Mut')

UUUUU_CLIP = c()
UUUUU_Ext = c()
UUUUU_SI = c()
UUUUU_Rank = c()

GGGGG_SI = c()
GGGGG_Rank = c()


DATE = format(Sys.Date(), "%Y%m%d")
if (!dir.exists(paste0(baseDir, '/Output/', DATE, '/'))) {
  dir.create(paste0(baseDir, '/Output/', DATE, '/'), recursive = TRUE, showWarnings = FALSE)
  print(paste("Created directory:", paste0(baseDir, '/Output/', DATE, '/')))
}

for (sample in samples) {
  # Read peaks and make GR object
  print(paste0('Working on ', sample, '....'))
  print(paste0('     Reading and processing peak data....'))
  
  if (sample == 'eCLIP_hnRNPC_WT') {
    Peak = Peak_eCLIP
  } else if (sample == 'K562_hnRNPC_WT') {
    Peak = Peak_K562
  } else if (sample == 'HepG2_hnRNPC_WT') {
    Peak = Peak_HepG2
  } else if (sample == 'HeLa_hnRNPC_WT') {
    Peak = Peak_HeLa
  } else if (sample == 'Huh75_hnRNPC_WT') {
    Peak = Peak_hnRNPC_WT
  } else if (sample == 'Huh75_hnRNPC_Mut') {
    Peak = Peak_hnRNPC_Mut
  } else if (sample == 'Huh75_RBM25_WT') {
    Peak = Peak_RBM25_WT
  } else if (sample == 'Huh75_RBM25_Mut') {
    Peak = Peak_RBM25_Mut
  }
  
  for (extension in extensions) {
    print(paste0('     Working on ', extension, 'nt extension....'))
    
    ######
    # Peak_GR = GRanges(seqnames = Rle(Peak$chr),
    #                   ranges = IRanges(start = Peak$start - extension, end = Peak$end),
    #                   strand = Rle(Peak$strand))
    
    Peak_GR_original = GRanges(seqnames = Rle(Peak$chr),
                               ranges = IRanges(start = Peak$start - extension, end = Peak$end),
                               strand = Rle(Peak$strand))
    original_width = width(Peak_GR_original)
    new_width = original_width + extension
    Peak_GR = resize(Peak_GR_original, width = new_width, fix = "end")
    ######
    
    # Get sequence of peaks and count K-mers
    Peak_Seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, Peak_GR, as.character = TRUE)
    Peak_Seqs = DNAStringSet(Peak_Seqs)
    Peak_KMers = data.frame(MOTIF = KMer_combinations, 
                            COUNT = vcountPDict(KMer_PDict, Peak_Seqs, collapse = 1))
    
    # Make background set of peaks that are random distance away from the peaks, but with the same size.
    # We will do random selection 100 times, count K-mers, and take average
    print(paste0('          Making and processing background data....'))
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
    colnames(AVG) = c('MOTIF', paste0('Iter_', c(1:iterations)))
    
    # Set the 100th average as the background K-mer counts:
    Bkg_KMer$COUNT = AVG[, length(colnames(AVG))]
    
    # First, calculate the background corrected K-mer counts:
    Peak_KMers_C = Peak_KMers
    Peak_KMers_C$COUNT = Peak_KMers_C$COUNT - Bkg_KMer$COUNT

    # Second, scale the normalized counts to scale of 0-to-1:
    KMer_Matches_Scaled = Peak_KMers_C
    KMer_Matches_Scaled = KMer_Matches_Scaled %>% mutate_if(is.numeric, min_max_norm, a = 1, b = exp(1)) %>% mutate_if(is.numeric, round, 10)
    KMer_Matches_Scaled$COUNT = log(KMer_Matches_Scaled$COUNT)
    
    # Third, calculate the specificity index for each K-mer:
    KMer_SpecificityIndex = KMer_Matches_Scaled
    KMer_SpecificityIndex = KMer_Matches_Scaled %>% mutate_if(is.numeric, SI_Calculation) %>% mutate_if(is.numeric, round, 10)
    
    
    UUUUU_CLIP = c(UUUUU_CLIP, sample)
    UUUUU_Ext = c(UUUUU_Ext, extension)
    UUUUU_SI = c(UUUUU_SI, as.numeric(formatC(KMer_SpecificityIndex$COUNT[KMer_SpecificityIndex$MOTIF=='TTTTT'], digits = 2, format = "f")))
    GGGGG_SI = c(GGGGG_SI, as.numeric(formatC(KMer_SpecificityIndex$COUNT[KMer_SpecificityIndex$MOTIF=='GGGGG'], digits = 2, format = "f")))
    
    ordered_SI = KMer_SpecificityIndex[order(-KMer_SpecificityIndex$COUNT), ]
    row.names(ordered_SI) = NULL
    UUUUU_Rank = c(UUUUU_Rank, as.numeric(rownames(ordered_SI[ordered_SI$MOTIF=='TTTTT', ])))
    GGGGG_Rank = c(GGGGG_Rank, as.numeric(rownames(ordered_SI[ordered_SI$MOTIF=='GGGGG', ])))
    
    hist_count = hist(KMer_Matches_Scaled$COUNT, breaks = seq(0, 1, 0.01), title = 'test')
    hist_count = data.frame(CNTs = hist_count$counts)
    
    print(paste0('     '))
    print(paste0('     Writing Kmer information....'))
    
    write.csv(Peak_KMers_C,
              paste0(baseDir, '/Output/', DATE, '/', sample, '_RelativeFrequency_', extension, 'ntExt_5mer.csv'), quote = F)
    write.csv(KMer_Matches_Scaled,
              paste0(baseDir, '/Output/', DATE, '/', sample, '_RelativeFrequencyNormalized_', extension, 'ntExt_5mer.csv'), quote = F)
    write.csv(KMer_SpecificityIndex,
              paste0(baseDir, '/Output/', DATE, '/', sample, '_SIperMotif_', extension, 'ntExt_5mer.csv'), quote = F)
    write.csv(hist_count,
              paste0(baseDir, '/Output/', DATE, '/', sample, '_Histogram_', extension, 'ntExt_5mer.csv'), quote = F)
    
  }
}

TopMotif_Data_Compiled = data.frame(CLIP = UUUUU_CLIP,
                                 EXTN = UUUUU_Ext,
                                 UUUUU_RANK = UUUUU_Rank,
                                 UUUUU_SIDX = UUUUU_SI,
                                 GGGGG_RANK = GGGGG_Rank,
                                 GGGGG_SIDX = GGGGG_SI
                                 )

write.csv(TopMotif_Data_Compiled,
          paste0(baseDir, '/Output/', DATE, '/TopMotif_SI_5mer_Optimization.csv'), quote = F)

################################################################################

## Stacked Bar
################################################################################
## Get Annotation counts:
countAnnotation = function(peak_matrix, annotation_column, new_column_name = NULL, annotation_to_skip = NULL, fraction = NULL) {
  temp = data.frame(table(peak_matrix[, ..annotation_column]), row.names = 1)
  if(!is.null(new_column_name)) {
    colnames(temp) = new_column_name
  }
  
  if(!is.null(annotation_to_skip)) {
    temp = temp[rownames(temp) != annotation_to_skip, , drop = FALSE]
  }
  
  if(!is.null(fraction)) {
    temp = temp/sum(temp)
  }
  
  return(temp)
}

## Fill Annotation counts if anything is mixing:
fillAnnotation = function(annotation_counts, annotation_list) {
  colnames((annotation_counts))
  
  temp = data.frame(Sample = numeric(length(annotation_list)))
  rownames(temp) = annotation_list
  temp2 = merge(temp, annotation_counts, by = "row.names", all = TRUE)
  temp2[is.na(temp2)] = 0 
  temp2$Sample = NULL
  
  rownames(temp2) = temp2$Row.names
  temp2 = temp2[ncRNA_List, -1, drop = FALSE]
  
  return(temp2)
}

## Plot Stacked bar:
plotStackedBar = function(annotation_counts, sample_list, sample_label, title, y_lim = NULL, y_tick = NULL) {
  plot = ggplot(annotation_counts %>% filter(Source %in% sample_list), aes(fill = Annotation, y=Freq, x=Source)) + 
    geom_bar(position='stack', stat='identity') +
    scale_x_discrete(labels = sample_label) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette = "Set3") +
    theme_bw() + 
    theme(axis.text = element_text(size=14), 
          axis.title = element_text(size=14, face = 'bold'), 
          legend.text = element_text(size=14))
  
  if(!is.null(y_lim)) {
    plot = plot + ylim(y_lim) 
  }
  
  if (!is.null(y_tick)) {
    plot = plot + scale_y_continuous(breaks = seq(0, y_lim[2], by=y_tick), limits=c(0, y_lim[2]))
  }
  
  return(plot)
}


CLIP_List = c('eCLIP', 'HeLa', 'hnRNPC_WT', 'hnRNPC_Mut', 'RBM25_WT', 'RBM25_Mut')
All_Annotation_List = c("5'UTR", "CDS", "3'UTR", "intron", 'rRNA', 'tRNA', 'ncRNA', 'TE', "Other")

Peak_eCLIP = Peak_eCLIP %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_eCLIP = Peak_eCLIP %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_HeLa = Peak_HeLa %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_HeLa = Peak_HeLa %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_RBM25_WT = Peak_RBM25_WT %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_RBM25_WT = Peak_RBM25_WT %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_RBM25_Mut = Peak_RBM25_Mut %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_RBM25_Mut = Peak_RBM25_Mut %>% filter(finalized_annotation %in% All_Annotation_List)

Counts_Peak_eCLIP = countAnnotation(Peak_eCLIP, 'finalized_annotation', 'eCLIP', fraction = TRUE)
Counts_Peak_eCLIP['tRNA', 'eCLIP'] = 0
Counts_Peak_HeLa = countAnnotation(Peak_HeLa, 'finalized_annotation', 'HeLa', fraction = TRUE)
Counts_Peak_hnRNPC_WT = countAnnotation(Peak_hnRNPC_WT, 'finalized_annotation', 'hnRNPC_WT', fraction = TRUE)
Counts_Peak_hnRNPC_Mut = countAnnotation(Peak_hnRNPC_Mut, 'finalized_annotation', 'hnRNPC_Mut', fraction = TRUE)
Counts_Peak_RBM25_WT = countAnnotation(Peak_RBM25_WT, 'finalized_annotation', 'RBM25_WT', fraction = TRUE)
Counts_Peak_RBM25_Mut = countAnnotation(Peak_RBM25_Mut, 'finalized_annotation', 'RBM25_Mut', fraction = TRUE)

PeakDistribution_combined = cbind(Counts_Peak_eCLIP, Counts_Peak_HeLa, 
                                  Counts_Peak_hnRNPC_WT, Counts_Peak_hnRNPC_Mut, 
                                  Counts_Peak_RBM25_WT, Counts_Peak_RBM25_Mut)
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>% gather(key = "Source", value = "Freq", CLIP_List) %>% dplyr::select(Source, Freq, Annotation)
PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = CLIP_List)
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = All_Annotation_List)

plotStackedBar(PeakDistribution_combined, CLIP_List, CLIP_List, 'CLIP')
################################################################################

## PWM
################################################################################
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(Biostrings)
library(seqLogo)

baseDir = '~/Desktop/Genomics/Specificity/Sequencing/Processed/IndividualPeaks'

numMotif = 10

## For hnRNPC WT
# data = read_csv(paste0(baseDir, '/Output/All_Peaks/For_PWM/hnRNPC_WT.csv'), col_names = T, show_col_types = F)
data = read_csv(paste0(baseDir, '/Output/', DATE, '/For_PWM/hnRNPC_WT.csv'), col_names = T, show_col_types = F)

data = data.frame(data[, c('MOTIF', 'COUNT')])
data = data[order(-data[, 'COUNT']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## For hnRNPC Mut
# data = read_csv(paste0(baseDir, '/Output/All_Peaks/For_PWM/hnRNPC_Mut.csv'), col_names = T, show_col_types = F)
data = read_csv(paste0(baseDir, '/Output/', DATE, '/For_PWM/hnRNPC_Mut.csv'), col_names = T, show_col_types = F)

data = data.frame(data[, c('MOTIF', 'COUNT')])
data = data[order(-data[, 'COUNT']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## For RBM25 WT
# data = read_csv(paste0(baseDir, '/Output/All_Peaks/For_PWM/RBM25_WT.csv'), col_names = T, show_col_types = F)
data = read_csv(paste0(baseDir, '/Output/', DATE, '/For_PWM/RBM25_WT.csv'), col_names = T, show_col_types = F)

data = data.frame(data[, c('MOTIF', 'COUNT')])
data = data[order(-data[, 'COUNT']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## For RBM25 Mut
# data = read_csv(paste0(baseDir, '/Output/All_Peaks/For_PWM/RBM25_Mut.csv'), col_names = T, show_col_types = F)
data = read_csv(paste0(baseDir, '/Output/', DATE, '/For_PWM/RBM25_Mut.csv'), col_names = T, show_col_types = F)

data = data.frame(data[, c('MOTIF', 'COUNT')])
data = data[order(-data[, 'COUNT']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
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

## Histogram
################################################################################
# hnRNPC_WT = read_csv(paste0(baseDir, '/Output/All_Peaks/For_PWM/hnRNPC_WT.csv'), col_names = T, show_col_types = F)
# hnRNPC_Mut = read_csv(paste0(baseDir, '/Output/All_Peaks/For_PWM/hnRNPC_Mut.csv'), col_names = T, show_col_types = F)
hnRNPC_WT = read_csv(paste0(baseDir, '/Output/', DATE, '/For_PWM/hnRNPC_WT.csv'), col_names = T, show_col_types = F)
hnRNPC_Mut = read_csv(paste0(baseDir, '/Output/', DATE, '/For_PWM/hnRNPC_Mut.csv'), col_names = T, show_col_types = F)

ggplot() +
  geom_histogram(data = hnRNPC_WT, aes(x = COUNT), binwidth = 0.01, fill = "darkgoldenrod1", alpha = 0.5) +
  geom_histogram(data = hnRNPC_Mut, aes(x = COUNT), binwidth = 0.01, fill = "cornflowerblue", alpha = 0.5) +
  labs(title = "hnRNPC", x = "Normalized Counts", y = "Frequency") +
  scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, by = 50)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))


# RBM25_WT = read_csv(paste0(baseDir, '/Output/All_Peaks/For_PWM/RBM25_WT.csv'), col_names = T, show_col_types = F)
# RBM25_Mut = read_csv(paste0(baseDir, '/Output/All_Peaks/For_PWM/RBM25_Mut.csv'), col_names = T, show_col_types = F)
RBM25_WT = read_csv(paste0(baseDir, '/Output/', DATE, '/For_PWM/RBM25_WT.csv'), col_names = T, show_col_types = F)
RBM25_Mut = read_csv(paste0(baseDir, '/Output/', DATE, '/For_PWM/RBM25_Mut.csv'), col_names = T, show_col_types = F)

ggplot() +
  geom_histogram(data = RBM25_WT, aes(x = COUNT), binwidth = 0.01, fill = "darkgoldenrod1", alpha = 0.5) +
  geom_histogram(data = RBM25_Mut, aes(x = COUNT), binwidth = 0.01, fill = "cornflowerblue", alpha = 0.5) +
  labs(title = "RBM25", x = "Normalized Counts", y = "Frequency") +
  scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, by = 50)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))

################################################################################

## MS Analysis
################################################################################
# hnRNPC_WT = read_csv(paste0(baseDir, '/Output/All_Peaks/For_PWM/hnRNPC_WT.csv'), col_names = T, show_col_types = F)
hnRNPC_WT = data.frame(hnRNPC_WT[, c('MOTIF', 'COUNT')])
colnames(hnRNPC_WT) = c('Motif', 'Score')
hnRNPC_WT = hnRNPC_WT %>% mutate(Motif = str_replace_all(Motif, "T", "U"))

# hnRNPC_Mut = read_csv(paste0(baseDir, '/Output/All_Peaks/For_PWM/hnRNPC_Mut.csv'), col_names = T, show_col_types = F)
hnRNPC_Mut = data.frame(hnRNPC_Mut[, c('MOTIF', 'COUNT')])
colnames(hnRNPC_Mut) = c('Motif', 'Score')
hnRNPC_Mut = hnRNPC_Mut %>% mutate(Motif = str_replace_all(Motif, "T", "U"))

# RBM25_WT = read_csv(paste0(baseDir, '/Output/All_Peaks/For_PWM/RBM25_WT.csv'), col_names = T, show_col_types = F)
RBM25_WT = data.frame(RBM25_WT[, c('MOTIF', 'COUNT')])
colnames(RBM25_WT) = c('Motif', 'Score')
RBM25_WT = RBM25_WT %>% mutate(Motif = str_replace_all(Motif, "T", "U"))

# RBM25_Mut = read_csv(paste0(baseDir, '/Output/All_Peaks/For_PWM/RBM25_Mut.csv'), col_names = T, show_col_types = F)
RBM25_Mut = data.frame(RBM25_Mut[, c('MOTIF', 'COUNT')])
colnames(RBM25_Mut) = c('Motif', 'Score')
RBM25_Mut = RBM25_Mut %>% mutate(Motif = str_replace_all(Motif, "T", "U"))

plot_MS(return_MS(hnRNPC_WT))
plot_MS(return_MS(hnRNPC_Mut))
plot_MS(return_MS(RBM25_WT))
plot_MS(return_MS(RBM25_Mut))

RBM25_WT_UUUUUAnchored = RBM25_WT
RBM25_WT_UUUUUAnchored$Score = RBM25_WT$Score + (1 - RBM25_WT$Score[RBM25_WT$Motif == "UUUUU"])

RBM25_Mut_UUUUUAnchored = RBM25_Mut
RBM25_Mut_UUUUUAnchored$Score = RBM25_Mut$Score + (1 - RBM25_Mut$Score[RBM25_Mut$Motif == "UUUUU"])

# plot_MS(return_MS(RBM25_WT_UUUUUAnchored, top_motif = "UUUUU"))
# plot_MS(return_MS(RBM25_Mut_UUUUUAnchored, top_motif = "UUUUU"))

plot_MS(return_MS(RBM25_WT, top_motif = "UUUUU"))
plot_MS(return_MS(RBM25_Mut, top_motif = "UUUUU"))


MS_hnRNPC_WT = return_MS(hnRNPC_WT)
MS_hnRNPC_WT = MS_hnRNPC_WT[, c("Pos1", "Pos2", "Pos3", "Pos4", "Pos5")] - min(MS_hnRNPC_WT[, c("Pos1", "Pos2", "Pos3", "Pos4", "Pos5")])
mean(unlist(MS_hnRNPC_WT))

MS_hnRNPC_Mut = return_MS(hnRNPC_Mut)
MS_hnRNPC_Mut = MS_hnRNPC_Mut[, c("Pos1", "Pos2", "Pos3", "Pos4", "Pos5")] - min(MS_hnRNPC_Mut[, c("Pos1", "Pos2", "Pos3", "Pos4", "Pos5")])
mean(unlist(MS_hnRNPC_Mut))

MS_RBM25_WT = return_MS(RBM25_WT, "UUUUU")
MS_RBM25_WT = MS_RBM25_WT[, c("Pos1", "Pos2", "Pos3", "Pos4", "Pos5")] - min(MS_RBM25_WT[, c("Pos1", "Pos2", "Pos3", "Pos4", "Pos5")])
mean(unlist(MS_RBM25_WT))

MS_RBM25_Mut = return_MS(RBM25_Mut, "UUUUU")
MS_RBM25_Mut = MS_RBM25_Mut[, c("Pos1", "Pos2", "Pos3", "Pos4", "Pos5")] - min(MS_RBM25_Mut[, c("Pos1", "Pos2", "Pos3", "Pos4", "Pos5")])
mean(unlist(MS_RBM25_Mut))

hnRNPC_FC = data.frame(Motif = hnRNPC_WT$Motif,
                       WT = hnRNPC_WT$Score,
                       Mut = hnRNPC_Mut$Score,
                       FC = hnRNPC_Mut$Score/hnRNPC_WT$Score)

RBM25_FC = data.frame(Motif = RBM25_WT$Motif,
                      WT = RBM25_WT$Score,
                      Mut = RBM25_Mut$Score,
                      FC = RBM25_Mut$Score/RBM25_WT$Score)

################################################################################

## Peak Scatter Plot
################################################################################

peaks2GR = function(peaks_df) {
  peaks_df = peaks_df[, c('chr', 'start', 'end', 'strand', 'name')]
  peaks_df$start = as.integer(peaks_df$start)
  peaks_df$end = as.integer(peaks_df$end)
  peaks_GR = GRanges(peaks_df)
  
  return(peaks_GR)
}

GR_Peak_eCLIP = peaks2GR(Peak_eCLIP)
GR_Peak_HeLa = peaks2GR(Peak_HeLa)
GR_Peak_hnRNPC_WT = peaks2GR(Peak_hnRNPC_WT)
GR_Peak_hnRNPC_Mut = peaks2GR(Peak_hnRNPC_Mut)
GR_Peak_RBM25_WT = peaks2GR(Peak_RBM25_WT)
GR_Peak_RBM25_Mut = peaks2GR(Peak_RBM25_Mut)

# hnRNPC_WT vs hnRNPC_Mut
overlaps = as.data.frame(findOverlaps(GR_Peak_hnRNPC_WT, GR_Peak_hnRNPC_Mut, minoverlap = 15))

overlaps_hnRNPC_WT = Peak_hnRNPC_WT[overlaps$queryHits]
overlaps_hnRNPC_Mut = Peak_hnRNPC_Mut[overlaps$subjectHits]

overlaps_RPM = data.frame(hnRNPC_WT = overlaps_hnRNPC_WT$hnRNPC_WT_avg_RPM,
                          hnRNPC_Mut = overlaps_hnRNPC_Mut$hnRNPC_Mut_avg_RPM)

RMSDist = sqrt(mean((overlaps_RPM$hnRNPC_Mut - overlaps_RPM$hnRNPC_WT)^2)) # for magnitude
MDist = mean(overlaps_RPM$hnRNPC_WT - overlaps_RPM$hnRNPC_Mut) # for direction
L2FC = log2(mean(overlaps_RPM$hnRNPC_Mut/overlaps_RPM$hnRNPC_WT))

# Calculate the *perpendicular* distance from the y=x line (as before)
custom_color_function = function(dist, pos) {
  if (dist < 0.001) { # Threshold for "near" (adjust as needed)
    return("black")
  } else {
    if (pos == "Above") {
      # Use pmax to avoid 0 index
      index <- pmax(1, pmin(100, round(dist * 100))) # Correct indexing
      return(colorRampPalette(c("black", "#ed6677"))(100)[index])
    } else {
      index <- pmax(1, pmin(100, round(dist * 100))) # Correct indexing
      return(colorRampPalette(c("black", "#258942"))(100)[index])
    }
  }
}

overlaps_RPM$distance = ((abs(overlaps_RPM$hnRNPC_Mut - overlaps_RPM$hnRNPC_WT) / sqrt(2)))
overlaps_RPM$position = ifelse(overlaps_RPM$hnRNPC_Mut > overlaps_RPM$hnRNPC_WT, "Above", "Below")
overlaps_RPM$color = mapply(custom_color_function, overlaps_RPM$distance, overlaps_RPM$position)

ggplot(overlaps_RPM, 
       aes(x = hnRNPC_WT, y = hnRNPC_Mut, color = color)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_color_identity() +  # Use the pre-calculated colors
  scale_x_log10(limits = c(1, 1000)) +
  scale_y_log10(limits = c(1, 1000)) +
  labs(title = paste0("hnRNPC WT vs Mut RPM Comparison (N: ", (nrow(overlaps_RPM)), ", RMSD: ", (RMSDist), ", MD: ", (MDist), ")"),
       x = "hnRNPC WT Average RPM",
       y = "hnRNPC Mut Average RPM") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))

overlaps = as.data.frame(findOverlaps(GR_Peak_RBM25_WT, GR_Peak_RBM25_Mut))
overlaps_RBM25_WT = Peak_RBM25_WT[overlaps$queryHits]
overlaps_RBM25_Mut = Peak_RBM25_Mut[overlaps$subjectHits]

overlaps_RPM = data.frame(RBM25_WT = overlaps_RBM25_WT$RBM25_WT_avg_RPM,
                          RBM25_Mut = overlaps_RBM25_Mut$RBM25_Mut_avg_RPM)

RMSDist = sqrt(mean((overlaps_RPM$RBM25_Mut - overlaps_RPM$RBM25_WT)^2)) # for magnitude
MDist = mean(overlaps_RPM$RBM25_WT - overlaps_RPM$RBM25_Mut) # for direction
L2FC = log2(mean(overlaps_RPM$RBM25_Mut/overlaps_RPM$RBM25_WT))

overlaps_RPM$distance = (abs(overlaps_RPM$RBM25_Mut - overlaps_RPM$RBM25_WT) / sqrt(2))
overlaps_RPM$position = ifelse(overlaps_RPM$RBM25_Mut > overlaps_RPM$RBM25_WT, "Above", "Below")
overlaps_RPM$color = mapply(custom_color_function, overlaps_RPM$distance, overlaps_RPM$position)

ggplot(overlaps_RPM, 
       aes(x = RBM25_WT, y = RBM25_Mut, color = color)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_color_identity() +  # Use the pre-calculated colors
  scale_x_log10(limits = c(1, 1000)) +
  scale_y_log10(limits = c(1, 1000)) +
  labs(title = paste0("RBM25 WT vs Mut RPM Comparison (N: ", (nrow(overlaps_RPM)), ", RMSD: ", (RMSDist), ", MD: ", (MDist), ")"),
       x = "RBM25 WT Average RPM",
       y = "RBM25 Mut Average RPM") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))

################################################################################

## SI Analysis
################################################################################
hnRNPC_WT$SI = hnRNPC_WT$COUNT/median(hnRNPC_WT$COUNT)
hnRNPC_Mut$SI = hnRNPC_Mut$COUNT/median(hnRNPC_Mut$COUNT)
hnRNPC_Mut2WT_SIRatio = data.frame(MOTIF = hnRNPC_WT$MOTIF,
                                  WT_SI = hnRNPC_WT$SI,
                                  Mut_SI = hnRNPC_Mut$SI,
                                  SI_Raio = (hnRNPC_Mut$SI+0.001)/(hnRNPC_WT$SI+0.001))

RBM25_WT$SI = RBM25_WT$COUNT/median(RBM25_WT$COUNT)
RBM25_Mut$SI = RBM25_Mut$COUNT/median(RBM25_Mut$COUNT)
RBM25_Mut2WT_SIRatio = data.frame(MOTIF = RBM25_WT$MOTIF,
                                  WT_SI = RBM25_WT$SI,
                                  Mut_SI = RBM25_Mut$SI,
                                  SI_Raio = (RBM25_Mut$SI+0.001)/(RBM25_WT$SI+0.001))

numMotif = 20

## For hnRNPC, top 20 SI changes
data = data.frame(hnRNPC_Mut2WT_SIRatio[, c('MOTIF', 'SI_Raio')])
data = data[order(-data[, 'SI_Raio']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## For hnRNPC, bottom 20 SI changes
data = data[order(data[, 'SI_Raio']), ]

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## For RBM25, top 20 SI changes
data = data.frame(RBM25_Mut2WT_SIRatio[, c('MOTIF', 'SI_Raio')])
data = data[order(-data[, 'SI_Raio']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## For RBM25, bottom 20 SI changes
data = data[order(data[, 'SI_Raio']), ]

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

K <- 5
NTs <- c("A", "C", "G", "T")

# Generate all 5-mers
KMer_polyT = expand.grid(replicate(K, NTs, simplify = FALSE))
KMer_polyT = apply(KMer_polyT, 1, paste, collapse = "")
KMer_polyT = KMer_polyT[grepl("TTT", KMer_polyT)]

RBM25_polyT = RBM25_Mut2WT_SIRatio %>% filter(MOTIF %in% KMer_polyT) %>% filter(MOTIF != "TTTTT")

################################################################################

