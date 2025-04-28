################################################################################
## eCLIP Analysis V5: 
## Written by Soon Yi
## Created: 2023-12-08
## Last Edited: 2025-04-26
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

return_MS_forAll = function(data) {
  
  unique_motifs_t = unique(data$MOTIF)
  
  results_list = vector("list", length(unique_motifs_t))
  names(results_list) = unique_motifs_t
  
  # --- Loop through each unique original motif (with T) ---
  for (motif_t in unique_motifs_t) {
    
    # 1. Generate true variants using the USER'S Motif_Variants function
    #    Using 'T' nucleotides. Applying unique() here for safety.
    variants_t = unique(Motif_Variants(motif_t, nucleotides = c('A', 'C', 'G', 'T')))
    
    if (length(variants_t) == 0) {
      results_list[[motif_t]] = NA_real_
      next
    }
    
    # 2. Find scores for these variants in the input data
    match_indices = match(variants_t, data$MOTIF)
    valid_indices = match_indices[!is.na(match_indices)]
    variant_scores = data$COUNT[valid_indices]
    
    # 3. Calculate 1 - score for each valid, non-NA score
    ms_scores_for_variants = 1 - variant_scores[!is.na(variant_scores)]
    
    # 4. Calculate average MS
    scalar_ms = NA_real_
    if (length(ms_scores_for_variants) > 0) {
      if(all(is.na(ms_scores_for_variants))) {
        scalar_ms = NA_real_
      } else {
        mean_val = mean(ms_scores_for_variants, na.rm = TRUE)
        scalar_ms = ifelse(is.nan(mean_val), NA_real_, mean_val)
      }
    }
    
    results_list[[motif_t]] = scalar_ms # Store result
  } # --- End loop ---
  
  # Create the final output data frame
  final_df = data.frame(
    MOTIF = names(results_list),
    MS = unlist(results_list),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  return(final_df)
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

DATE = format(Sys.Date(), "%Y%m%d")

extensions = c(0, 5, 10, 25, 50, 100)

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
    for (extension in extensions) {
      # Read peaks and make GR object
      print(paste0('Working on ', RBP, '....'))
      print(paste0('     Reading and processing peak data....'))
      file_ID = sampleTable_eCLIP[sampleTable_eCLIP$RBP == RBP & sampleTable_eCLIP$Cell == cell_line, ]$eCLIP
      Peak = fread(paste0(baseDir, 'eCLIP_peak/raw/', file_ID,  '.bed'))
      Peak = Peak %>% mutate(V1 = case_when(V1 == 'chrMT' ~ 'chrM', TRUE ~ as.character(V1)))
      Peak = Peak %>% filter(V1 %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 
                                       'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
      
      #######
      # Peak_GR = GRanges(seqnames = Rle(Peak$V1),
      #                   ranges = IRanges(start = Peak$V2 - extension, end = Peak$V3),
      #                   strand = Rle(Peak$V6))
      
      Peak_GR_original = GRanges(seqnames = Rle(Peak$V1),
                                 ranges = IRanges(start = Peak$V2, end = Peak$V3),
                                 strand = Rle(Peak$V6))
      original_width = width(Peak_GR_original)
      new_width = original_width + extension
      Peak_GR = resize(Peak_GR_original, width = new_width, fix = "end")
      #######
      
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
      
      # Set the 100th average as the background K-mer counts:
      Bkg_KMer$COUNT = AVG[, length(colnames(AVG))]
      
      # First, calculate the background corrected K-mer counts:
      Peak_KMers_C = Peak_KMers
      Peak_KMers_C$COUNT = Peak_KMers_C$COUNT - Bkg_KMer$COUNT
      # Peak_KMers_C$COUNT = Peak_KMers_C$COUNT / Bkg_KMer$COUNT
      
      # Second, scale the normalized counts to scale of 0-to-1:
      KMer_Matches_Scaled = Peak_KMers_C
      # KMer_Matches_Scaled = KMer_Matches_Scaled %>% mutate_if(is.numeric, min_max_norm) %>% mutate_if(is.numeric, round, 10)
      KMer_Matches_Scaled = KMer_Matches_Scaled %>% mutate_if(is.numeric, min_max_norm, a = 1, b = exp(1)) %>% mutate_if(is.numeric, round, 10)
      KMer_Matches_Scaled$COUNT = log(KMer_Matches_Scaled$COUNT)
      
      # Third, calculate the specificity index for each K-mer:
      KMer_SpecificityIndex = KMer_Matches_Scaled
      KMer_SpecificityIndex = KMer_Matches_Scaled %>% mutate_if(is.numeric, SI_Calculation) %>% mutate_if(is.numeric, round, 10)
      
      # Fourth, calculate MS for each K-mer:
      KMer_MutationalSensitivity = return_MS_forAll(KMer_Matches_Scaled)
      
      print(paste0('     Writing all Kmer information....'))
      
      Bkg_Output_Dir = paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', extension, 'nt_extension/', cell_line, '/0_Bkg/')
      Raw_Output_Dir = paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', extension, 'nt_extension/', cell_line, '/1_Raw/')
      Scaled_Output_Dir = paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', extension, 'nt_extension/', cell_line, '/2_Scaled/')
      SI_Output_Dir = paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', extension, 'nt_extension/', cell_line, '/3_SpecificityIndex/')
      MS_Output_Dir = paste0(baseDir, 'eCLIP_peak/Output/' , DATE, '/', extension, 'nt_extension/', cell_line, '/4_MutationalSensitivity/')
      
      
      if (!dir.exists(Bkg_Output_Dir)) {
        dir.create(Bkg_Output_Dir, recursive = TRUE, showWarnings = FALSE)
        print(paste("Created directory:", Bkg_Output_Dir))
      }
      if (!dir.exists(Raw_Output_Dir)) {
        dir.create(Raw_Output_Dir, recursive = TRUE, showWarnings = FALSE)
        print(paste("Created directory:", Raw_Output_Dir))
      }
      if (!dir.exists(Scaled_Output_Dir)) {
        dir.create(Scaled_Output_Dir, recursive = TRUE, showWarnings = FALSE)
        print(paste("Created directory:", Scaled_Output_Dir))
      }
      if (!dir.exists(SI_Output_Dir)) {
        dir.create(SI_Output_Dir, recursive = TRUE, showWarnings = FALSE)
        print(paste("Created directory:", SI_Output_Dir))
      }
      if (!dir.exists(MS_Output_Dir)) {
        dir.create(MS_Output_Dir, recursive = TRUE, showWarnings = FALSE)
        print(paste("Created directory:", MS_Output_Dir))
      }
      
      write.csv(AVG,
                paste0(Bkg_Output_Dir, RBP, '_', K, 'mer.csv'), quote = F)
      write.csv(Peak_KMers_C, 
                paste0(Raw_Output_Dir, RBP, '_', K, 'mer.csv'), quote = F)
      write.csv(KMer_Matches_Scaled, 
                paste0(Scaled_Output_Dir, RBP, '_', K, 'mer.csv'), quote = F)
      write.csv(KMer_SpecificityIndex, 
                paste0(SI_Output_Dir, RBP, '_', K, 'mer.csv'), quote = F)
      write.csv(KMer_MutationalSensitivity, 
                paste0(MS_Output_Dir, RBP, '_', K, 'mer.csv'), quote = F)
    }
    }
}
################################################################################
















