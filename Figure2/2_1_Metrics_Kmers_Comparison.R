################################################################################
## Compare IS and MS between K-mers 
## Written by Soon Yi
## Created: 2023-11-11
## Last Edited: 2023-11-12
## Figure 2
################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

# Set some basic parameters and custom functions:
################################################################################
baseDir = '~/Desktop/Genomics/Specificity/Data_Final/'

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


##


# Read RBNS file:
################################################################################
RBNS_4mer = read_csv(paste0(baseDir, 'RBNS/RBNS_singleRBD_normalized_4mer.csv'), col_names = T, show_col_types = F)
RBNS_5mer = read_csv(paste0(baseDir, 'RBNS/RBNS_singleRBD_normalized_5mer.csv'), col_names = T, show_col_types = F)
RBNS_6mer = read_csv(paste0(baseDir, 'RBNS/RBNS_singleRBD_normalized_6mer.csv'), col_names = T, show_col_types = F)
RBNS_7mer = read_csv(paste0(baseDir, 'RBNS/RBNS_singleRBD_normalized_7mer.csv'), col_names = T, show_col_types = F)

RBPs_4mer = colnames(RBNS_4mer)[2:ncol(RBNS_4mer)]
RBPs_5mer = colnames(RBNS_5mer)[2:ncol(RBNS_5mer)]
RBPs_6mer = colnames(RBNS_6mer)[2:ncol(RBNS_6mer)]
RBPs_7mer = colnames(RBNS_7mer)[2:ncol(RBNS_7mer)]
################################################################################

# Calculate Specificity Indexes for RBPs:
################################################################################
SpecificityIndex = as.data.frame(matrix(NA, nrow = 26, ncol = 5))
colnames(SpecificityIndex) = c('RBP', '4mer', '5mer', '6mer', '7mer')
SpecificityIndex$RBP = RBPs_5mer

Ks = c(4, 5, 6, 7)

for (K in Ks) {
  RBNS_Data = get(paste0('RBNS_', K, 'mer'))
  RBP_Data = get(paste0('RBPs_', K, 'mer'))
  for (RBP in RBP_Data) {
    SpecificityIndex[which(SpecificityIndex$RBP == RBP), paste0(K, 'mer')] = max(unlist(RBNS_Data[, RBP])) / median(unlist(RBNS_Data[, RBP]))
  }
}
################################################################################

# Plot Specificity Index across different K-mers:
################################################################################
long_data = SpecificityIndex %>%
  pivot_longer(cols = c('4mer', '5mer', '6mer', '7mer'), 
               names_to = 'mer', 
               values_to = 'value')

ggplot(long_data, aes(x = mer, y = value, group = RBP, color = RBP)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "K-Mer", y = "Specificity Index", title = "SI per RBP across Different K-mers") +
  scale_y_continuous(limits = c(0, 275), 
                     breaks = seq(0, 300, by = 25))  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
################################################################################

# Plot Specificity Index rank across different K-mers:
################################################################################
long_data_rank = SpecificityIndex %>%
  gather(key = "mer", value = "value", -RBP) %>%
  group_by(mer) %>%
  mutate(rank = rank(-value)) %>%
  ungroup()

long_data_rank[which(is.na(long_data_rank$value)), ]$rank = NA

# In Heatmap Form:
ranked_df = long_data_rank %>%
  select(RBP, mer, rank) %>%
  pivot_wider(names_from = mer, values_from = rank, values_fill = NA) %>%
  arrange(RBP)

heatmap_data = ranked_df %>% column_to_rownames(var = 'RBP')
avg_ranks = rowMeans(heatmap_data, na.rm = TRUE)
ordered_heatmap_data = heatmap_data[order(avg_ranks), ]

SI_data_for_heatmap = SpecificityIndex %>% column_to_rownames(var = 'RBP') %>% as.matrix()
SI_data_for_heatmap = round(SI_data_for_heatmap, 1)
SI_data_for_heatmap = SI_data_for_heatmap[rownames(ordered_heatmap_data), ]

pheatmap(ordered_heatmap_data, 
         display_numbers = SI_data_for_heatmap,
         cluster_rows = FALSE,
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)),
         na_col = "grey",
         number_color = 'black')

################################################################################

# Specificity Index Comparison Between Kmers
################################################################################
## 4 v 5
r = cor(SpecificityIndex$'4mer', SpecificityIndex$'5mer', use = 'complete.obs')
ggplot(SpecificityIndex, aes(x = `4mer`, y = `5mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "4mer Specificity Index", y = "5mer Specificity Index", title = "SI Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 100, by = 20)) +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 100, by = 20)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## 4 v 6
r = cor(SpecificityIndex$'4mer', SpecificityIndex$'6mer', use = 'complete.obs')
ggplot(SpecificityIndex, aes(x = `4mer`, y = `6mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "4mer Specificity Index", y = "6mer Specificity Index", title = "SI Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 125), breaks = seq(0, 125, by = 25)) +
  scale_x_continuous(limits = c(0, 125), breaks = seq(0, 125, by = 25)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## 4 v 7
r = cor(SpecificityIndex$'4mer', SpecificityIndex$'7mer', use = 'complete.obs')
ggplot(SpecificityIndex, aes(x = `4mer`, y = `7mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "4mer Specificity Index", y = "7mer Specificity Index", title = "SI Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 275), breaks = seq(0, 300, by = 25)) +
  scale_x_continuous(limits = c(0, 275), breaks = seq(0, 300, by = 25)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## 5 v 6
r = cor(SpecificityIndex$'5mer', SpecificityIndex$'6mer', use = 'complete.obs')
ggplot(SpecificityIndex, aes(x = `5mer`, y = `6mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "5mer Specificity Index", y = "6mer Specificity Index", title = "SI Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 125), breaks = seq(0, 150, by = 25)) +
  scale_x_continuous(limits = c(0, 125), breaks = seq(0, 150, by = 25)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## 5 v 7
r = cor(SpecificityIndex$'5mer', SpecificityIndex$'7mer', use = 'complete.obs')
ggplot(SpecificityIndex, aes(x = `5mer`, y = `7mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "5mer Specificity Index", y = "7mer Specificity Index", title = "SI Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 275), breaks = seq(0, 275, by = 25)) +
  scale_x_continuous(limits = c(0, 275), breaks = seq(0, 275, by = 25)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## 6 v 7
r = cor(SpecificityIndex$'6mer', SpecificityIndex$'7mer', use = 'complete.obs')
ggplot(SpecificityIndex, aes(x = `6mer`, y = `7mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "6mer Specificity Index", y = "7mer Specificity Index", title = "SI Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 275), breaks = seq(0, 275, by = 25)) +
  scale_x_continuous(limits = c(0, 275), breaks = seq(0, 275, by = 25)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
################################################################################


##


# Calculate Mutational Sensitivity for RBPs:
################################################################################
MutationalSensitivity = as.data.frame(matrix(NA, nrow = 26, ncol = 5))
colnames(MutationalSensitivity) = c('RBP', '4mer', '5mer', '6mer', '7mer')
MutationalSensitivity$RBP = RBPs_5mer

Ks = c(4, 5, 6, 7)
NTs = c('A', 'U', 'C', 'G')

for (K in Ks) {
  RBNS_Data = get(paste0('RBNS_', K, 'mer'))
  RBP_Data = get(paste0('RBPs_', K, 'mer'))
  for (RBP in RBP_Data) {
    Motif_Top = RBNS_Data$Motif[which(RBNS_Data[, RBP] == max(RBNS_Data[, RBP]))]
    Motif_Vars = Motif_Variants(Motif_Top, NTs)
    MSs = c()
    for (idx in c(1:length(Motif_Vars))) {
      Motif_Var = Motif_Vars[idx]
      MSs[idx] = 1 - RBNS_Data[which(RBNS_Data$Motif == Motif_Var), ][[RBP]][1]
    }
    MutationalSensitivity[which(MutationalSensitivity$RBP == RBP), paste0(K, 'mer')] = mean(MSs)
  }
}

################################################################################

# Plot Mutational Sensitivity across different K-mers:
################################################################################
long_data = MutationalSensitivity %>%
  pivot_longer(cols = c('4mer', '5mer', '6mer', '7mer'), 
               names_to = 'mer', 
               values_to = 'value')

ggplot(long_data, aes(x = mer, y = value, group = RBP, color = RBP)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "K-Mer", y = "Mutational Sensitivity", title = "MS per RBP across Different K-mers") +
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.1))  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
################################################################################

# Plot Mutational Sensitivity rank across different K-mers:
################################################################################
long_data_rank = MutationalSensitivity %>%
  gather(key = "mer", value = "value", -RBP) %>%
  group_by(mer) %>%
  mutate(rank = rank(-value)) %>%
  ungroup()

long_data_rank[which(is.na(long_data_rank$value)), ]$rank = NA

# In Heatmap Form:
ranked_df = long_data_rank %>%
  select(RBP, mer, rank) %>%
  pivot_wider(names_from = mer, values_from = rank, values_fill = NA) %>%
  arrange(RBP)

heatmap_data = ranked_df %>% column_to_rownames(var = 'RBP')
avg_ranks = rowMeans(heatmap_data, na.rm = TRUE)
ordered_heatmap_data = heatmap_data[order(avg_ranks), ]

MS_data_for_heatmap = MutationalSensitivity %>% column_to_rownames(var = 'RBP') %>% as.matrix()
MS_data_for_heatmap = round(MS_data_for_heatmap, 2)
MS_data_for_heatmap = MS_data_for_heatmap[rownames(ordered_heatmap_data), ]

pheatmap(ordered_heatmap_data, 
         display_numbers = MS_data_for_heatmap,
         cluster_rows = FALSE,
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)),
         na_col = "grey",
         number_color = 'black')

################################################################################

# Mutational Sensitivity Comparison Between Kmers:
################################################################################
## 4 v 5
r = cor(MutationalSensitivity$'4mer', MutationalSensitivity$'5mer', use = 'complete.obs')
ggplot(MutationalSensitivity, aes(x = `4mer`, y = `5mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "4mer Mutational Sensitivity", y = "5mer Mutational Sensitivity", title = "MS Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 1.5, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## 4 v 6
r = cor(MutationalSensitivity$'4mer', MutationalSensitivity$'6mer', use = 'complete.obs')
ggplot(MutationalSensitivity, aes(x = `4mer`, y = `6mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "4mer Mutational Sensitivity", y = "6mer Mutational Sensitivity", title = "MS Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 1.5, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## 4 v 7
r = cor(MutationalSensitivity$'4mer', MutationalSensitivity$'7mer', use = 'complete.obs')
ggplot(MutationalSensitivity, aes(x = `4mer`, y = `7mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "4mer Mutational Sensitivity", y = "7mer Mutational Sensitivity", title = "MS Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 1.5, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## 5 v 6
r = cor(MutationalSensitivity$'5mer', MutationalSensitivity$'6mer', use = 'complete.obs')
ggplot(MutationalSensitivity, aes(x = `5mer`, y = `6mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "5mer Mutational Sensitivity", y = "6mer Mutational Sensitivity", title = "MS Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 1.5, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## 5 v 7
r = cor(MutationalSensitivity$'5mer', MutationalSensitivity$'7mer', use = 'complete.obs')
ggplot(MutationalSensitivity, aes(x = `5mer`, y = `7mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "5mer Mutational Sensitivity", y = "7mer Mutational Sensitivity", title = "MS Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 1.5, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

## 6 v 7
r = cor(MutationalSensitivity$'6mer', MutationalSensitivity$'7mer', use = 'complete.obs')
ggplot(MutationalSensitivity, aes(x = `6mer`, y = `7mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "6mer Mutational Sensitivity", y = "7mer Mutational Sensitivity", title = "MS Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 1.5, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
################################################################################


##


## SI vs MS Between Kmers:
################################################################################
SI_MS_4 = data.frame(RBP = SpecificityIndex$RBP, SI = SpecificityIndex$'4mer', MS = MutationalSensitivity$'4mer')
SI_MS_5 = data.frame(RBP = SpecificityIndex$RBP, SI = SpecificityIndex$'5mer', MS = MutationalSensitivity$'5mer')
SI_MS_6 = data.frame(RBP = SpecificityIndex$RBP, SI = SpecificityIndex$'6mer', MS = MutationalSensitivity$'6mer')
SI_MS_7 = data.frame(RBP = SpecificityIndex$RBP, SI = SpecificityIndex$'7mer', MS = MutationalSensitivity$'7mer')

r = cor(SI_MS_4$SI, SI_MS_4$MS, use = 'complete.obs')
ggplot(SI_MS_4, aes(x = SI, y = MS, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "4mer Specificity Index", y = "4mer Mutational Sensitivity", title = "4-mer: SI vs MS") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 40), breaks = seq(0, 100, by = 5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

r = cor(SI_MS_5$SI, SI_MS_5$MS, use = 'complete.obs')
ggplot(SI_MS_5, aes(x = SI, y = MS, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "5mer Specificity Index", y = "5mer Mutational Sensitivity", title = "5-mer: SI vs MS") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 100, by = 10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

r = cor(SI_MS_6$SI, SI_MS_6$MS, use = 'complete.obs')
ggplot(SI_MS_6, aes(x = SI, y = MS, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "6mer Specificity Index", y = "6mer Mutational Sensitivity", title = "6-mer: SI vs MS") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 160), breaks = seq(0, 160, by = 20)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

r = cor(SI_MS_7$SI, SI_MS_7$MS, use = 'complete.obs')
ggplot(SI_MS_7, aes(x = SI, y = MS, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "7mer Specificity Index", y = "7mer Mutational Sensitivity", title = "7-mer: SI vs MS") +
  annotate("text", x = Inf, y = Inf, label = bquote("R" == .(format(r, digits = 2))), hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 320), breaks = seq(0, 320, by = 40)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
################################################################################


##


## Heatmap of correlations:
################################################################################
SI_Data = SpecificityIndex[, c('4mer', '5mer', '6mer', '7mer')]
colnames(SI_Data) = c('SI_4', 'SI_5', 'SI_6', 'SI_7')
MS_Data = MutationalSensitivity[, c('4mer', '5mer', '6mer', '7mer')]
colnames(MS_Data) = c('MS_4', 'MS_5', 'MS_6', 'MS_7')

# Generate SI heatmap
KMer_Data = SI_Data
CorrMatrix = cor(KMer_Data,  use = "pairwise.complete.obs")
CorrMatrix = matrix(round(CorrMatrix,2), nrow = ncol(KMer_Data))
rownames(CorrMatrix) = colnames(KMer_Data)
colnames(CorrMatrix) = colnames(KMer_Data)

color_palette = (colorRampPalette(brewer.pal(9, "GnBu"))(100))
# color_palette = (colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
breaks = seq(0.5, 1, length.out = length(color_palette) + 1)
pheatmap(CorrMatrix, cluster_rows = F, cluster_cols = F, display_numbers = T, color = color_palette, breaks = breaks)

# Generate MS heatmap
KMer_Data = MS_Data
CorrMatrix = cor(KMer_Data,  use = "pairwise.complete.obs")
CorrMatrix = matrix(round(CorrMatrix,2), nrow = ncol(KMer_Data))
rownames(CorrMatrix) = colnames(KMer_Data)
colnames(CorrMatrix) = colnames(KMer_Data)

color_palette = (colorRampPalette(brewer.pal(9, "GnBu"))(100))
# color_palette = (colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
breaks = seq(0.5, 1, length.out = length(color_palette) + 1)
pheatmap(CorrMatrix, cluster_rows = F, cluster_cols = F, display_numbers = T, color = color_palette, breaks = breaks)
################################################################################


