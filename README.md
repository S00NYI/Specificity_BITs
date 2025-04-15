## RBP Specificity Analysis Scripts:
- Author: Soon Yi
- Last updated: 2025-04-15
-------------------------------------------------------------------------------------------------------------------
![alt text](https://github.com/S00NYI/Specificity_BITs/blob/main/Analysis_Outline.png?raw=true)
-------------------------------------------------------------------------------------------------------------------

#### 0_1_RBNS_DataParsing.R
  - RBNS data parsing based on the RBNS sample list.

#### 0_2_RBNS_DataParsing.R
  - IS and MS calculation across K-mers for RBNS data.

#### 1_1_RBNS_to_Distribution.R
  - Generate normalized distribution plot for selected RBPs.

#### 1_2_RBNS_to_PWM_Logo.R
  - Generate PWM and sequence logo for selected RBPs.

#### 1_3_RBNS_to_PWM_Comparison.R
  - Generate score based on PWM and compare to RBNS.

#### 2_1_Metrics_Kmers_Comparison.R
  - Compare IS/MS values across K-mers for RBPs with single RBD purified.
  
#### 3_1_eCLIP_Analysis.R
  - eCLIP data processing based on the eCLIP sample list.
   
#### 3_2_eCLIP_Analysis_Plot.R
  - eCLIP data analysis (AS/A-MS calculation, compare to IS/MS).

#### 4_1_RRM_Swap_Analysis.R
  - RRM swap CLIP analysis.

#### 5_1_RNA_RBP_Modeling.py
  - RNA-RBP binding simulation using model RBP data.
   
#### 5_2_RNA_RBP_Modeling.R
  - Simulation result processing for figure output.

#### 6_1_in_vitro_iCLIP_modeling.py
  - Modeling the in vitro iCLIP results with RNACompete data.

