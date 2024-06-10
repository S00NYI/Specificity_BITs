################################################################################
## RNA-RBP Modeling 
## Written by Soon Yi
## Created: 2024-05-15
## Last Edited: 2024-06-06
################################################################################


## Load libraries
################################################################################
import os
import pandas as pd
import numpy as np
import random
from scipy import integrate
from sympy import symbols, expand, collect, Poly
import matplotlib.pyplot as plt
import seaborn as sns
################################################################################

## Set basic parameters and functions:
################################################################################
baseDir = '~/in_vitro_CLIP/'

pd.set_option('display.float_format', lambda x: '%.5f' % x)

def featureScale(X, MAX, MIN):
    if (MAX <= MIN):
        raise ValueError("MAX cannot be lower than or equal to MIN.")
    X = np.array(X)  
    X_min = np.nanmin(X)
    X_max = np.nanmax(X)
    return MIN + ((X - X_min) * (MAX - MIN)) / (X_max - X_min)

def processFASTA(PATH, NAME):
    ID = ''
    SEQ = ''
    
    # Read the FASTA file
    with open(os.path.expanduser(PATH + 'target_RNA/' + NAME + '.fa'), 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract chromosome number and calculate length
                header_info = line.split(' ')[1]  # Get the first part of the header
                ID = header_info.split('=')[1].split(':')[0]  # Extract chromosome number
                start, end = map(int, header_info.split('=')[1].split(':')[1].split('-'))  # Extract start and end positions
            else:
                SEQ += line.strip()  # Add sequence
    
    LEN = len(SEQ)

    # Return the processed information
    return [ID, LEN, NAME, SEQ.replace('T', 'U')]

def generateRNA(L, A=0.25, G=0.25, C=0.25, U=0.25):
    
    if A + G + C + U != 1:
        raise ValueError("Input fractions A, G, C, and U must sum up to 1")

    NTs = {
        'A': A,
        'G': G,
        'C': C,
        'U': U
    }

    seq = [NT for NT, percentage in NTs.items() for _ in range(int(L * percentage))]
    random.shuffle(seq)
    seq = ''.join(seq)
    
    return seq

def allBind(N, x0, P0n_vals, Kn_vals):
    # N:        Number of RBP(s).
    # x0:       Initial concentration of RNA.
    # P0n_vals: List of initial protein concentrations.
    # Kn_vals:  List of Kd values for each protein.
    #
    # Here, we solve for equilibrium concentration of RNA, X.
    # X can be utilized to calculate equilibrium concentrations of RBPs. 

    import sympy as sp

    if (len(P0n_vals) != N):
        raise ValueError("Number of initial protein concentrations provided does not match number of proteins.")
    if (len(Kn_vals) != N):
        raise ValueError("Number of Kd values provided does not match number of proteins.")

    # Set up the symbols for the equation:
    X, x, n = sp.symbols('X x n', real=True)
    P0n = sp.symbols('P0n')
    Kn = sp.symbols('Kn')

    # Define the equation:
    eq = X - x + sp.Sum((P0n * X) / (Kn + X), (n, 1, N))

    # Expand the summation manually
    eq_expand = sp.Add(*[eq.subs({P0n: P0n_val, Kn: Kn_val}) for P0n_val, Kn_val in zip(P0n_vals, Kn_vals)])

    # Substitute the value of 'x'
    eq_subs_x = eq_expand.subs({x: x0})

    # Simplify step by step
    eq_collect = sp.collect(eq_subs_x, X)  # Collect X terms
    eq_simplify = sp.simplify(eq_collect)  # Simplify the expression

    eq_final = sp.simplify(eq_simplify)
    eq_final = sp.denom(eq_final) * eq_final

    X_roots = np.roots(sp.Poly(eq_final, X).all_coeffs())
    X_eq = X_roots[X_roots >= 0][X_roots[X_roots >= 0] <= x0][0]

    return X_eq


## Process RNACompete Data
################################################################################
RNCMPT_org = pd.read_table(baseDir + 'RNACompete/z_scores.txt')
RBP_List = pd.read_table(baseDir + 'RNACompete/rbp_list.txt')

RNCMPT = RNCMPT_org[['7mer'] + [col for col in RNCMPT_org.columns if col.endswith('_setAB')]]
RNCMPT.columns = ['Motif'] + [col.replace('_z_setAB', '') for col in RNCMPT.columns[1:]]
RNCMPT = RNCMPT[['Motif'] + list(RBP_List.RNACompete)]

RNCMPT_norm = RNCMPT.copy()
RNCMPT_norm.iloc[:, 1:] = RNCMPT_norm.iloc[:, 1:].apply(lambda x: featureScale(x, 1, 0))

# RNCMPT_norm.to_csv(baseDir + 'RNACompete/normed_7mer.csv')

RNCMPT_norm['PTBP1'] = RNCMPT_norm[RBP_List[RBP_List.RBP == 'PTBP1'].RNACompete].mean(axis = 1)
RNCMPT_norm = RNCMPT_norm.drop(columns = RBP_List[RBP_List.RBP == 'PTBP1'].RNACompete)

RNCMPT_norm['HuR'] = RNCMPT_norm[RBP_List[RBP_List.RBP == 'HuR'].RNACompete].mean(axis = 1)
RNCMPT_norm = RNCMPT_norm.drop(columns = RBP_List[RBP_List.RBP == 'HuR'].RNACompete)

RNCMPT_norm.columns = ['Motif', 'U2AF2', 'HNRNPC', 'MBNL1', 'RBM24', 'RBM41', 'SNRPA', 'KHDRBS1', 'PCBP1', 'PTBP1', 'HuR']

RBP_IS = pd.DataFrame({
    'RBP': (1/RNCMPT_norm.iloc[:, 1:11].median()).index,
    'IS': (1/RNCMPT_norm.iloc[:, 1:11].median()).values
})


RNCMPT['PTBP1'] = RNCMPT[RBP_List[RBP_List.RBP == 'PTBP1'].RNACompete].mean(axis = 1)
RNCMPT = RNCMPT.drop(columns = RBP_List[RBP_List.RBP == 'PTBP1'].RNACompete)

RNCMPT['HuR'] = RNCMPT[RBP_List[RBP_List.RBP == 'HuR'].RNACompete].mean(axis = 1)
RNCMPT = RNCMPT.drop(columns = RBP_List[RBP_List.RBP == 'HuR'].RNACompete)

RNCMPT.columns = ['Motif', 'U2AF2', 'HNRNPC', 'MBNL1', 'RBM24', 'RBM41', 'SNRPA', 'KHDRBS1', 'PCBP1', 'PTBP1', 'HuR']


## Set up query RNA:
################################################################################
# query:          RNA sequence to test.
# L:              length of the query sequence.
# K:              K-mer
# W:              number of possible K-mers in the query sequence.
# query_W:        all possible K-mers within the query sequence.
# unique_query_W: all unique K-mers within the query sequence.

K = 7

RNA_Target = 'PTBP2'

query = 'GTACTGACCTATATTTTATTTTGTTTTTGTTCCCCAATTCCTTATTTTTT TCTTCTGCATTGCTGTTTCCCTTCCCCATTTCATCCTTTTCCCTGTGTGT TCACCTTCCCTTTCCTTGTCTTTTCCCAAATGCCCATTCCCTTCCCTGTC TTATCCTTTATTTTCCTTGTCCTTGTCTTCATTCCCTGTCTCCATTCCCT A'
query = query.replace('T', 'U').replace(' ', '')
L = len(query)

# L = len(query)
W = L - K + 1

query_W = [query[i:i+K] for i in range(W)]

unique_query_W = pd.DataFrame({'motif': list(set(query_W))})
unique_query_W['cnt'] = unique_query_W['motif'].apply(lambda x: query_W.count(x))

## Set up RBP:
################################################################################
## Used for Figures
U2AF2_KaH = 1/100
U2AF2_KaL = 0.00001
PTBP1_KaH = 1/300
PTBP1_KaL = 0.00001
HNRNPC_KaH = 1/700
HNRNPC_KaL = 0.00001
##

RBP_Ka = pd.DataFrame({
    'Motif': RNCMPT_norm.Motif,
    'U2AF2': RNCMPT_norm.U2AF2,
    'PTBP1': RNCMPT_norm.PTBP1,
    'HNRNPC': RNCMPT_norm.HNRNPC
})
RBP_Ka.U2AF2 = featureScale(RBP_Ka.U2AF2, MAX=U2AF2_KaH, MIN=U2AF2_KaL)
RBP_Ka.PTBP1 = featureScale(RBP_Ka.PTBP1, MAX=PTBP1_KaH, MIN=PTBP1_KaL)
RBP_Ka.HNRNPC = featureScale(RBP_Ka.HNRNPC, MAX=HNRNPC_KaH, MIN=HNRNPC_KaL)

## Visualize affinity distributions for the selected RBPs
################################################################################
plt.figure(figsize = (10, 6))
sns.histplot(RBP_Ka.U2AF2, kde = False, color = 'red', alpha = 0.5, label = 'U2AF2', bins = 250)
sns.histplot(RBP_Ka.PTBP1, kde = False, color = 'blue', alpha = 0.5, label = 'PTBP1', bins = 250)
sns.histplot(RBP_Ka.HNRNPC, kde = False, color = 'black', alpha = 0.5, label = 'HNRNPC', bins = 250)
plt.xlabel('Pseudo Ka')
plt.ylabel('Frequency')
plt.legend()
plt.show()

## Set up simulation for U2AF2:
################################################################################
N = 1
U2AF20 = 500

P0 = [U2AF20]
R0 = 10

## Run simulation for U2AF2:
U2AF2 = pd.DataFrame({
   'motif': RBP_Ka.Motif,
   'Kd': 1/RBP_Ka.U2AF2,
   })

simulation = pd.DataFrame({
   'win_num': range(W),
   'motif': query_W,
   'Req': [None] * W,
   'U2AF2Req': [None] * W
   })

## Loop through each window:
for idx in range(len(unique_query_W)):
  window = unique_query_W.loc[idx, 'motif']
  appearance = unique_query_W.loc[idx, 'cnt']
  eR0 = R0 * appearance
  Kd_U2AF2 = U2AF2.loc[U2AF2.motif == window, 'Kd'].values[0]

  Kd = [Kd_U2AF2]
  Req = allBind(N, eR0, P0, Kd)
  
  U2AF2Req = (U2AF20 * Req)/(Kd_U2AF2 + Req)
 
  simulation.loc[simulation.motif == window, 'Req'] = Req/appearance
  simulation.loc[simulation.motif == window, 'U2AF2Req'] = U2AF2Req/appearance

simulation = simulation.sort_values(by = ["win_num"])

## Calculate per position values:
simulation_per_pos = pd.DataFrame({
    'pos': range(L),
    'nt': list(query),
    'Req': np.zeros(L),
    'U2AF2Req': np.zeros(L),
    'U2AF2/U2AF20': np.zeros(L)
})
update_counts = np.zeros(L)

for idx in range(W):
    update_counts[idx:idx+K] += 1
    simulation_per_pos.loc[idx:idx+K-1, 'Req'] = simulation_per_pos.loc[idx:idx+K-1, 'Req'] + simulation.Req[idx]
    simulation_per_pos.loc[idx:idx+K-1, 'U2AF2Req'] = simulation_per_pos.loc[idx:idx+K-1, 'U2AF2Req'] + simulation.U2AF2Req[idx]

simulation_per_pos['Req'] = simulation_per_pos['Req'] / np.maximum(update_counts, 1)
simulation_per_pos['U2AF2Req'] = simulation_per_pos['U2AF2Req'] / np.maximum(update_counts, 1)

## Heatmap

plot_data_U2AF2 = pd.DataFrame({
    'pos': simulation_per_pos.pos,
    'U2AF2Req': (U2AF20*simulation_per_pos.U2AF2Req/integrate.trapezoid(simulation_per_pos.U2AF2Req,simulation_per_pos.pos)) / (U2AF20/len(query)),
})

plot_data_U2AF2 = plot_data_U2AF2.melt(id_vars='pos', value_vars=['U2AF2Req'], var_name='variable', value_name='value')

heatmap_data = plot_data_U2AF2.pivot(index='pos', columns='variable', values='value')

# Create the heatmap
plt.figure(figsize=(50, 2))
sns.heatmap(heatmap_data.T, cmap='mako', cbar_kws={'label': 'U2AF2 Proportion Fold Over PTBP1'}, vmin=0.5, vmax=2.5)

plt.xlabel('Position')
plt.ylabel('U2AF2 Proportion Fold Over PTBP1')
plt.title('RNA-Protein Complex Proportion by Position with Nucleotide Labels (U2AF20: ' + str(U2AF20) + 'nM)')
plt.xticks(simulation_per_pos['pos'], simulation_per_pos['nt'])

plt.show()

## Set up simulation for U2AF2 and PTBP1:
################################################################################
N = 2
U2AF20 = 500
PTBP10 = 200

P0 = [U2AF20, PTBP10]
R0 = 10

## Run simulation for U2AF2 and PTBP1:
U2AF2 = pd.DataFrame({
   'motif': RBP_Ka.Motif,
   'Kd': 1/RBP_Ka.U2AF2,
   })


PTBP1 = pd.DataFrame({
   'motif': RBP_Ka.Motif,
   'Kd': 1/RBP_Ka.PTBP1,
   })

simulation = pd.DataFrame({
   'win_num': range(W),
   'motif': query_W,
   'Req': [None] * W,
   'U2AF2Req': [None] * W,
   'PTBP1Req': [None] * W
   })


## Loop through each window:
for idx in range(len(unique_query_W)):
  window = unique_query_W.loc[idx, 'motif']
  appearance = unique_query_W.loc[idx, 'cnt']
  eR0 = R0 * appearance
  
  Kd_U2AF2 = U2AF2.loc[U2AF2.motif == window, 'Kd'].values[0]
  Kd_PTBP1 = PTBP1.loc[PTBP1.motif == window, 'Kd'].values[0]

  Kd = [Kd_U2AF2, Kd_PTBP1]
  Req = allBind(N, eR0, P0, Kd)
  
  U2AF2Req = (U2AF20 * Req)/(Kd_U2AF2 + Req)
  PTBP1Req = (PTBP10 * Req)/(Kd_PTBP1 + Req)
 
  simulation.loc[simulation.motif == window, 'Req'] = Req/appearance
  simulation.loc[simulation.motif == window, 'U2AF2Req'] = U2AF2Req/appearance
  simulation.loc[simulation.motif == window, 'PTBP1Req'] = PTBP1Req/appearance

simulation = simulation.sort_values(by = ["win_num"])

## Calculate per position values:
simulation_per_pos = pd.DataFrame({
    'pos': range(L),
    'nt': list(query),
    'Req': np.zeros(L),
    'U2AF2Req': np.zeros(L),
    'U2AF2/U2AF20': np.zeros(L),
    'PTBP1Req': np.zeros(L),
    'PTBP1/PTBP10': np.zeros(L)
})
update_counts = np.zeros(L)

for idx in range(W):
    update_counts[idx:idx+K] += 1
    simulation_per_pos.loc[idx:idx+K-1, 'Req'] = simulation_per_pos.loc[idx:idx+K-1, 'Req'] + simulation.Req[idx]
    simulation_per_pos.loc[idx:idx+K-1, 'U2AF2Req'] = simulation_per_pos.loc[idx:idx+K-1, 'U2AF2Req'] + simulation.U2AF2Req[idx]
    simulation_per_pos.loc[idx:idx+K-1, 'PTBP1Req'] = simulation_per_pos.loc[idx:idx+K-1, 'PTBP1Req'] + simulation.PTBP1Req[idx]

simulation_per_pos['Req'] = simulation_per_pos['Req'] / np.maximum(update_counts, 1)
simulation_per_pos['U2AF2Req'] = simulation_per_pos['U2AF2Req'] / np.maximum(update_counts, 1)
simulation_per_pos['PTBP1Req'] = simulation_per_pos['PTBP1Req'] / np.maximum(update_counts, 1)

## Calculate and normalize probability:
simulation_per_pos['U2AF2/R0'] = (simulation_per_pos['U2AF2Req'])/R0
simulation_per_pos['PTBP1/R0'] = (simulation_per_pos['PTBP1Req'])/R0

## Proportion in heatmap:
plot_data = pd.DataFrame({
    'pos': simulation_per_pos.pos,
    'PTBP1Req': (PTBP10*simulation_per_pos.PTBP1Req/integrate.trapezoid(simulation_per_pos.PTBP1Req,simulation_per_pos.pos))/ (PTBP10/len(query)),
    'U2AF2Req': (U2AF20*simulation_per_pos.U2AF2Req/integrate.trapezoid(simulation_per_pos.U2AF2Req,simulation_per_pos.pos))/ (U2AF20/len(query))
})

# plot_data = plot_data.melt(id_vars='pos', value_vars=['U2AF2Req'], var_name='variable', value_name='value')

plot_data['U2AF20'] = (plot_data.U2AF2Req / plot_data.PTBP1Req)
plot_data = plot_data.melt(id_vars='pos', value_vars=['U2AF20'], var_name='variable', value_name='value')

heatmap_data = plot_data.pivot(index='pos', columns='variable', values='value')

# Create the heatmap
plt.figure(figsize=(50, 2))
sns.heatmap(heatmap_data.T, cmap='mako', cbar_kws={'label': 'U2AF2 Proportion Fold Over PTBP1'}, vmin=0, vmax=2.5)

plt.xlabel('Position')
plt.ylabel('U2AF2 Proportion Fold Over PTBP1')
plt.title('RNA-Protein Complex Proportion by Position with Nucleotide Labels (U2AF20: ' + str(U2AF20) + 'nM; PTBP10: ' + str(PTBP10) + 'nM)')
plt.xticks(simulation_per_pos['pos'], simulation_per_pos['nt'])

plt.show()

## Set up simulation for U2AF2 and HNRNPC:
################################################################################
N = 2
U2AF20 = 500
HNRNPC0 = 200

P0 = [U2AF20, HNRNPC0]
R0 = 10

## Run simulation for U2AF2 and HNRNPC:
U2AF2 = pd.DataFrame({
   'motif': RBP_Ka.Motif,
   'Kd': 1/RBP_Ka.U2AF2,
   })


HNRNPC = pd.DataFrame({
   'motif': RBP_Ka.Motif,
   'Kd': 1/RBP_Ka.HNRNPC,
   })

simulation = pd.DataFrame({
   'win_num': range(W),
   'motif': query_W,
   'Req': [None] * W,
   'U2AF2Req': [None] * W,
   'HNRNPCReq': [None] * W
   })


## Loop through each window:
for idx in range(len(unique_query_W)):
  window = unique_query_W.loc[idx, 'motif']
  appearance = unique_query_W.loc[idx, 'cnt']
  eR0 = R0 * appearance
  
  Kd_U2AF2 = U2AF2.loc[U2AF2.motif == window, 'Kd'].values[0]
  Kd_HNRNPC = HNRNPC.loc[HNRNPC.motif == window, 'Kd'].values[0]

  Kd = [Kd_U2AF2, Kd_HNRNPC]
  Req = allBind(N, eR0, P0, Kd)
  
  U2AF2Req = (U2AF20 * Req)/(Kd_U2AF2 + Req)
  HNRNPCReq = (HNRNPC0 * Req)/(Kd_HNRNPC + Req)
 
  simulation.loc[simulation.motif == window, 'Req'] = Req/appearance
  simulation.loc[simulation.motif == window, 'U2AF2Req'] = U2AF2Req/appearance
  simulation.loc[simulation.motif == window, 'HNRNPCReq'] = HNRNPCReq/appearance

simulation = simulation.sort_values(by = ["win_num"])

## Calculate per position values:
simulation_per_pos = pd.DataFrame({
    'pos': range(L),
    'nt': list(query),
    'Req': np.zeros(L),
    'U2AF2Req': np.zeros(L),
    'U2AF2/U2AF20': np.zeros(L),
    'HNRNPCReq': np.zeros(L),
    'HNRNPC/HNRNPC0': np.zeros(L)
})
update_counts = np.zeros(L)

for idx in range(W):
    update_counts[idx:idx+K] += 1
    simulation_per_pos.loc[idx:idx+K-1, 'Req'] = simulation_per_pos.loc[idx:idx+K-1, 'Req'] + simulation.Req[idx]
    simulation_per_pos.loc[idx:idx+K-1, 'U2AF2Req'] = simulation_per_pos.loc[idx:idx+K-1, 'U2AF2Req'] + simulation.U2AF2Req[idx]
    simulation_per_pos.loc[idx:idx+K-1, 'HNRNPCReq'] = simulation_per_pos.loc[idx:idx+K-1, 'HNRNPCReq'] + simulation.HNRNPCReq[idx]

simulation_per_pos['Req'] = simulation_per_pos['Req'] / np.maximum(update_counts, 1)
simulation_per_pos['U2AF2Req'] = simulation_per_pos['U2AF2Req'] / np.maximum(update_counts, 1)
simulation_per_pos['HNRNPCReq'] = simulation_per_pos['HNRNPCReq'] / np.maximum(update_counts, 1)

## Calculate and normalize probability:
simulation_per_pos['U2AF2/R0'] = (simulation_per_pos['U2AF2Req'])/R0
simulation_per_pos['HNRNPC/R0'] = (simulation_per_pos['HNRNPCReq'])/R0

plot_data = pd.DataFrame({
    'pos': simulation_per_pos.pos,
    'HNRNPCReq': (HNRNPC0*simulation_per_pos.HNRNPCReq/integrate.trapezoid(simulation_per_pos.HNRNPCReq,simulation_per_pos.pos))/ (HNRNPC0/len(query)),
    'U2AF2Req': (U2AF20*simulation_per_pos.U2AF2Req/integrate.trapezoid(simulation_per_pos.U2AF2Req,simulation_per_pos.pos))/ (U2AF20/len(query))
})

## Proportion in heatmap:
plot_data['HNRNPC0'] = plot_data.U2AF2Req / plot_data.HNRNPCReq
plot_data = plot_data.melt(id_vars='pos', value_vars=['HNRNPC0'], var_name='variable', value_name='value')

heatmap_data = plot_data.pivot(index='pos', columns='variable', values='value')

# Create the heatmap
plt.figure(figsize=(50, 2))
sns.heatmap(heatmap_data.T, cmap='mako', cbar_kws={'label': 'U2AF2 Proportion Fold Over HNRNPC'}, vmin=0, vmax=2.5)

plt.xlabel('Position')
plt.ylabel('U2AF2 Proportion Fold Over HNRNPC')
plt.title('RNA-Protein Complex Proportion by Position with Nucleotide Labels (U2AF20: ' + str(U2AF20) + 'nM; HNRNPC: ' + str(PTBP10) + 'nM)')
plt.xticks(simulation_per_pos['pos'], simulation_per_pos['nt'])

plt.show()
