################################################################################
## RNA-RBP Modeling 
## Written by Soon Yi
## Created: 2023-12-11
## Last Edited: 2023-12-12
################################################################################


## Load libraries
################################################################################
import pandas as pd
import numpy as np
import random
from scipy.integrate import trapz
from sympy import symbols, expand, collect, Poly
import matplotlib.pyplot as plt
import seaborn as sns
################################################################################

## Set basic parameters and functions:
################################################################################
baseDir = '~/Desktop/Genomics/Specificity/Data_Final/Model_RBP/'

pd.set_option('display.float_format', lambda x: '%.5f' % x)

def featureScale(X, MAX, MIN):
    X = np.array(X)  
    X_min = np.nanmin(X)
    X_max = np.nanmax(X)
    return MIN + ((X - X_min) * (MAX - MIN)) / (X_max - X_min)

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
L = 100

query = generateRNA(L=L, 
                    A=0.2,
                    G=0.2,
                    C=0.2,
                    U=0.4)

query = 'UAGCGGUGCGAUUGGCCCGUGGACCGCGUUUUUGCACUCAUCGUUUCGCACUAAGUACAUAUAGUUGCGACAAAGCCGCUUAUGAGUUGGGGGUAUAUUC'

# L = len(query)
W = L - K + 1

query_W = [query[i:i+K] for i in range(W)]

unique_query_W = pd.DataFrame({'motif': list(set(query_W))})
unique_query_W['cnt'] = unique_query_W['motif'].apply(lambda x: query_W.count(x))
################################################################################

## Set up model RBPs:
################################################################################
Ka_H = 100
Ka_L = 10

# Read in the model RBPs
ModelHH = pd.read_csv(baseDir + 'Model_HH.csv').sort_values(by='Motif')
ModelHL = pd.read_csv(baseDir + 'Model_HL.csv').sort_values(by='Motif')
ModelLH = pd.read_csv(baseDir + 'Model_LH.csv').sort_values(by='Motif')
ModelLL = pd.read_csv(baseDir + 'Model_LL.csv').sort_values(by='Motif')

# Combine the models into a single DataFrame
Model_Ka = pd.DataFrame({
    'Motif': ModelHH['Motif'],
    'HH': ModelHH['HH'],
    'HL': ModelHL['HL'],
    'LH': ModelLH['LH'],
    'LL': ModelLL['LL']
})

# Apply feature scaling
Model_Ka['HH'] = featureScale(Model_Ka['HH'], MAX=Ka_H, MIN=1)
Model_Ka['HL'] = featureScale(Model_Ka['HL'], MAX=Ka_H, MIN=1)
Model_Ka['LH'] = featureScale(Model_Ka['LH'], MAX=Ka_L, MIN=1)
Model_Ka['LL'] = featureScale(Model_Ka['LL'], MAX=Ka_L, MIN=1)
################################################################################

## Set up simulation for four RBPs:
################################################################################
N = 4
HH0 = 100
HL0 = 100
LH0 = 100
LL0 = 100

P0 = [HH0, HL0, LH0, LL0]
R0 = 10

HH = pd.DataFrame({
   'motif': Model_Ka.Motif,
   'Kd': 1/Model_Ka['HH'],
   })

HL = pd.DataFrame({
   'motif': Model_Ka.Motif,
   'Kd': 1/Model_Ka['HL'],
   })

LH = pd.DataFrame({
   'motif': Model_Ka.Motif,
   'Kd': 1/Model_Ka['LH'],
   })

LL = pd.DataFrame({
   'motif': Model_Ka.Motif,
   'Kd': 1/Model_Ka['LL'],
   })

simulation = pd.DataFrame({
   'win_num': range(W),
   'motif': query_W,
   'Req': [None] * W,
   'HHReq': [None] * W,
   'HLReq': [None] * W,
   'LHReq': [None] * W,
   'LLReq': [None] * W
   })

## Loop through each window:
for idx in range(len(unique_query_W)):
  window = unique_query_W.loc[idx, 'motif']
  appearance = unique_query_W.loc[idx, 'cnt']
  eR0 = R0 * appearance
  Kd_HH = HH.loc[HH['motif'] == window, 'Kd'].values[0]
  Kd_HL = HL.loc[HH['motif'] == window, 'Kd'].values[0]
  Kd_LH = LH.loc[HH['motif'] == window, 'Kd'].values[0]
  Kd_LL = LL.loc[HH['motif'] == window, 'Kd'].values[0]
  
  Kd = [Kd_HH, Kd_HL, Kd_LH, Kd_LL]
  Req = allBind(N, eR0, P0, Kd)
  
  HHReq = (HH0 * Req)/(Kd_HH + Req)
  HLReq = (HL0 * Req)/(Kd_HL + Req)
  LHReq = (LH0 * Req)/(Kd_LH + Req)
  LLReq = (LL0 * Req)/(Kd_LL + Req)

  simulation.loc[simulation.motif == window, 'Req'] = Req/appearance
  simulation.loc[simulation.motif == window, 'HHReq'] = HHReq/appearance
  simulation.loc[simulation.motif == window, 'HLReq'] = HLReq/appearance
  simulation.loc[simulation.motif == window, 'LHReq'] = LHReq/appearance
  simulation.loc[simulation.motif == window, 'LLReq'] = LLReq/appearance

simulation = simulation.sort_values(by = ["win_num"])


simulation = simulation[simulation.columns[2:]].applymap(lambda x: x.real)

## Calculate per position values:
simulation_per_pos = pd.DataFrame({
    'pos': range(L),
    'nt': list(query),
    'Req': np.zeros(L),
    'HHReq': np.zeros(L),
    'HLReq': np.zeros(L),
    'LHReq': np.zeros(L),
    'LLReq': np.zeros(L),
    'HH/HH0': np.zeros(L),
    'HL/HL0': np.zeros(L),
    'LH/LH0': np.zeros(L),
    'LL/LL0': np.zeros(L),
})
update_counts = np.zeros(L)

for idx in range(W):
    update_counts[idx:idx+K] += 1
    simulation_per_pos.loc[idx:idx+K-1, 'Req'] = simulation_per_pos.loc[idx:idx+K-1, 'Req'] + simulation.Req[idx]
    simulation_per_pos.loc[idx:idx+K-1, 'HHReq'] = simulation_per_pos.loc[idx:idx+K-1, 'HHReq'] + simulation.HHReq[idx]
    simulation_per_pos.loc[idx:idx+K-1, 'HLReq'] = simulation_per_pos.loc[idx:idx+K-1, 'HLReq'] + simulation.HLReq[idx]
    simulation_per_pos.loc[idx:idx+K-1, 'LHReq'] = simulation_per_pos.loc[idx:idx+K-1, 'LHReq'] + simulation.LHReq[idx]
    simulation_per_pos.loc[idx:idx+K-1, 'LLReq'] = simulation_per_pos.loc[idx:idx+K-1, 'LLReq'] + simulation.LLReq[idx]
    

simulation_per_pos['Req'] = simulation_per_pos['Req'] / np.maximum(update_counts, 1)
simulation_per_pos['HHReq'] = simulation_per_pos['HHReq'] / np.maximum(update_counts, 1)
simulation_per_pos['HLReq'] = simulation_per_pos['HLReq'] / np.maximum(update_counts, 1)
simulation_per_pos['LHReq'] = simulation_per_pos['LHReq'] / np.maximum(update_counts, 1)
simulation_per_pos['LLReq'] = simulation_per_pos['LLReq'] / np.maximum(update_counts, 1)

## Calculate and normalize probability:
simulation_per_pos['HH/R0'] = (simulation_per_pos['HHReq'])/R0
simulation_per_pos['HL/R0'] = (simulation_per_pos['HLReq'])/R0
simulation_per_pos['LH/R0'] = (simulation_per_pos['LHReq'])/R0
simulation_per_pos['LL/R0'] = (simulation_per_pos['LLReq'])/R0

## Make dataframe for plotting:
################################################################################
plot_data = pd.DataFrame({
    'pos': simulation_per_pos.pos,
    'HHReq': (simulation_per_pos.HHReq/R0)/sum(simulation_per_pos.HHReq/R0),
    'HLReq': (simulation_per_pos.HLReq/R0)/sum(simulation_per_pos.HLReq/R0),
    'LHReq': (simulation_per_pos.LHReq/R0)/sum(simulation_per_pos.LHReq/R0),
    'LLReq': (simulation_per_pos.LLReq/R0)/sum(simulation_per_pos.LLReq/R0)
})

plot_data = plot_data.melt(id_vars='pos', value_vars=['HHReq', 'HLReq', 'LHReq', 'LLReq'], var_name='variable', value_name='value')

plt.figure(figsize=(20, 2.5))
sns.lineplot(data=plot_data, x='pos', y='value', hue='variable')

plt.ylim(0, 1)
plt.xticks(simulation_per_pos['pos'], simulation_per_pos['nt'])

plt.xlabel('Position')
plt.ylabel('[PReq]/[R0]')
plt.title('Four RBPs Probability by Position with Nucleotide Labels')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
plt.show()
