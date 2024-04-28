import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

base_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/010-slxn_qc/'

# Function to calculate sequence identity
def calculate_similarity(seq1, seq2):
    alignment = pairwise2.align.globalxx(seq1, seq2)
    top_aln = alignment[0]
    align_length = max(len(seq1), len(seq2))
    identity = top_aln[2] / align_length
    return identity

overhang = 'ATTCCGCA'
c_overhang = 'TAAGGCGT'
rc_overhang = 'TGCGGAAT'

c_overhang = rc_overhang

# Getting substrings of lengths 4, 5, 6, and 7
substrings_4 = [c_overhang[i:i+4] for i in range(len(c_overhang) - 3)]
substrings_5 = [c_overhang[i:i+5] for i in range(len(c_overhang) - 4)]
substrings_6 = [c_overhang[i:i+6] for i in range(len(c_overhang) - 5)]
substrings_7 = [c_overhang[i:i+7] for i in range(len(c_overhang) - 6)]

substrings = substrings_4+substrings_5+substrings_6+substrings_7

all_fractions = []

for round_num in range(1,9):
    print(round_num)
    
    exact = 0
    close = 0
    
    seqs = pd.read_csv(base_path+'filtered_round'+str(round_num)+'.csv')
    for seq in list(seqs.trimmed_seq):
        if c_overhang in seq:
            exact+=1
        else:
            for sub in substrings:
                if sub in seq:
                    close+=1
                    break
    all_fractions.append([exact/len(seqs), close/len(seqs)])

# Plotting
rounds = range(1, 9)  # 8 rounds
conditions = ['Complete overhang reverse complement', 'Partial overhang reverse complement', ]

plt.figure(figsize=(10, 6))  

for i in range(len(all_fractions[0])):
    plt.plot(rounds, [row[i] for row in all_fractions], label=conditions[i], marker = 'o', linewidth = 3)

plt.xlabel('Round', fontsize = 15)
plt.ylabel('Fraction of reads', fontsize = 15)
plt.legend(fontsize = 15)
plt.grid(True)
plt.show()

pd.DataFrame(all_fractions).to_csv('/Users/zoeweiss/Desktop/evolution_manuscript/proj-population_dynamics_ligase/scripts/040_overhang/overhang_fractions.csv')
