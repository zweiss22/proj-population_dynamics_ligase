import pandas as pd
from collections import Counter
from Bio import pairwise2
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import entropy


# Function to calculate sequence identity
def calculate_similarity(seq1, seq2):
    alignment = pairwise2.align.globalxx(seq1, seq2)
    top_aln = alignment[0]
    align_length = max(len(seq1), len(seq2))
    identity = top_aln[2] / align_length
    return identity

#base_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/010-slxn_qc/'
#seqs = pd.read_csv(base_path+'filtered_round8.csv')

seqs = pd.read_csv('/Users/zoeweiss/Desktop/evolution_manuscript/raw_data/doped_slxn/zip files/filtered_doped_round3.csv')
#seqs = list(seqs[seqs.Ranked_Cluster == 1].Sequence)

seqs = list(seqs.trimmed_seq)

RS1 = 'GAATGCTGCCAACCGTGCGGGCTAATTGGCAGACTGAGCT'

counted_seqs = Counter(seqs)
seqs_g10 = list({key: value for key, value in counted_seqs.items() if value > 10}.keys())


# Initialize counts for each position and each nucleotide
count_matrix = np.zeros((len(seqs_g10[0]), 4))

# Iterate over each sequence to count nucleotide occurrences at each position
for seq in seqs_g10:
    for i, nucleotide in enumerate(seq):
        if nucleotide == 'A':
            count_matrix[i, 0] += 1
        elif nucleotide == 'C':
            count_matrix[i, 1] += 1
        elif nucleotide == 'G':
            count_matrix[i, 2] += 1
        elif nucleotide == 'T':
            count_matrix[i, 3] += 1

# Calculate percentage
total_counts = count_matrix.sum(axis=1)
percentage_matrix = count_matrix / total_counts[:, None] * 100

# Plot heatmap
plt.figure(figsize=(23, 2))
plt.imshow(percentage_matrix.T, cmap='Purples', aspect='auto')
plt.colorbar(label='Percentage')
plt.xlabel('Position', fontsize = 15)
plt.ylabel('Nucleotide', fontsize = 15)
plt.xticks(range(len(seqs_g10[0])), range(1, len(seqs_g10[0]) + 1), fontsize = 15)
plt.yticks(range(4), ['A', 'C', 'G', 'T'], fontsize = 15)
plt.show()


# Calculate entropy for each position
entropy_values = [entropy(count_matrix[i, :]) for i in range(count_matrix.shape[0])]

# Find the indices of the top 5 most conserved spots
most_conserved_indices = sorted(range(len(entropy_values)), key=lambda i: entropy_values[i])[:20]

# Find the indices of the top 5 most variable spots
most_variable_indices = sorted(range(len(entropy_values)), key=lambda i: entropy_values[i], reverse=True)[:20]

print("Top 5 most conserved spots:")
for i in most_conserved_indices:
    print("Position:", i+1, "Entropy:", entropy_values[i])

print("\nTop 5 most variable spots:")
for i in most_variable_indices:
    print("Position:", i+1, "Entropy:", entropy_values[i])