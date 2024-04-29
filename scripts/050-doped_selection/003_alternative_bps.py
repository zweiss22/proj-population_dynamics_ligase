import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
base_path = '/Users/zoeweiss/Desktop/evolution_manuscript/raw_data/doped_slxn/zip files/'
doped_r3 = pd.read_csv(base_path+'filtered_doped_round3.csv')
bps = [[6,32],[7,31],[8,30],[9,29],[10,28],[11,27],[12,26]]
seqs = list(doped_r3.trimmed_seq)
import matplotlib.pyplot as plt
import numpy as np

# Initialize the subplot
fig, axs = plt.subplots(2, 4, figsize=(12, 6))

# Iterate through each base pair
for bp, ax in zip(range(len(bps)), axs.flat):
    bp_combos = []
    for seq in range(len(seqs)):
        bp_combos.append([seqs[seq][bps[bp][0] - 1], seqs[seq][bps[bp][1] - 1]])

    # Count the occurrences of each unique combination
    combo_counts = {}
    total_combos = len(bp_combos)
    for combo in bp_combos:
        combo_key = tuple(combo)
        combo_counts[combo_key] = combo_counts.get(combo_key, 0) + 1

    # Calculate the fraction of each combination
    combo_fractions = {combo: count / total_combos for combo, count in combo_counts.items()}

    # Create a matrix to represent the heatmap
    bases = ['A', 'C', 'G', 'T']
    heatmap_matrix = np.zeros((len(bases), len(bases)))
    for i, base1 in enumerate(bases):
        for j, base2 in enumerate(bases):
            combo_count = combo_counts.get((base1, base2), 0)
            #heatmap_matrix[i, j] = combo_count / total_combos
            heatmap_matrix[i, j] = np.log(combo_count)

    # Plot the heatmap in the current subplot
    im = ax.imshow(heatmap_matrix, cmap='Purples', interpolation='nearest')
    ax.set_xticks(np.arange(len(bases)))
    ax.set_yticks(np.arange(len(bases)))
    ax.set_xticklabels(bases)
    ax.set_yticklabels(bases)
    ax.set_xlabel('Second Base')
    ax.set_ylabel('First Base')
    ax.set_title('Base Pair ' + str(bp + 1))
    
    # Add colorbar for each subplot
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    #cbar.set_label('Fraction of Base Pair Combinations')

# Adjust layout
plt.tight_layout()

# Show the plot
plt.show()
