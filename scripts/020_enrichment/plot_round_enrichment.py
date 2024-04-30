import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

filtered_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/010-slxn_qc/'

filtered_seqs = pd.read_csv(filtered_path+'filtered_round'+str(1)+'.csv')
filtered_seqs['round'] = 1
all_filtered_seqs = filtered_seqs

for round_num in range(2,9):
    filtered_seqs = pd.read_csv(filtered_path+'filtered_round'+str(round_num)+'.csv')
    filtered_seqs['round'] = round_num
    all_filtered_seqs = pd.concat([all_filtered_seqs,filtered_seqs])

fig, axs = plt.subplots(2, 4, figsize=(15, 8))
fig.subplots_adjust(hspace=0.5, wspace=0.3)
axs = axs.flatten()

for main_round in range(1, 9):
    print(main_round)
    top_10_seqs = [i[0] for i in Counter(list(all_filtered_seqs[all_filtered_seqs['round'] == main_round].trimmed_seq)).most_common(10)]

    all_counts = []
    for seq in top_10_seqs:
        counts = []
        for round_num in range(1, 9):
            round_seqs = list(all_filtered_seqs[all_filtered_seqs['round'] == round_num].trimmed_seq)
            counts.append(round_seqs.count(seq))
        all_counts.append(counts)

    im = axs[main_round-1].imshow(np.log(1 + np.array(all_counts)), cmap='Purples', interpolation='nearest')
    axs[main_round-1].set_title('Round ' + str(main_round))
    axs[main_round-1].set_xlabel('Round')
    axs[main_round-1].set_ylabel('Top Sequences')
    fig.colorbar(im, ax=axs[main_round-1], label='Log counts')

plt.savefig('/Users/zoeweiss/Desktop/evolution_manuscript/proj-population_dynamics_ligase/data/enrichment_by_round.png')
plt.show()
