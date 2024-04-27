import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

filtered_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/010-slxn_qc/'
clustered_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/020-slxn_clusters/'

final_clusters = pd.read_csv(clustered_path + 'clusters_round'+str(8)+'.csv')

all_peak_counts = []

for cluster in range(1,9):
    print(cluster)
    peak = list(final_clusters[final_clusters.Ranked_Cluster == cluster].Sequence)[0]
    peak_counts_list = []
    for round_num in range(1,9):
        print(round_num)
        filtered_seqs = pd.read_csv(filtered_path+'filtered_round'+str(round_num)+'.csv')
        peak_counts = list(filtered_seqs.trimmed_seq).count(peak)
        peak_counts_list.append(peak_counts/len(filtered_seqs))
    all_peak_counts.append(peak_counts_list)

pd.DataFrame(all_peak_counts).to_csv('/Users/zoeweiss/Desktop/evolution_manuscript/proj-population_dynamics_ligase/data/peak_enrichment_fraction.csv')

# Plotting each row as a line
plt.figure(figsize=(10, 6))

for index, row in pd.DataFrame(all_peak_counts).iterrows():
    plt.plot(row, label=f'Family {1+index}')

# Adding labels and title
plt.xlabel('Round', fontsize=18)
plt.ylabel('Fraction abundance', fontsize=18)
plt.xticks(np.arange(8), np.arange(1, 9), fontsize=18)

plt.yticks(fontsize=18)
plt.legend(fontsize=18)

plt.savefig('/Users/zoeweiss/Desktop/evolution_manuscript/proj-population_dynamics_ligase/data/family_enrichment.png')
