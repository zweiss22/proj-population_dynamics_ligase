import pandas as pd
import os
import subprocess

base_path = '/Users/zoeweiss/Desktop/evolution_manuscript/2ai_slxn/Sequences/'

filtered_combined = []
for round_num in range(1,9):
    csv_data = pd.read_csv(base_path+'filtered_round'+str(round_num)+'.csv')
    filtered_combined.append(list(csv_data.trimmed_seq))
    
rounds = filtered_combined

#Cluster sequences in each group using Clustal Omega
for i, round_num in enumerate(rounds):
    # Prepare a temporary FASTA file for the cluster
    with open(f'temp_slxn_cluster_{i+1}.fasta', 'w') as f:
        for j, seq in enumerate(list(set(round_num))):
            f.write(f'>seq{j+1}\n{seq}\n')

# Run Clustal Omega on the round
round_num = 1
for round_num in range(1,9):
    cmd = 'clustalo -i temp_slxn_cluster_'+str(round_num)+'.fasta --clustering-out=cluster_'+str(round_num)+'.aux  --cluster-size=500'
    subprocess.run(cmd, shell=True)


