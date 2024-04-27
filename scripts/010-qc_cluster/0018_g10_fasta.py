import pandas as pd
from collections import Counter

base_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/010-slxn_qc/'

round_num = 1

for round_num in range(1,9):
    g10 = pd.read_csv(base_path + 'filtered_round'+str(round_num)+'_g10_seqs.csv')
    with open(base_path+'temp_slxn_g10_round'+str(round_num)+'.fasta', 'w') as f:
                j = 0
                for seq in list(g10['0']):
                    f.write(f'>seq{j+1}\n{seq}\n')
                    j+=1