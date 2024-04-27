import pandas as pd
from collections import Counter

base_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/010-slxn_qc/'

def keep_elements_at_least_10_times(lst):
    # Count occurrences of each element
    counts = Counter(lst)
    
    # Filter elements that appear at least 10 times
    result = [elem for elem, count in counts.items() if count >= 10]
    
    return result

for round_num in range(1,9):
    print(round_num)
    filtered_seqs = pd.read_csv(base_path + 'filtered_round'+str(round_num)+'.csv')
    result = keep_elements_at_least_10_times(filtered_seqs.trimmed_seq)
    pd.DataFrame(result).to_csv(base_path + 'filtered_round'+str(round_num)+'_g10_seqs.csv')
