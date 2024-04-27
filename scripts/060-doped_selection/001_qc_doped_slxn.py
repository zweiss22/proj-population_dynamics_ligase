import pandas as pd
from collections import Counter

base_path = '/Users/zoeweiss/Desktop/evolution_manuscript/doped_slxn/zip files/'
out_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/'

# Define start and end sequences
start_seq = 'ACGGACAGCG'
end_seq = 'CGCTGTCCTTTTTTGGCTAAGGGACCTACCG'

def calculate_error_rate(qual_string):
    """
    Function to calculate error rate from quality string.
    """
    # ASCII to Phred score mapping
    phred_score = lambda q: ord(q) - 33
    
    # Convert each quality score to a Phred score
    scores = [phred_score(q) for q in qual_string]
    
    # Calculate the error rate
    error_rate = sum(10 ** (-q / 10) for q in scores) / len(scores)
    
    return error_rate

def reverse_complement(seq):
    """
    Function to compute the reverse complement of a DNA sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(seq))

# Function to trim sequences
def trim_sequence(seq):
    start_index = seq.find(start_seq) + len(start_seq)
    end_index = seq.find(end_seq)
    if start_index >= 0 and end_index >= 0:
        return seq[start_index:end_index]
    else:
        return None
    
for round_num in range(1,4):
    print(round_num)
    path = base_path+'/doped21-r'+str(round_num)+'_S'+str(9+round_num)+'_L001_R1_001.fastq.gz'

    # Read the data
    data = pd.read_csv(path, sep='\n', header=None).values.reshape(-1, 4)
    data = pd.DataFrame(data, columns=['read_id', 'seq', '+', 'qual'])

    # Calculate error rates
    data['error_rate'] = data['qual'].apply(calculate_error_rate)

    # Filter data based on error rate
    filtered_data = data[data['error_rate'] < 0.01]
    filtered_data = filtered_data.drop(columns=['error_rate'])

    # Take the RC
    filtered_data['reverse_complement'] = filtered_data['seq'].apply(reverse_complement)

    # Filter sequences based on start and end patterns
    filtered_data = filtered_data[filtered_data['reverse_complement'].str.startswith(start_seq) & 
                                  filtered_data['reverse_complement'].str.endswith(end_seq)]

    # Apply trim_sequence function to 'seq' column
    filtered_data['trimmed_seq'] = filtered_data['reverse_complement'].apply(trim_sequence)

    # Drop rows where the trimming failed (i.e., start or end sequence not found)
    filtered_data = filtered_data.dropna(subset=['trimmed_seq'])

    filtered_data['round']=round_num
    
    filtered_data.to_csv(out_path+'filtered_doped_round'+str(round_num)+'.csv')
