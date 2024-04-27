import pandas as pd

raw_path = '/Users/zoeweiss/Desktop/evolution_manuscript/raw_data/2ai_slxn/Sequences/'
filtered_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/010-slxn_qc/'
clusters_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/020-slxn_clusters/'

# Initialize lists to store metrics
total_sequences_list = []
high_quality_sequences_list = []
unique_sequences_list = []
diversity_percentage_list = []
g10_sequence_list = []
number_clusters_list = []

round_num = 8

raw_data = pd.read_csv(raw_path+f'r{round_num}_S{round_num}_L001_R1_001.fastq.gz', sep='\n', header=None).values.reshape(-1, 4)
filtered_data = pd.read_csv(filtered_path+'filtered_round'+str(round_num)+'.csv')
g10_data = pd.read_csv(filtered_path+'filtered_round'+str(round_num)+'_g10_seqs.csv')
clustering_data = pd.read_csv(clusters_path+'clusters_round'+str(round_num)+'.csv')

# Loop through rounds 1 to 8
for round_num in range(1, 9):
    print(round_num)
    
    # Read raw and filtered data
    raw_data = pd.read_csv(raw_path+f'r{round_num}_S{round_num}_L001_R1_001.fastq.gz', sep='\n', header=None).values.reshape(-1, 4)
    filtered_data = pd.read_csv(filtered_path+'filtered_round'+str(round_num)+'.csv')
    g10_data = pd.read_csv(filtered_path+'filtered_round'+str(round_num)+'_g10_seqs.csv')
    

    # Calculate metrics
    total_sequences = len(raw_data)
    high_quality_sequences = len(filtered_data)
    unique_sequences = filtered_data['trimmed_seq'].nunique()
    diversity_percentage = (unique_sequences / high_quality_sequences) * 100
    g10_sequences = len(g10_data)

    # Append metrics to lists
    total_sequences_list.append(total_sequences)
    high_quality_sequences_list.append(high_quality_sequences)
    unique_sequences_list.append(unique_sequences)
    diversity_percentage_list.append(diversity_percentage)
    g10_sequence_list.append(g10_sequences)
    
    if round_num>5:
        clustering_data = pd.read_csv(clusters_path+'clusters_round'+str(round_num)+'.csv')
        number_clusters = max(clustering_data.Cluster)
        number_clusters_list.append(number_clusters)
    else:
        number_clusters_list.append('-')


# Create DataFrame
metrics_df = pd.DataFrame({
    'Round': list(range(1, 9)),
    'Total sequences': total_sequences_list,
    'High quality sequences': high_quality_sequences_list,
    'Unique sequences': unique_sequences_list,
    'Diversity (%)': diversity_percentage_list,
    'Sequences with >10 reads': g10_sequence_list,
    'Number of sequence families': number_clusters_list
})

# Display DataFrame
print(metrics_df)

metrics_df.to_csv(filtered_path+'qc_table.csv')
