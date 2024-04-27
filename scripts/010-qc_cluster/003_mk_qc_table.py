import pandas as pd

raw_path = '/Users/zoeweiss/Desktop/evolution_manuscript/2ai_slxn/Sequences/'
filtered_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/'

# Initialize lists to store metrics
total_sequences_list = []
high_quality_sequences_list = []
unique_sequences_list = []
diversity_percentage_list = []

# Loop through rounds 1 to 8
for round_num in range(1, 9):
    print(round_num)
    
    # Read raw and filtered data
    raw_data = pd.read_csv(raw_path+f'r{round_num}_S{round_num}_L001_R1_001.fastq.gz', sep='\n', header=None).values.reshape(-1, 4)
    filtered_data = pd.read_csv(in_path+f'filtered_round{round_num}.csv')

    # Calculate metrics
    total_sequences = len(raw_data)
    high_quality_sequences = len(filtered_data)
    unique_sequences = filtered_data['trimmed_seq'].nunique()
    diversity_percentage = (unique_sequences / total_sequences) * 100

    # Append metrics to lists
    total_sequences_list.append(total_sequences)
    high_quality_sequences_list.append(high_quality_sequences)
    unique_sequences_list.append(unique_sequences)
    diversity_percentage_list.append(diversity_percentage)

# Create DataFrame
metrics_df = pd.DataFrame({
    'Round': list(range(1, 9)),
    'Total sequences': total_sequences_list,
    'High quality sequences': high_quality_sequences_list,
    'Unique sequences': unique_sequences_list,
    'Diversity (%)': diversity_percentage_list
})

# Display DataFrame
print(metrics_df)
metrics_df.to_csv(filtered_path+'qc_table.csv')

