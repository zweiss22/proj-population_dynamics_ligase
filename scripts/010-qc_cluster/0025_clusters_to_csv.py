import pandas as pd
base_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/010-slxn_qc/'
clusters_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/020-slxn_clusters/'
filtered_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/slxn_qc/'

round_num = 8

for round_num in range(6,9):

    seqs = open(base_path+'temp_slxn_g10_round'+str(round_num)+'.fasta').readlines()
    # Read the contents of the aux file
    with open(clusters_path+'cluster_'+str(round_num)+'.aux', 'r') as file:
        aux_data = file.readlines()

    # Initialize lists to store sequence indices and their corresponding clusters
    indices = []
    clusters = []

    # Parse the lines in the aux file to extract indices and clusters
    for line in aux_data:
        parts = line.split()
        index = int(parts[6])  # Extract sequence index
        cluster = int(parts[1].split(':')[0])  # Extract cluster
        indices.append(index)
        clusters.append(cluster)

    # Create a DataFrame from the lists
    df = pd.DataFrame({'Sequence Index': indices, 'Cluster': clusters})

    seqs_df = pd.DataFrame({'Sequence Index': [int(elem.replace(">seq", "").rstrip("\n")) for elem in seqs[::2]], 'Sequence': seqs[1::2]})
    assigned_clusters = seqs_df.merge(df, on = 'Sequence Index')

    filtered_seqs = pd.read_csv(base_path+'filtered_round'+str(round_num)+'.csv')

    # Count the occurrences of each sequence and create a new DataFrame with ordered cluster numbers
    count_df = pd.DataFrame({"Sequence": filtered_seqs.trimmed_seq})['Sequence'].value_counts().reset_index()
    count_df.columns = ['Sequence', 'Counts']

    assigned_clusters['Sequence'] = assigned_clusters['Sequence'].str.rstrip('\n')
    assigned_clusters = assigned_clusters.merge(count_df, on = 'Sequence')
    assigned_clusters = assigned_clusters.sort_values(by='Counts', ascending=False)

    assigned_clusters['new_cluster'] = assigned_clusters.groupby('Cluster')['Counts'].rank(ascending=False, method='first').astype(int)


    assigned_clusters.to_csv(clusters_path+'clusters_round'+str(round_num)+'.csv')