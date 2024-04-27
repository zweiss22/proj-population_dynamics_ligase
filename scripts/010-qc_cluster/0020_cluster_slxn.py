import subprocess

base_path = '/Users/zoeweiss/Desktop/evolution_manuscript/results/'

for round_num in range(1,9):
    print(round_num)
    cmd = 'clustalo -i '+base_path+'010-slxn_qc/temp_slxn_g10_round'+str(round_num)+'.fasta --clustering-out='+base_path+'020-slxn_clusters/cluster_'+str(round_num)+'.aux  --cluster-size=500'
    subprocess.run(cmd, shell=True)
