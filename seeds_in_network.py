import pandas as pd

# network assumed to be tsv, gene gene score

# data/networks/entrez_filtered_brain_0.3.tsv
net_path = input("Enter network file path: ")

# 
seed_path = input("Enter seed file path: ")
out_path = input("Enter the output path/file name: ")

df = pd.read_csv(net_path, delimiter='\t')

genes = set()

for gene in df.iloc[:, 0]:
    genes.add(str(gene))

for gene in df.iloc[:, 1]:
    genes.add(str(gene))

fout = open(out_path, "w")
with open(seed_path) as seeds:
    for seed in seeds:
        seed = seed.strip()
        if seed in genes:
            fout.write(seed + '\n')

