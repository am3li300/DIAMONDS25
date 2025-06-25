import pandas as pd

# network assumed to be tsv, gene gene score
path = input("Enter network file path: ")
df = pd.read_csv(path, delimiter='\t')

unique_genes = set()

for gene in df.iloc[:, 0]:
    if gene not in unique_genes:
        unique_genes.add(gene)

for gene in df.iloc[:, 1]:
    if gene not in unique_genes:
        unique_genes.add(gene)

print(f"There are {len(unique_genes)} genes in the network")