import pandas as pd

# /Users/dkyee/Desktop/Homo_sapiens.gene_info
gene_info_file = input("Enter path to Homo_sapiens.gene_info: ")
gene_info_df = pd.read_csv(gene_info_file, sep='\t', low_memory=False)

# Build mapping from gene symbol to Entrez Gene ID
symbol_to_entrez = dict(zip(gene_info_df['Symbol'], gene_info_df['GeneID']))

seed_file = input("Enter path to STRING_to_gene_symbols text file: ")
seeds = open(seed_file)
outfile = open("entrez_seeds.tsv", 'w')

for line in seeds:
    STRING, gene_symbol = line.split()
    entrez_id = symbol_to_entrez[gene_symbol]
    outfile.write(str(entrez_id) + '\n')

print("Saved as entrez_seeds.tsv")







# [entrez gene id 1][entrez gene id 2][posterior prob.]