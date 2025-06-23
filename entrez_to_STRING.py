import pandas as pd

# /Users/dkyee/Desktop/Homo_sapiens.gene_info
gene_info_file = input("Enter path to Homo_sapiens.gene_info: ")
gene_info_df = pd.read_csv(gene_info_file, sep='\t', low_memory=False)

# Build mapping from Entrez Gene ID to Symbol
entrez_to_symbol = dict(zip(gene_info_df['GeneID'], gene_info_df['Symbol']))

# output/og/SZ_tissue_0.3_threshold.txt
rankings_file = input("Enter path to rankings file: ")
rankings = open(rankings_file)
mapped_symbols = []

for line in rankings:
    entrezID, score = line.split()
    entrezID = int(entrezID)
    if entrezID in entrez_to_symbol:
        gene_symbol = entrez_to_symbol[entrezID]
        mapped_symbols.append(((entrez_to_symbol[entrezID] if entrezID in entrez_to_symbol else "NULL"), score))

# ADAGIO outputs rankings in reverse order
mapped_symbols.reverse()

# Load alias file (STRING protein aliases)
# ../PPI Networks/Human/Data/9606__protein_aliases.tsv
alias_file_name = input("Enter STRING alias file path/name: ")
alias_df = pd.read_csv(alias_file_name, sep="\t", header=None, names=["protein_id", "alias", "source"])

# symbol_to_string = dict(zip(alias_df["alias"], alias_df["protein_id"]))
symbol_to_string = dict(zip(
    alias_df[alias_df["source"].str.contains("Gene_Symbol|Ensembl_HGNC", na=False)]["alias"],
    alias_df[alias_df["source"].str.contains("Gene_Symbol|Ensembl_HGNC", na=False)]["protein_id"]
))

outfile = open(input("Enter path to outfile: "), 'w')
for symbol, score in mapped_symbols:
    outfile.write(f"{symbol_to_string[symbol] if symbol in symbol_to_string else "NULL"} {score}\n")

print("Done")