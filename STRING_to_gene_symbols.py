import pandas as pd

def main():
	# Load alias file
	alias_df = pd.read_csv("9606.protein.aliases.v12.0.txt", sep="\t", header=None, names=["protein_id", "alias", "source"])

	# Load nodes_to_expand
	with open("nodes_to_expand.txt") as f:
    		seed_ids = [line.strip() for line in f if line.strip()]

	# Map to HGNC gene symbols
	symbol_map = alias_df[(alias_df["protein_id"].isin(seed_ids)) & (alias_df["source"].str.contains("Ensembl_HGNC|BLAST_UniProt_AC|Ensembl_Gene"))]
	# You can filter for a preferred source if you want only HGNC, for example.

	# Find gene symbols for each protein_id
	gene_dict = symbol_map.groupby("protein_id")["alias"].apply(list).to_dict()

	with open("STRING_to_gene_symbols.txt", 'w') as file:
		for pid in seed_ids:
    			file.write(pid + ': ' + ", ".join(gene_dict.get(pid, [])) + '\n\n')

	print("Saved STRING to gene symbol mapping to STRING_to_gene_symbols.txt")


if __name__ == "__main__":
	main()
