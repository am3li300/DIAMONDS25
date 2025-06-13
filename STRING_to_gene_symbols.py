import pandas as pd

def main():
	# Load alias file
	alias_file_name = input("Enter STRING alias file path/name: ")
	alias_df = pd.read_csv(alias_file_name, sep="\t", header=None, names=["protein_id", "alias", "source"])

	# Load nodes file
	# ../PPI Networks/Human/
	file_name = input("Enter tsv/txt file path/name: ")
	seed_ids = []
	with open(file_name, 'r') as f:
		for line in f:
			protein = line.strip().split()[0]
			seed_ids.append(protein)


	# Map to HGNC gene symbols
	symbol_map = alias_df[
		(alias_df["protein_id"].isin(seed_ids)) &
		(alias_df["source"].str.contains("Ensembl_HGNC|BLAST_UniProt_AC|Ensembl_Gene"))
	]

	# Find gene symbols for each protein_id
	gene_dict = symbol_map.groupby("protein_id")["alias"].apply(list).to_dict()

	with open("STRING_to_gene_symbols.txt", 'w') as file:
		for pid in seed_ids:
			file.write(pid + ': ' + ", ".join(gene_dict.get(pid, [])) + '\n\n')

	print("Saved STRING to gene symbol mapping to STRING_to_gene_symbols.txt")

if __name__ == "__main__":
	main()
