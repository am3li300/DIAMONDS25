import pandas as pd

def main():
	# Load alias file (STRING protein aliases)
	# ../PPI Networks/Human/Data/9606__protein_aliases.tsv
	alias_file_name = input("Enter STRING alias file path/name: ")
	alias_df = pd.read_csv(alias_file_name, sep="\t", header=None, names=["protein_id", "alias", "source"])

	# Load nodes file (your list of STRING protein IDs)
	# data/20_data_drug_schizophrenia.txt
	file_name = input("Enter tsv/txt file path/name: ")
	seed_ids = []
	with open(file_name, 'r') as f:
		for line in f:
			protein = line.strip().split()[0]
			seed_ids.append(protein)

	# Load official gene info (from NCBI)
	# /Users/dkyee/Desktop/Homo_sapiens.gene_info
	gene_info_file = input("Enter path to Homo_sapiens.gene_info: ")
	gene_info_df = pd.read_csv(gene_info_file, sep='\t', comment='#', header=0, low_memory=False)

	# Build a set for fast lookup
	canonical_symbols = set(gene_info_df['Symbol'])

	# Filter alias file for symbols labeled as "Entrez"
	symbol_map = alias_df[
		(alias_df["protein_id"].isin(seed_ids)) &
		(alias_df["source"].str.contains("Entrez"))
	]

	# Group symbols by STRING protein ID
	gene_dict = symbol_map.groupby("protein_id")["alias"].apply(list).to_dict()

	# For each STRING protein, keep the first alias that matches the NCBI canonical symbol list
	with open("STRING_to_gene_symbols.txt", 'w') as file:
		for pid in seed_ids:
			aliases = gene_dict.get(pid, [])
			# Find the first alias that's in the canonical set
			mapped_symbol = next((symbol for symbol in aliases if symbol in canonical_symbols), None)
			if mapped_symbol:
				# file.write(f"{pid} {mapped_symbol}\n")
				file.write(f"{mapped_symbol}\n")
			else:
				print(f"No valid symbol found for {pid}")

	print("Saved STRING to gene symbol mapping to STRING_to_gene_symbols.txt")

if __name__ == "__main__":
	main()
