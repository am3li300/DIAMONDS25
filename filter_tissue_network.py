# /Users/dkyee/Desktop/brain_humanbase_edges.tsv
input_file = input("Enter path to the brain network edges file: ")

# filtered_brain_1.tsv
output_file = input("Enter path for output filtered file: ")
threshold = float(input("Enter posterior probability threshold (0-1): "))

kept = 0
total = 0

with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
    for line in fin:
        parts = line.strip().split('\t')
        gene1, gene2, prob = parts
        prob = float(prob)
        total += 1
        if prob >= threshold:
            fout.write(f"{gene1}\t{gene2}\t{prob}\n")
            kept += 1

print(f"Kept {kept} out of {total} edges with threshold {threshold}.")