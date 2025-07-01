from collections import defaultdict

file = open(input("Enter fasta sequence file: "))
mapping = defaultdict(string)
gene = ""
sequence = []
for line in file:
    if line[0] == '>':
        if gene:
            mapping[gene] = "".join(sequence)
            sequence = []

        x = "".join(line.split('|')[2])
        humanIndx = x.index("_HUMAN")
        gene = x[:humanIndx]


    else:
        sequence.append(line.strip())

mapping[gene] = "".join(sequence)