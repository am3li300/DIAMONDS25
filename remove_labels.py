with open(input("Enter input file path: "), "r") as fin:
    fout = open(input("Enter output file path: "), "w")
    for line in fin:
        gene, score, label = line.split(sep=' ')
        fout.write(f"{gene}\t{score}\n")
    fout.close()
