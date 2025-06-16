from file_read_backwards import FileReadBackwards

disease = input("rheumatoid_arthritis or schizophrenia: ")
filePath = "output/validation_rankings/" + disease + f"/{"SZ" if disease == "schizophrenia" else "RA"}_cross_validation_"

for i in range(25):
    with FileReadBackwards(f"{filePath}{str(i)}.out", encoding="utf-8") as infile_seeds:
        infile_nonseeds = open(f"cross_validation/{disease}/{disease if disease == "schizophrenia" else "RA"}_non_seeds_{str(i)}.txt")
        nonseeds = {node[:-1] for node in infile_nonseeds}
        outfile = open(f"output/validation_output_labels/{disease}/{"SZ" if disease == "schizophrenia" else "RA"}_validation_labels_{str(i)}.txt", 'w')
        outfile.write("Gene Score Label\n")
        for line in infile_seeds:
            line = line.split()
            node, score = line[0], line[1][:-1]
            outfile.write(f"{node} {score} {'1' if node in nonseeds else '0'}\n")

        outfile.close()