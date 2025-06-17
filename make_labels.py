from file_read_backwards import FileReadBackwards

partition_flag = int(input("Enter 1 for 2-fold partitioning, 0 for 3-fold partitioning: "))
prefix = '' if partition_flag else "3_"

disease = input("rheumatoid_arthritis or schizophrenia: ")
filePath = "output/validation_rankings/" + disease + f"/{prefix}{"SZ" if disease == "schizophrenia" else "RA"}_cross_validation_"

for i in range(int(input("Enter number of files: "))):
    with FileReadBackwards(f"{filePath}{str(i)}.out", encoding="utf-8") as infile_ranking:
        infile_nonseeds = open(f"cross_validation/{disease}/{prefix}{disease if disease == "schizophrenia" else "RA"}_non_seeds_{str(i)}.txt")
        nonseeds = {node[:-1] for node in infile_nonseeds}

        outfile = open(f"output/validation_output_labels/{disease}/{prefix}{"SZ" if disease == "schizophrenia" else "RA"}_validation_labels_{str(i)}.txt", 'w')
        outfile.write("Gene Score Label\n")

        for line in infile_ranking:
            line = line.split()
            node, score = line[0], line[1][:-1]
            outfile.write(f"{node} {score} {'1' if node in nonseeds else '0'}\n")

        outfile.close()