from file_read_backwards import FileReadBackwards
import os

# cross_validation/validation_rankings/
directory_path = input("Enter directory path for input ranking: ")

# cross_validation/validation_output_labels/
outfolder_path = input("Enter output folder path: ")
reverse_flag = int(input("Enter 1 to read file in reverse, 0 otherwise: "))

# Does not iterate through directory in order
for fileName in os.listdir(directory_path):
    i = int(fileName[-5])

    # Build file paths
    if reverse_flag:
        with FileReadBackwards(os.path.join(directory_path, fileName)) as infile_ranking:
            infile_nonseeds = open("cross_validation/partitioning/schizophrenia_STRING/3_schizophrenia_non_seeds_{0}.txt".format(i))
            infile_seeds = open("cross_validation/partitioning/schizophrenia_STRING/3_schizophrenia_new_seeds_{0}.txt".format(i))
            outfile = open("{0}/SZ_validation_labels_{1}".format(outfolder_path, i), 'w')

            # Prepare seed/nonseed sets
            nonseeds = {line.strip() for line in infile_nonseeds}
            seeds = {line.strip() for line in infile_seeds}
            outfile.write("Gene Score Label\n")

            for line in infile_ranking:
                node, score = line.split()[:2]
                if node in seeds:
                    continue
                outfile.write(f"{node} {score} {'1' if node in nonseeds else '0'}\n")

            # Close files
            infile_ranking.close()
            infile_nonseeds.close()
            infile_seeds.close()
            outfile.close()

    else:
        infile_ranking = open(os.path.join(directory_path, fileName))
        infile_nonseeds = open("cross_validation/partitioning/schizophrenia_STRING/3_schizophrenia_non_seeds_{0}.txt".format(i))
        infile_seeds = open("cross_validation/partitioning/schizophrenia_STRING/3_schizophrenia_new_seeds_{0}.txt".format(i))
        outfile = open("{0}/SZ_validation_labels_{1}".format(outfolder_path, i), 'w')

        # Prepare seed/nonseed sets
        nonseeds = {line.strip() for line in infile_nonseeds}
        seeds = {line.strip() for line in infile_seeds}
        outfile.write("Gene Score Label\n")

        for line in infile_ranking:
            node, score = line.split()[:2]
            if node in seeds:
                continue
            outfile.write(f"{node} {score} {'1' if node in nonseeds else '0'}\n")

        # Close files
        infile_ranking.close()
        infile_nonseeds.close()
        infile_seeds.close()
        outfile.close()