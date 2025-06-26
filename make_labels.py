from file_read_backwards import FileReadBackwards

"""

for i in range(int(input("Enter number of files: "))):
    # with FileReadBackwards(f"{filePath}{str(i)}.out", encoding="utf-8") as infile_ranking:

    # cross_validation/validation_rankings/SZ_humanbase/3_SZ_cross_validation_0.out
    infile_ranking = open(input(f"Enter file path to gene rankings file {i}: "))

    # cross_validation/partitioning/SZ_humanbase/3_schizophrenia_non_seeds_0.txt
    infile_nonseeds = open(input(f"Enter file path to nonseeds file {i}: "))
    nonseeds = {node[:-1] for node in infile_nonseeds}

    # cross_validation/partitioning/SZ_humanbase/3_schizophrenia_new_seeds_0.txt
    infile_seeds = open(input(f"Enter file path to seeds file {i}: "))
    seeds = {node[:-1] for node in infile_seeds}

    # cross_validation/validation_output_labels/SZ_humanbase/SZ_validation_labels_0.txt
    outfile = open(input(f"Enter file path to output file: "), 'w')
    outfile.write("Gene Score Label\n")

    for line in infile_ranking:
        line = line.split()
        node, score = line[0], line[1]
        if node in seeds:
            continue
            
        outfile.write(f"{node} {score} {'1' if node in nonseeds else '0'}\n")

    outfile.close()

"""

# User enters base path, prefix, and number of files
base_path = input("Enter base path (e.g., cross_validation): ").rstrip('/')
dataset = input("Enter dataset name (e.g., SZ_humanbase): ")
prefix = input("Enter folds number (e.g., 3): ")
n_files = int(input("Enter number of files: "))
reverse_flag = int(input("Enter 1 to read file in reverse, 0 otherwise: "))

for i in range(n_files):
    # Build file paths
    if reverse_flag:
        with FileReadBackwards(f"{base_path}/validation_rankings/{dataset}/{prefix}_SZ_cross_validation_{i}.out") as infile_ranking:
            infile_nonseeds = open(f"{base_path}/partitioning/{dataset}/{prefix}_schizophrenia_non_seeds_{i}.txt")
            infile_seeds = open(f"{base_path}/partitioning/{dataset}/{prefix}_schizophrenia_new_seeds_{i}.txt")
            outfile = open(f"{base_path}/validation_output_labels/{dataset}/SZ_validation_labels_{i}.txt", 'w')

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
        infile_ranking = open(f"{base_path}/validation_rankings/{dataset}/{prefix}_SZ_cross_validation_{i}.out")
        infile_nonseeds = open(f"{base_path}/partitioning/{dataset}/{prefix}_schizophrenia_non_seeds_{i}.txt")
        infile_seeds = open(f"{base_path}/partitioning/{dataset}/{prefix}_schizophrenia_new_seeds_{i}.txt")
        outfile = open(f"{base_path}/validation_output_labels/{dataset}/SZ_validation_labels_{i}.txt", 'w')

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