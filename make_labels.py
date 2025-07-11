from file_read_backwards import FileReadBackwards
import os

def make_labels_for_file():
    ranking_file = open(input("Enter file path for input ranking: "))

    # data/seed_nodes/20_data_drug_schizophrenia.txt
    seed_file = open(input("Enter file path for seed nodes: "))
    seeds = {line.strip() for line in seed_file}
    outfile = open(input("Enter file path for output file: "), 'w')
    for line in ranking_file:
        gene, score = line.split()
        outfile.write("{0} {1} {2}\n".format(gene, score, 1 if gene in seeds else 0))
        
    ranking_file.close()
    seed_file.close()
    outfile.close()

def main():
    if int(input("Enter 1 for single-file label, 0 for cross-validation: ")):
        make_labels_for_file()
        return

    # cross_validation/validation_rankings/
    directory_path = input("Enter directory path for input rankings: ")

    # cross_validation/validation_output_labels/
    outfolder_path = input("Enter output folder path: ")
    reverse_flag = int(input("Enter 1 to read file in reverse, 0 otherwise: "))

    infile_path = input("Enter path to partitionings: ")
    num_folds = input("Enter number of folds: ")
    disease = input("Enter disease: ")

    # Does not iterate through directory in order
    for fileName in os.listdir(directory_path):
        i = int(fileName[-5]) if not fileName[-6].isdigit() else int(fileName[-6])*10 + int(fileName[-5])

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
            infile_nonseeds = open(f"{infile_path}/{num_folds}_{disease}_non_seeds_{i}.txt")
            infile_seeds = open(f"{infile_path}/{num_folds}_{disease}_new_seeds_{i}.txt")
            outfile = open(f"{outfolder_path}/{disease}_validation_labels_{i}", 'w')

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

if __name__ == "__main__":
    main()