from random import shuffle

def write_to_file(file_path, nodes):
    with open(file_path, 'w') as f:
        for node in nodes:
            f.write(node + '\n')

    f.close()

def clustering_partition(clusters, gene_groups, path, i, disease, bipartition=True):
    all_seeds = []
    all_nonseeds = []
    for j in range(len(clusters)):
        g = gene_groups[j][:]
        shuffle(g)

        n = len(g)
        all_nonseeds.extend(g[:n//2] if bipartition else g[:n//3])
        all_seeds.extend(g[n//2:] if bipartition else g[n//3:])

    prefix = "2_" if bipartition else "3_"
    write_to_file(path + "/" + prefix + disease + "_new_seeds_" + str(i) + ".txt", all_seeds)
    write_to_file(path + "/" + prefix + disease + "_non_seeds_" + str(i) + ".txt", all_nonseeds)

def partition(path, i, nodes, disease, bipartition=True):
    shuffle(nodes)
    n = len(nodes)
    nonseeds = nodes[:n//2] if bipartition else nodes[:n//3]
    seeds = nodes[n//2:] if bipartition else nodes[n//3:]
    
    prefix = "2_" if bipartition else "3_"
    write_to_file(path + "/" + prefix + disease + "_new_seeds_" + str(i) + ".txt", seeds)
    write_to_file(path + "/" + prefix + disease + "_non_seeds_" + str(i) + ".txt", nonseeds)

def p(clusters=None, gene_groups=None):
    two_or_three = int(input("Enter 1 for 2-partitioning, 0 for 3-partitioning: "))
    disease = input("Enter disease: ")

    # cross_validation/partitioning/
    # SZ_STRING_louvain_clustering
    path = input("Enter folder path to save output: ")

    # data/seed_nodes/20_data_drug_schizophrenia.txt
    fileName = input("Enter seed file name: ")
    f = open(fileName, 'r')
    nodes = [ line.strip() for line in f ]

    offset = int(input("Enter offset for file number: "))
    
    if clusters is None or gene_groups is None:
        for i in range(int(input("Enter number of partitions: "))):
            partition(path, i+offset, nodes, disease, two_or_three)

    else:
        for i in range(int(input("Enter number of clustering partitions: "))):
            clustering_partition(clusters, gene_groups, path, i+offset, disease, two_or_three)

    print("All done!")

""" 
if __name__ == "__main__":
    p()
"""