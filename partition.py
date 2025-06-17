from random import shuffle

def partition(path, i, nodes, disease, bipartition=True):
    shuffle(nodes)
    n = len(nodes)
    nonseeds = nodes[:n//2] if bipartition else nodes[:n//3]
    seeds = nodes[n//2:] if bipartition else nodes[n//3:]
    
    prefix = "2_" if bipartition else "3_"
    with open(path + "/" + prefix + disease + "_new_seeds_" + str(i) + ".txt", 'w') as f1:
        for node in seeds:
            f1.write(node + '\n')
    f1.close()

    with open(path + "/" + prefix + disease + "_non_seeds_" + str(i) + ".txt", 'w') as f2:
        for node in nonseeds:
            f2.write(node + '\n')
    f2.close()

def main():
    two_or_three = int(input("Enter 1 for 2-partitioning, 0 for 3-partitioning: "))
    disease = input("Enter disease: ")

    # cross_validation/schizophrenia
    path = input("Enter folder path to save output: ")

    # data/20_data_drug_schizophrenia.txt
    fileName = input("Enter seed file name: ")
    f = open(fileName, 'r')
    nodes = [ line[:-1] for line in f ]

    for i in range(int(input("Enter number of partitions: "))):
        partition(path, i, nodes, disease, two_or_three)

    print("All done!")
    
if __name__ == "__main__":
    main()