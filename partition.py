from random import shuffle

def partition(path, i, nodes, disease):
    shuffle(nodes)
    n = len(nodes)
    nonseeds = nodes[:n//2]
    seeds = nodes[n//2:]
    
    with open(path + "/" + disease + "_new_seeds_" + str(i) + ".txt", 'w') as f1:
        for node in seeds:
            f1.write(node + '\n')
    f1.close()

    with open(path + "/" + disease + "_non_seeds_" + str(i) + ".txt", 'w') as f2:
        for node in nonseeds:
            f2.write(node + '\n')
    f2.close()

def main():
    disease = input("Enter disease: ")
    path = input("Enter folder path to save output: ")
    fileName = input("Enter seed file name: ")
    f = open(fileName, 'r')
    nodes = [ line[:-1] for line in f ]
    for i in range(25):
        partition(path, i, nodes)

    print("All done!")
    
if __name__ == "__main__":
    main()