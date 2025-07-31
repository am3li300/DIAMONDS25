def main():
    # data/seed_nodes/
    seeds_file = open(input("Enter seed file path: "))
    seeds = {gene.strip() for gene in seeds_file}

    # output/
    ranking_file = open(input("Enter ranking file path: "))

    k = int(input("Enter top k value: "))
    assert k > 0, "k must be greater than 0"

    top_k_output = open(input("Enter output file path: "), 'w')

    line_number = 0
    top_k_output.write("Gene Score Original_Placement\n")

    for line in ranking_file:
        gene, score = line.split()
        if gene not in seeds:
            top_k_output.write("{0} {1} {2}\n".format(gene, score, line_number))
            k -= 1
            if not k:
                break

        line_number += 1

if __name__ == "__main__":
    main()