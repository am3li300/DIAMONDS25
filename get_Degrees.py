import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

# data/networks/STRING_protein_links_parsed.tsv
edges_list = input("Enter path to edges: ")
edges = []
file = open(edges_list)

for line in file:
    g1, g2, score = line.split()
    edges.append((g1, g2, score))

graph = nx.Graph()
graph.add_weighted_edges_from(edges)
"""
# output/SZ/main_program_SZ (Run 1).txt
ranked_list = input("Enter path to rankings: ")
outpath = input("Enter path for output file: ")
outfile = open(outpath, "w")
rankings = open(ranked_list, "r")

for line in rankings:
    gene, score = line.split()
    outfile.write(gene + '\t' + str(graph.degree[gene]) + '\n')

rankings.close()
outfile.close()
"""

degrees = [degree for node, degree in graph.degree() if degree <= 400]
print(len([degree for node, degree in graph.degree() if degree > 1000]))
plt.hist(degrees)
plt.xlabel("Degree")
plt.ylabel("Frequency")
plt.show()
