import networkx as nx

edges_list = input("Enter path to edges: ")
edges = []
file = open(edges_list)

for line in file:
    g1, g2, score = line.split()
    edges.append((g1, g2, score))

graph = nx.Graph
graph.add_weighted_edges_from(edges)

ranked_list = input("Enter path to rankings: ")

outpath = input("Enter path for output file: ")

outfile = open(outpath, "w")
outfile.write()