import pandas as pd
import networkx as nx

df = pd.read_csv("9606__protein_links.tsv", sep=' ')
df = df[(df["combined_score"] >= 700) & ((df["experiments"] >= 1) | (df["database"] >= 1))]
df = df[["protein1", "protein2", "combined_score"]]

graph = nx.from_pandas_edgelist(df = df, source = 'protein1', target = 'protein2')
component = list(nx.connected_components(graph))[0]
subgraph = graph.subgraph(component)
print("Comp edges:", subgraph.number_of_edges())

print(graph.number_of_edges())
"""
df.to_csv('9606__protein_links_parsed.txt', sep=' ', index=False)
print("Saved dataframe to 9606__protein_links_parsed.txt")
"""
