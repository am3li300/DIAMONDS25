import pandas as pd
import numpy as np
import networkx as nx
from scipy.spatial.distance import pdist, squareform

def cw(A, **kwargs):
    assert isinstance(A, np.ndarray) and (len(A.shape) == 2) and (A.shape[0] == A.shape[1])
    weighted   = kwargs.get("weighted", False)
    normalized = kwargs.get("normalize", True)
    n, _ = A.shape
    if not weighted:
        A = np.where(A > 0, 1, 0)
    A1 = np.where(A > 0, 1, 0)
    C  = A @ A1
    C  = C + C.T
    if normalized:
        d1_2 = np.sqrt(A @ np.ones((n, 1)))
        C    = (C / d1_2) / d1_2.T
    if "rescale" in kwargs:
        C = C / np.max(C)
    return C

def compute_dsd_embedding(A, t=-1, gamma=1, is_normalized=True):
    n, _ = A.shape
    e = np.ones((n, 1))
    d = A @ e
    P = A/d
    W = (1 / np.sum(d)) * np.dot(e, (e * d).T)
    P1 = gamma * (P - W)
    X = np.linalg.pinv(np.eye(n, n) - P1)
    if t > 0:
        P1_t = np.eye(n, n) - np.linalg.matrix_power(P1, t)
        X = np.matmul(X, P1_t)
    if not is_normalized:
        return X
    return (X / np.sqrt(d).T)

def glide(A, alpha=0.1, beta=1000, delta=1, gamma=1, normalize_dsd=False, **kwargs):
    assert (A is not None) and (len(A.shape) == 2) and (A.shape[0] == A.shape[1])
    n, _ = A.shape
    L = None
    Ge = None
    Gd = None
    if ("G_emb" in kwargs) or ("L_sc" in kwargs):
        if "G_emb" in kwargs:
            Ge = kwargs["G_emb"]
        if "L_sc" in kwargs:
            L = kwargs["L_sc"]
    elif ("localf" in kwargs) or ("globalf" in kwargs):
        if "localf" in kwargs:
            local = kwargs["localf"]
            L = local(A, **kwargs)
        if "globalf" in kwargs:
            global_ = kwargs["globalf"]
            Ge = global_(A, **kwargs)
    L = cw(A, **kwargs)
    if Ge is None:
        Ge = compute_dsd_embedding(A, gamma=gamma, is_normalized=normalize_dsd)
    Gd = squareform(pdist(Ge)) + np.eye(n, n)
    R = np.exp(alpha / (1 + beta * Gd)) * L + (delta / Gd)
    return R

def read_ppi_network(file_path):
    df = pd.read_csv(file_path, sep=" ")
    return df

def create_ppi_network(ppi_edges):
    g = nx.from_pandas_edgelist(ppi_edges, source="protein1", target="protein2", edge_attr="combined_score")
    return g

def reweight_graph_with_glide(G, nodes_to_expand=None, alpha=0.1, beta=1000, delta=1, gamma=1, normalize_dsd=False):
    node_list = list(G.nodes())
    node_to_index = {node: i for i, node in enumerate(node_list)}
    index_to_node = {i: node for node, i in node_to_index.items()}
    A = nx.to_scipy_sparse_array(G, weight='combined_score')
    A = A.toarray().astype(float)
    print(f"Adjacency matrix A: shape {A.shape}, nonzeros {np.count_nonzero(A)}, min {A.min()}, max {A.max()}")
    glide_weights = glide(A, alpha=alpha, beta=beta, delta=delta, gamma=gamma, normalize_dsd=normalize_dsd)
    print(f"GLIDE weights: min {glide_weights.min()}, max {glide_weights.max()}, sum {glide_weights.sum()}")
    existing_edges_mask = A > 0
    new_weights = glide_weights * existing_edges_mask
    print(f"New weights: min {new_weights.min()}, max {new_weights.max()}, sum {new_weights.sum()}, nonzero {np.count_nonzero(new_weights)}")
    for i, j in zip(*existing_edges_mask.nonzero()):
        if i < j:
            node1 = index_to_node[i]
            node2 = index_to_node[j]
            G[node1][node2]['weight'] = float(new_weights[i, j])
    if nodes_to_expand:
        for node in nodes_to_expand:
            if node in G:
                node_idx = node_to_index[node]
                original_degree = G.degree(node)
                potential_edges = [(node, index_to_node[j], glide_weights[node_idx, j])
                                   for j in range(len(G))
                                   if j != node_idx and not G.has_edge(node, index_to_node[j])]
                top_edges = sorted(potential_edges, key=lambda x: x[2], reverse=True)[:original_degree]
                for source, target, weight in top_edges:
                    G.add_edge(source, target, weight=float(weight))
    total_weight = sum(d['weight'] for u, v, d in G.edges(data=True))
    print(f"Total edge weight after reweighting: {total_weight}")
    if total_weight == 0:
        raise RuntimeError("Your reweighted graph has no edge weights. Debug your GLIDE output.")
    return G

def extract_disease_genes(G, alias_file="9606.protein.aliases.v12.0.txt", tsv_file="schizophrenia-DOID_5419-genes-2025-06-07.tsv"):
    df = pd.read_csv(tsv_file, sep='\t')
    human_genes = set(df["Gene Symbol"].dropna().unique())
    print(f"Number of unique human schizophrenia genes: {len(human_genes)}")
    alias_df = pd.read_csv(alias_file, sep="\t", header=None, names=["protein_id", "alias", "source"])
    mapping = alias_df[alias_df["alias"].isin(human_genes)][["alias", "protein_id"]].drop_duplicates()
    seed_protein_ids = mapping["protein_id"].unique().tolist()
    print(f"Found {len(seed_protein_ids)} seed STRING IDs for {len(human_genes)} gene symbols.")
    nodes_to_expand = [node for node in seed_protein_ids if node in G]
    if not nodes_to_expand:
        raise ValueError("None of the mapped seed proteins are present in the PPI network!")
    print("Using nodes_to_expand:", nodes_to_expand[:10], "..." if len(nodes_to_expand) > 10 else "")

    def save_nodes_to_txt():
        with open("nodes_to_expand.txt", 'w') as file:
            for node in nodes_to_expand:
                file.write(node + '\n')

        print("Saved nodes to expand results to nodes_to_expand.txt")

    #save_nodes_to_txt()
    return nodes_to_expand

# --- MAIN PIPELINE ---

# 1. Load full PPI
ppi_df = read_ppi_network("9606.protein.physical.links.v12.0.txt")
ppi_df = ppi_df[ppi_df["combined_score"] >= 600]
G_full = create_ppi_network(ppi_df)
print(f"Graph: {G_full.number_of_nodes()} nodes, {G_full.number_of_edges()} edges")

# 2. Extract mapped schizophrenia protein IDs as seeds
# nodes_to_expand = extract_disease_genes(G_full)
nodes_to_expand = ["9606.ENSP00000403888", "9606.ENSP00000354511", "9606.ENSP00000287842", "9606.ENSP00000332549", "9606.ENSP00000354859", "9606.ENSP00000303252", "9606.ENSP00000341680", "9606.ENSP00000451828", "9606.ENSP00000392423", "9606.ENSP00000380878", "9606.ENSP00000342235", "9606.ENSP00000381382", "9606.ENSP00000489407"]
# other relevant schizophrenia genes: N/A

"""
# 3. Build 1-hop subgraph around seed proteins
neighbor_nodes = set(nodes_to_expand)
for node in nodes_to_expand:
    if node in G_full:
        neighbor_nodes.update(G_full.neighbors(node))
G_sub = G_full.subgraph(neighbor_nodes).copy()
print(f"Subgraph: {G_sub.number_of_nodes()} nodes, {G_sub.number_of_edges()} edges")
"""

# 4. Run reweighting and ranking pipeline
G_reweighted = reweight_graph_with_glide(G_full, nodes_to_expand=nodes_to_expand)

for node in nodes_to_expand:
    print(f"{node} degree: {G_reweighted.degree(node)}")
    for neighbor in list(G_reweighted.neighbors(node))[:5]:
        print(f"  edge {node}-{neighbor} weight: {G_reweighted[node][neighbor]['weight']}")

pr = nx.pagerank(G_reweighted, personalization={node: 1 for node in nodes_to_expand})
pr_df = pd.Series(pr, name="score").sort_values(ascending=False)
pr_df.to_csv("pagerank_results.tsv", sep="\t", header=True)
print("Saved PageRank results to pagerank_results.tsv")
