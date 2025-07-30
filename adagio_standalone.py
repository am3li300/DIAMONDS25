import pandas as pd
import numpy as np
import networkx as nx
from scipy.spatial.distance import pdist, squareform
from time import time

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
    """
    Extracts schizophrenia-associated seed proteins present in the input PPI graph.

    This function performs the following steps:
    1. Reads a TSV file of disease-associated genes (provided by a database).
    2. Extracts the unique set of gene symbols associated with the disease.
    3. Loads a STRING alias mapping file, which maps gene symbols and other aliases to STRING protein identifiers.
    4. Filters the alias mappings to retain only those that match the disease-associated gene symbols.
    5. From the filtered mappings, it selects the unique STRING protein IDs.
    6. It checks which of these protein IDs are actually present in the provided PPI network `G`.
    7. Returns the list of STRING protein IDs that both match disease-associated genes and exist in the graph.

    The function also prints out:
    - The number of unique disease gene symbols.
    - The number of successfully mapped STRING IDs.
    - A preview of the final list of seed nodes (up to 10 shown).

    Parameters:
    ----------
    G : networkx.Graph
        The protein-protein interaction graph where STRING protein IDs are the node identifiers.

    alias_file : str
        Path to the STRING alias mapping file (e.g., "9606.protein.aliases.v12.0.txt").
        This file should contain mappings between STRING protein IDs and gene symbols/aliases.

    tsv_file : str
        Path to a TSV file containing disease-associated genes (e.g., schizophrenia).
        This file should include at least a column named "Gene Symbol".

    Returns:
    -------
    nodes_to_expand : List[str]
        A list of STRING protein IDs that are both associated with the disease of interest and
        are present in the given PPI network `G`. These serve as seed nodes for the ADAGIO program.
    """

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
# data/STRING_protein_links_parsed.txt
fileName = input("Enter PPI Network file path/name: ")
ppi_df = read_ppi_network(fileName)
G_full = create_ppi_network(ppi_df)
print(f"Graph: {G_full.number_of_nodes()} nodes, {G_full.number_of_edges()} edges")

# 2. Establish seed nodes
choice = int(input("Enter '0' to use alias and tsv file mapping, '1' to use pre-defined seed nodes, '2' to use txt file: "))
if not choice:
    aliasFileName = input("Enter alias file path/name for gene symbol mapping: ")
    tsvFileName = input("Enter tsv file path/name for known disease-associated genes: ")
    nodes_to_expand = extract_disease_genes(G_full, aliasFileName, tsvFileName)

elif choice == 1:
    nodes_to_expand = [node for node in input("Enter space-separated seed nodes: ").split() if node in G_full]
    print(f"Found {len(nodes_to_expand)} genes in the network.")

else:
    # ../PPI Networks/Human/Data/20_data_schizophrenia.txt
    fileName = input("Enter seed file path/name: ")
    file = open(fileName, 'r')
    nodes_to_expand = []
    for line in file:
        node = line[:-1] if line[-1] == '\n' else line
        if node in G_full:
            nodes_to_expand.append(node)

    print(f"Found {len(nodes_to_expand)} genes in the network.")
    
def cross_validation():
    from rand import randint
    indx = randint(0, len(nodes_to_expand)-1)
    print(f"Removing {nodes_to_expand[indx]} from seeds.")
    nodes_to_expand.pop(indx)

#cross_validation()

# 4. Run reweighting and ranking pipeline
start = time()
G_reweighted = reweight_graph_with_glide(G_full, nodes_to_expand=nodes_to_expand)

for node in nodes_to_expand:
    print(f"{node} degree: {G_reweighted.degree(node)}")
    for neighbor in list(G_reweighted.neighbors(node))[:5]:
        print(f"  edge {node}-{neighbor} weight: {G_reweighted[node][neighbor]['weight']}")

pr = nx.pagerank(G_reweighted, personalization={node: 1 for node in nodes_to_expand})
pr_df = pd.Series(pr, name="score").sort_values(ascending=False)
pr_df.to_csv("pagerank_results.tsv", sep="\t", header=True)

end = time()
print("Saved PageRank results to pagerank_results.tsv")
print("'start': {0}\n'end': {1}\n'total time': {2}".format(start, end, end-start))


"""
glide-compute \
  --network path/to/ppi.tsv \
  --output  glide.pkl            # or .tsv if asking for neighbours
  --get-glide-neighbors \
  --glide-neighbors-k 50         # top-50 per node (change as needed)
  --neighbors-return-format dataframe
"""

# output/disease_cluster/schizophrenia_cluster_0.txt
# STRING network
"""
python3 main.py \
  --network '../data/networks/STRING_protein_links_parsed.tsv' \
  --genelist '../data/seed_nodes/schizophrenia_drug.txt' 
"""

# cross validation
"""
python3 main.py \
  --network '../data/networks/STRING_protein_links_parsed.tsv' \
  --genelist '../cross_validation/partitions/schizophrenia_STRING/3_schizophrenia_new_seeds_0.txt' \
  --out '../cross_validation/rankings/SZ_STRING_supervised/3_schizophrenia_cross_validation_0.out'

python3 main.py \
  --network '../data/networks/STRING_protein_links_parsed.tsv' \
  --genelist '../cross_validation/partitions/schizophrenia_STRING/3_schizophrenia_new_seeds_0.txt' \
  --out '../lenore_clustering_test'

python3 main.py \
  --network '../data/networks/STRING_protein_links_parsed.tsv' \
  --genelist '../cross_validation/partitions/diabetes_set_STRING/3_diabetes_new_seeds_0.txt' \
  --out '../DB_test'

"""

# small dataset
"""
python3 main.py \
  --network '../data/networks/STRING_sample_links.tsv' \
  --genelist '../data/seed_nodes/small_dataset.txt' \
  --out d.out 
"""

