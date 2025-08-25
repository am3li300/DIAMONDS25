import networkx as nx
import markov_clustering as mc
from scipy.sparse import csr_matrix
import numpy as np


############################## Internal Helpers ################################

def build_steiner_tree(full_graph, disease_genes):
        w = nx.get_edge_attributes(full_graph, 'weight')
        max_w = max(w.values())
        nx.set_edge_attributes(full_graph, {e: max_w - wt + 1 for e, wt in w.items()}, name='cost')
        
        idx = {n:i for i,n in enumerate(full_graph.nodes())}
        edges = np.fromiter((i for uv in full_graph.edges()
                        for i in (idx[uv[0]], idx[uv[1]])),
                        dtype=np.int32).reshape(-1,2)

        costs = np.fromiter((full_graph[u][v]['cost'] for u,v in full_graph.edges()), dtype=np.float64)
        prizes = np.zeros(len(idx), dtype=np.float64)
        prizes[[idx[g] for g in disease_genes]] = 1e6      # force inclusion

        
        nodes_kept, edges_kept = pcst_fast.pcst_fast(
                edges,      # numpy int32 [:,2]
                prizes,     # numpy float64 [:]
                costs,      # numpy float64 [:]
                -1,               # root index (-1 = unrooted)
                1,                # num_clusters
                'strong',                # pruning: 0=none, 1=simple, 2=strong
                0)                # verbosity (0 = silent)

        edges = edges.astype(np.int64, copy=False)
        edge_pairs = edges[edges_kept] 
        inv_idx = {i: n for n, i in idx.items()}

        steiner_edges = [(inv_idx[u], inv_idx[v]) for u, v in edge_pairs]

        return full_graph.edge_subgraph(steiner_edges).copy()

def markov_clustering(graph, node_list):
        matrix_array = nx.to_scipy_sparse_array(graph, nodelist=node_list)
        matrix = csr_matrix(matrix_array)

        # inflation = 2; inflation determines granularity
        # allergy: 1.3; SZ & RA: 1.2
        res = mc.run_mcl(matrix, inflation=1.2) 
        return mc.get_clusters(res)


################################### Step 1 #####################################

def cluster_disease_genes(full_graph, disease_genes):
        steiner = build_steiner_tree(full_graph, disease_genes)
        assert disease_genes <= set(steiner.nodes), "Some disease genes not in the steiner tree!"

        # add edges between disease genes
        for gene in disease_genes:
                for neighbor in full_graph.neighbors(gene):
                        if neighbor in disease_genes:
                                steiner.add_edge(gene, neighbor, weight=full_graph[gene][neighbor]['weight'])

        node_list = list(steiner.nodes())
        clustering = markov_clustering(steiner, node_list)

        disease_clusters = []
        for c in clustering:
                disease_clusters.append({node_list[indx] for indx in c})

        return disease_clusters

################################## Step 2 ######################################

def merge_cluster_rankings(rankings, disease_genes):
        def assign_label(score, threshold):
                return score / max_score if score >= threshold else score

        def get_threshold(ranking, min_n=20, fallback_q=75, visualize=False):
                # valley-finding heuristic

                # Parameters
                # ----------
                # ranking : list[(gene, score)]
                # min_n   : int        # require at least this many scores to attempt valley search
                # fallback_q : float   # percentile to return when valley search fails
                
                # --- 0. collect scores and exit early if we have nothing ---
                scores = np.asarray([s for gene, s in ranking if 0 < s <= 1.0 and gene not in disease_genes], dtype=float)
                scores = scores[np.isfinite(scores)]
                if scores.size == 0:
                        return 0

                # --- 1. skip valley-finding on tiny clusters; use quantile fallback ---
                if scores.size < min_n:
                        return float(np.percentile(scores, fallback_q))

                # --- 2. histogram & smoothing (unchanged except for adaptive sigma) ---
                counts, bin_edges = np.histogram(scores, bins=20)
                log_counts = np.where(counts > 0, np.log10(counts), np.nan)

                sigma = max(1, log_counts.size // 50)           # light adaptive smoothing
                smoothed = gaussian_filter1d(log_counts, sigma=sigma)

                minima = argrelextrema(smoothed, np.less)[0]
                maxima = argrelextrema(smoothed, np.greater)[0]

                if visualize:
                        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
                        plt.figure(figsize=(8, 5))
                        plt.bar(bin_centers, log_counts, width=np.diff(bin_edges), alpha=0.4, label="Histogram")
                        plt.plot(bin_centers, smoothed, label="Smoothed")
                        plt.scatter(bin_centers[minima], smoothed[minima], color='red', label='Minima', zorder=3)
                        plt.scatter(bin_centers[maxima], smoothed[maxima], color='green', label='Maxima', zorder=3)
                        plt.xlabel("Score")
                        plt.ylabel("Frequency (log scale)")
                        plt.title("Histogram + Smoothed Curve (for valley detection)")
                        plt.legend()
                        plt.show()

                # --- 3. find valley with widest *score-space* span ---
                widest_span = 0.0
                best_min = None
                for m in minima:
                        left = maxima[maxima < m]
                        right = maxima[maxima > m]
                        if left.size == 0 or right.size == 0:
                                continue

                        span = bin_edges[right.min() + 1] - bin_edges[left.max()]
                        if span > widest_span:
                                widest_span = span
                                best_min = m

                # --- 4. return threshold or fallback ---
                if best_min is not None:
                        return 0.5 * (bin_edges[best_min] + bin_edges[best_min + 1])

                return float(np.percentile(scores, fallback_q))

        
        thresholds = [get_threshold(ranking, 20, 75, False) for ranking in rankings]

        final_scores = {}
        for i, ranking in enumerate(rankings):
                if not ranking:
                        print("ranking {0} is empty".format(i))
                        continue

                max_score = ranking[0][1]
                # print("Threshold for cluster {0}:".format(i), thresholds[i])
                for gene, score in ranking:
                        final_scores[gene] = final_scores.get(gene, 0) + assign_label(score, thresholds[i])

        return sorted(list(final_scores.items()), key=lambda x: -x[1])

################################### Step 3 #####################################

def double_merge(og_ranking, cluster_ranking):
        assert len(og_ranking) == len(cluster_ranking), "Rankings are not the same length"
        
        genes = set()
        final_ranking = []

        def add_gene(curr_line, curr_gene):
                if curr_gene not in genes:
                        final_ranking.append(curr_line)
                        genes.add(curr_gene)


        for i in range(0, len(og_ranking)):
                add_gene(og_ranking[i], og_ranking[i][0])
                add_gene(cluster_ranking[i], cluster_ranking[i][0])
                
        return final_ranking

