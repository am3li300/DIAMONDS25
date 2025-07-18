from t_map.garbanzo.edgelist import EdgeListGarbanzo
from t_map.feta.glide import ADAGIO
from t_map.gene.gene import Gene

import argparse
from time import time
from heapq import heapify, heappush, heappop

import networkx as nx
from scipy.sparse import csr_matrix
import markov_clustering as mc
from local_community_detection import greedy_source_expansion

import igraph as ig
import multiprocessing

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from partition import *

import pcst_fast
import numpy as np
from joblib import Parallel, delayed

import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.signal import argrelextrema


parser = argparse.ArgumentParser()
parser.add_argument('--network', '-n', type=str, required=True, help="Path to edgelist, assumes tab separated.")
parser.add_argument('--genelist', '-g', type=str, required=True, help="Path to genelist, assumes genes are newline separated and in the network.")
parser.add_argument('--out', '-o', type=str, default="adagio.out",help="Path to output results")


def get_disease_genes(genelist_path):
        disease_genes = []
        with open(genelist_path) as f:
            for line in f.readlines():
                disease_genes.append(line.strip())

        return disease_genes

############################ Clustering ###############################

def merge_cluster_rankings(rankings):
        heap = []
        for i, ranking in enumerate(rankings):
                if ranking:
                        heap.append((-ranking[0][1], ranking[0][0].name, i, 0))
 
        heapify(heap)
        res = []
        while heap:
                score, gene, indx, jndx = heappop(heap)
                res.append((Gene(name=gene), -score))
                if jndx < len(rankings[indx])-1:
                        heappush(heap, (-rankings[indx][jndx+1][1], rankings[indx][jndx+1][0].name, indx, jndx+1))

        return res
         
def get_cluster_rankings(clusters, full_graph, disease_genes):
        n = len(clusters)
        gene_groups = [[] for _ in range(n)] # store the seed nodes for each cluster

        # store seed nodes as gene objects in a list from genelist_path
        for gene in disease_genes:
                for i in range(n):
                        if gene in clusters[i]:
                                gene_groups[i].append(Gene(name=gene))

        """
        # generate cross validation partitions
        p(clusters, [[g.name for g in genes] for genes in gene_groups])
        return []
        """
        # generate rankings for each cluster
        rankings = []
        for i, cluster in enumerate(clusters):
                if (not gene_groups[i]):
                        continue

                sub_graph = full_graph.subgraph(cluster).copy()

                # only keep seeds with at least one edge (not isolated)
                seeds = [g for g in gene_groups[i] if sub_graph.degree(g.name) > 0]
                if not seeds:
                        continue

                model = ADAGIO()
                model.setup(sub_graph)
                model.set_add_edges_amount(20)
                rankings.append(sorted(list(model.prioritize(seeds, sub_graph)), key=lambda x: -x[1]))

        return rankings

def clustering(network_path, genelist_path, algorithm="louvain"):
        full_graph = nx.read_weighted_edgelist(network_path)
        disease_genes = get_disease_genes(genelist_path)

        if algorithm == "louvain":
                clusters = nx.community.louvain_communities(full_graph)

        elif algorithm == "markov":
                indices_to_nodes = list(full_graph.nodes())

                matrix_array = nx.to_scipy_sparse_array(full_graph, nodelist=indices_to_nodes)
                matrix = csr_matrix(matrix_array)

                # expansion = 2, inflation = 2
                res = mc.run_mcl(matrix)

                clustering = mc.get_clusters(res)
                clusters = []
                for c in clustering:
                        clusters.append({indices_to_nodes[indx] for indx in c})

        elif algorithm == "walktrap":
                indices_to_nodes = list(full_graph.nodes())
                nodes_to_indices = {node: i for i, node in enumerate(indices_to_nodes)}

                edges = [(nodes_to_indices[u], nodes_to_indices[v]) for u, v in full_graph.edges()]
                weights = [full_graph[u][v]["weight"] for u, v in full_graph.edges()]

                iGraph = ig.Graph(edges=edges, directed=False)
                iGraph.vs["name"] = indices_to_nodes
                iGraph.es["weight"] = weights

                # steps = 3-5
                wtrap = iGraph.community_walktrap(steps=4, weights="weight")
                clustering = list(wtrap.as_clustering())
                clusters = []

                for c in clustering:
                        clusters.append({indices_to_nodes[indx] for indx in c})

        else:
                return []

        disease_genes = get_disease_genes(genelist_path)
        rankings = get_cluster_rankings(clusters, full_graph, disease_genes)
        return merge_cluster_rankings(rankings)


####################### Supervised Clustering ###############################

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

def run_adagio(full_graph, disease_genes, disease_clusters):
        set_up_start = time()
        model_all = ADAGIO()
        model_all.setup(full_graph) # change to steiner for quick testing
        model_all.set_add_edges_amount(20)

        print("Set up time:", time()-set_up_start)

        # score clusters in parallel
        def score_cluster(cluster):
                s = time()
                seeds = [Gene(name=g) for g in cluster if g in disease_genes]
                if not seeds:
                        return []
                print("Number of seeds:", len(seeds))
                cluster_ranking = sorted(list(model_all.prioritize(seeds, full_graph)), key=lambda x:-x[1])
                print("Time of one run:", time()-s)
                return cluster_ranking

        return Parallel(n_jobs=2)(delayed(score_cluster)(cl) for cl in disease_clusters)

def merge_supervised_cluster_rankings(rankings):
        max_score = max(ranking[0][1] for ranking in rankings)

        def assign_label(score, threshold):
                return score / max_score if score >= threshold else score

        def get_threshold(ranking, min_n=20, fallback_q=75, visualize=True):
                """
                Minimal but safer rewrite of the original valley-finding heuristic.

                Parameters
                ----------
                ranking : list[(gene, score)]
                min_n   : int        # require at least this many scores to attempt valley search
                fallback_q : float   # percentile to return when valley search fails
                """
                # --- 0. collect scores and exit early if we have nothing ---
                scores = np.asarray([s for _, s in ranking if s <= 1.0], dtype=float)
                scores = scores[np.isfinite(scores)]
                if scores.size == 0:
                        return 0

                # --- 1. skip valley-finding on tiny clusters; use quantile fallback ---
                if scores.size < min_n:
                        return float(np.percentile(scores, fallback_q))

                # --- 2. histogram & smoothing (unchanged except for adaptive sigma) ---
                counts, bin_edges = np.histogram(scores, bins="auto")
                if counts.size < 3:                         # need at least 3 bins
                        return float(np.percentile(scores, fallback_q))

                sigma = max(1, counts.size // 50)           # light adaptive smoothing
                smoothed = gaussian_filter1d(counts, sigma=sigma)

                minima = argrelextrema(smoothed, np.less)[0]
                maxima = argrelextrema(smoothed, np.greater)[0]

                if visualize:
                        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
                        plt.figure(figsize=(8, 5))
                        plt.bar(bin_centers, counts, width=np.diff(bin_edges), alpha=0.4, label="Histogram")
                        plt.plot(bin_centers, smoothed, label="Smoothed")
                        plt.scatter(bin_centers[minima], smoothed[minima], color='red', label='Minima', zorder=3)
                        plt.scatter(bin_centers[maxima], smoothed[maxima], color='green', label='Maxima', zorder=3)
                        plt.xlabel("Score")
                        plt.ylabel("Frequency")
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

        
        threshold = sum(get_threshold(ranking) for ranking in rankings) / len(rankings)
        print("Threshold:", threshold)

        final_scores = {}
        for ranking in rankings:
                for gene, score in ranking:
                        final_scores[gene] = final_scores.get(gene, 0) + assign_label(score, threshold)  # find general or adaptive threshold for other diseases

        return sorted(list(final_scores.items()), key=lambda x: -x[1])

        """
        max_score = max(ranking[0][1] for ranking in rankings)

        for ranking in rankings:
                for i in range(len(ranking)):
                        ranking[i] = (ranking[i][0], pow(ranking[i][1] / max_score, 4))

        from collections import defaultdict
        max_rankings = defaultdict(int)

        for ranking in rankings:
               for gene, score in ranking:
                        max_rankings[gene] += score 
        
        return sorted(list(max_rankings.items()), key=lambda x: -x[1])
        """


def supervised_clustering(network_path, genelist_path):
        disease_genes = set(get_disease_genes(genelist_path))
        full_graph = nx.read_weighted_edgelist(network_path)

       
        steiner = build_steiner_tree(full_graph, disease_genes)
        assert disease_genes <= set(steiner.nodes), "Some disease genes not in the steiner tree!"
        # print(nx.number_connected_components(steiner), len(steiner.edges), len(steiner.nodes))

        # add edges between disease genes
        for gene in disease_genes:
                for neighbor in full_graph.neighbors(gene):
                        if neighbor in disease_genes:
                                steiner.add_edge(gene, neighbor, weight=full_graph[gene][neighbor]['weight'])
                                        
        disease_clusters = nx.community.louvain_communities(steiner, resolution=0.2) # resolution determines num of clusters
        print("Number of clusters:", len(disease_clusters))

        rankings = run_adagio(full_graph, disease_genes, disease_clusters) # change to steiner for quick testing
        """
        disease = input("Enter disease: ")
        for i, ranking in enumerate(rankings):
                with open(f"../output/disease_cluster/{disease}_cluster_{i}.txt", "w") as fout:
                        for gene, score in ranking:
                                fout.write(f"{gene.name}\t{score}\n")

        """
        return merge_supervised_cluster_rankings(rankings)
        

def main(network_path: str, genelist_path: str, out_path: str="adagio.out"):
        start = time()

        """
        Baseline ADAGIO - constant k for disease nodes only (k = 20)
        """
        # graph = EdgeListGarbanzo(network_path, genelist_path)
        # print(len(graph.graph.edges))
        # model = ADAGIO()
        # model.setup(graph.graph)
        # model.set_add_edges_amount(20) # this will add edges to the graph
        # sorted(list(model.prioritize(graph.genes, graph.graph)), key=lambda x: x[1], reverse=True)

        """
        adaptive k for disease nodes only
        """
        # predictions = sorted(list(model.david_prioritize_2(graph.genes, graph.graph)), key=lambda x: x[1], reverse=True)

        """
        adaptive k for all nodes
        """
        # predictions = sorted(list(model.david_prioritize(1000, graph.graph, True)), key=lambda x: x[1], reverse=True)

        """
        constant k for all nodes (k = 20)
        """
        # graph = EdgeListGarbanzo(network_path, genelist_path)
        # print(len(graph.graph.edges))
        # model = ADAGIO()
        # model.setup(graph.graph)
        # model.set_add_edges_amount(20) # this will add edges to the graph
        # predictions = sorted(list(model.prioritize(graph.genes, graph.graph)), key=lambda x: x[1], reverse=True)


        """
        unsupervised clustering - louvain, markov, walktrap
        """
        # predictions = clustering(network_path, genelist_path, algorithm="louvain")

        """ 
        supervised clustering
        """
        predictions = supervised_clustering(network_path, genelist_path)


        with open(out_path, "w") as f:
                for gene, score in predictions:
                        f.write(f"{gene.name}\t{score}\n")

        end = time()
        print("Total time:", end-start)
        

if __name__ == "__main__":
        args = parser.parse_args()
        main(args.network, args.genelist, args.out)
        




"""
def gse_helper(graph, gene):
        return greedy_source_expansion(graph, source=gene)

def graveyard():
        def jaccard(set1, set2):
                union = set1.union(set2)
                intersection = set1.intersection(set2)
                return float(len(intersection))/len(union) if union else 0

        full_graph = nx.read_weighted_edgelist(network_path)
        disease_genes = get_disease_genes(genelist_path)
   
        gse_args = [(full_graph, gene) for gene in disease_genes]

        with multiprocessing.Pool() as pool:
                clusters = pool.starmap(gse_helper, gse_args)

        # get jaccard matrix of clusters
        n = len(clusters)
        jmat = [[0]*n for _ in range(n)]
        for i in range(n - 1):
                jmat[i][i] = 1
                for j in range(i + 1, n):
                        jmat[i][j] = jmat[j][i] = jaccard(clusters[i], clusters[j])

        for i in range(n):
                jmat[n - 1][i] = jmat[i][n - 1]
        
        # merge clusters based on jaccard scores

        rankings = get_cluster_rankings(clusters, full_graph, disease_genes)
        ranking_dict = {}

        for ranking in rankings:
                for gene, score in ranking:
                        ranking_dict[gene.name] = max(ranking_dict.get(gene.name, 0), score)
        
        return sorted([[Gene(x), ranking_dict[x]] for x in ranking_dict], key=lambda x: -x[1])


        rankings = []
        total = time()
        for cluster in disease_clusters:
                s = time()
                seeds = [Gene(name=g) for g in cluster if g in disease_genes]
                if not seeds:
                        continue

                res = list(model_all.prioritize(seeds, full_graph))
                print("Time to finish run:", time()-s)
                rankings.append(res)
        
        print("Time for all runs:", time()-total)

def gse_helper(graph, gene):
        return greedy_source_expansion(graph, source=gene)
"""