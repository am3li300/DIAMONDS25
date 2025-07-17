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
from scipy.signal import find_peaks


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

        def assign_label(score, threshold=0.003):
                return score / max_score if score >= threshold else score

        def get_threshold(ranking):
                """
                Compute a score threshold from a 2D list-like `ranking` whose elements are [Gene(), float_score].
                Steps:
                1. Extract scores < SCORE_CUTOFF.
                2. Histogram + Gaussian smooth.
                3. Find local maxima & minima in smoothed counts.
                4. Build all valleys (max->min->max triplets).
                5. Select widest valley on score axis.
                6. Derive several representative cut values; choose one via heuristic (sharp valley -> argmin, else midpoint).

                Returns
                -------
                dict with keys (some may be None if insufficient data):
                        scores              : np.ndarray of used scores
                        bin_edges            : np.ndarray
                        bin_centers          : np.ndarray
                        score_counts         : np.ndarray (raw histogram)
                        smoothed_counts      : np.ndarray (smoothed histogram)
                        valleys              : list of valley dicts (see below)
                        widest_valley        : valley dict for widest or None
                        cut_value            : chosen numeric threshold or None

                Each valley dict contains:
                        min_bin, left_max_bin, right_max_bin
                        width, depth, flatness
                        min_count, left_max_count, right_max_count
                        left_edge, right_edge
                        min_bin_center, valley_midpoint, argmin_center, quadratic_argmin, inverse_weighted_center
                        cut_value  (valley-level heuristic pick; the overall function returns the widest valley's cut_value)

                Edit the CONFIG constants below to tune behavior.
                """

                # -----------------
                # CONFIG (edit as needed)
                # -----------------
                SCORE_CUTOFF = 1       # only include scores < this
                BINS = 50              # histogram bins
                SIGMA = 2              # Gaussian smoothing sigma (in bins)
                PROMINENCE = None      # pass to find_peaks to filter noise (None = no filter)
                RANGE = None           # (min,max) for histogram; e.g., (0,1). None = auto from data.
                FLOOR_FRAC = 0.05      # fraction of valley depth defining "floor" for flatness calc
                FLATNESS_FRAC = 0.1    # if std(floor) < FLATNESS_FRAC * depth -> sharp valley -> use argmin

                # -----------------
                # Nested helpers
                # -----------------
                def _extract_scores(ranking, cutoff):
                        return np.array([row[1] for row in ranking if row[1] < cutoff], dtype=float)

                def _hist_and_smooth(scores, bins, sigma, range_):
                        counts, edges = np.histogram(scores, bins=bins, range=range_)
                        smoothed = gaussian_filter1d(counts.astype(float), sigma=sigma)
                        centers = 0.5 * (edges[:-1] + edges[1:])
                        return counts, edges, smoothed, centers

                def _flanking_maxima(m, max_idx):
                        left = max_idx[max_idx < m]
                        right = max_idx[max_idx > m]
                        if len(left) == 0 or len(right) == 0:
                                return None, None
                        return left[-1], right[0]

                def _valley_stats(m, left_max, right_max, edges, smoothed, centers):
                        left_edge = edges[left_max]
                        right_edge = edges[right_max + 1]  # upper bound of right max bin
                        width = right_edge - left_edge

                        lc = smoothed[left_max]
                        rc = smoothed[right_max]
                        mc = smoothed[m]
                        depth = min(lc, rc) - mc  # conservative

                        sl = slice(left_max, right_max + 1)
                        local_counts = smoothed[sl]
                        local_centers = centers[sl]

                        # representatives
                        min_bin_center = 0.5 * (edges[m] + edges[m + 1])
                        valley_midpoint = 0.5 * (left_edge + right_edge)
                        argmin_center = local_centers[np.argmin(local_counts)]

                        # quadratic fit
                        if len(local_centers) >= 3:
                                a, b, c = np.polyfit(local_centers, local_counts, deg=2)
                                quad_min = -b / (2 * a) if a != 0 else argmin_center
                        else:
                                quad_min = argmin_center

                        # inverse-count weighted
                        counts_clip = np.clip(local_counts, 1e-12, None)
                        weights = 1.0 / counts_clip
                        inv_weight_center = np.sum(weights * local_centers) / np.sum(weights)

                        # flatness
                        if depth > 0:
                                low_mask = local_counts < (mc + FLOOR_FRAC * depth)
                        else:
                                low_mask = np.ones_like(local_counts, dtype=bool)

                        low_counts = local_counts[low_mask]
                        flatness = np.std(low_counts) if low_counts.size > 1 else 0.0

                        # valley-level cut heuristic
                        if depth > 0 and flatness < FLATNESS_FRAC * depth:
                                cut_value = argmin_center
                        else:
                                cut_value = valley_midpoint

                        return {
                        "min_bin": int(m),
                        "left_max_bin": int(left_max),
                        "right_max_bin": int(right_max),
                        "width": float(width),
                        "depth": float(depth),
                        "flatness": float(flatness),
                        "min_count": float(mc),
                        "left_max_count": float(lc),
                        "right_max_count": float(rc),
                        "left_edge": float(left_edge),
                        "right_edge": float(right_edge),
                        "min_bin_center": float(min_bin_center),
                        "valley_midpoint": float(valley_midpoint),
                        "argmin_center": float(argmin_center),
                        "quadratic_argmin": float(quad_min),
                        "inverse_weighted_center": float(inv_weight_center),
                        "cut_value": float(cut_value),
                        }

                # -----------------
                # Main body
                # -----------------
                scores = _extract_scores(ranking, SCORE_CUTOFF)
                if scores.size == 0:
                        return {
                        "scores": scores,
                        "cut_value": None,
                        "widest_valley": None,
                        "valleys": [],
                        "bin_edges": None,
                        "bin_centers": None,
                        "score_counts": None,
                        "smoothed_counts": None,
                        }

                score_counts, bin_edges, smoothed_counts, bin_centers = _hist_and_smooth(
                        scores, BINS, SIGMA, RANGE
                )

                # maxima / minima
                max_idx, _ = find_peaks(smoothed_counts, prominence=PROMINENCE)
                min_idx, _ = find_peaks(-smoothed_counts, prominence=PROMINENCE)

                valleys = []
                for m in min_idx:
                        left_max, right_max = _flanking_maxima(m, max_idx)
                        if left_max is None or right_max is None:
                                continue
                        v = _valley_stats(m, left_max, right_max, bin_edges, smoothed_counts, bin_centers)
                        valleys.append(v)

                widest_valley = max(valleys, key=lambda d: d["width"]) if valleys else None
                cut_value = widest_valley["cut_value"] if widest_valley is not None else None

                return {
                        "scores": scores,
                        "bin_edges": bin_edges,
                        "bin_centers": bin_centers,
                        "score_counts": score_counts,
                        "smoothed_counts": smoothed_counts,
                        "valleys": valleys,
                        "widest_valley": widest_valley,
                        "cut_value": cut_value,
                }
        
        threshold = sum(get_threshold(ranking)["cut_value"] for ranking in rankings) / len(rankings)
        print("Threshold:", threshold)

        final_scores = {}
        for ranking in rankings:
                for gene, score in ranking:
                        final_scores[gene] = final_scores.get(gene, 0) + assign_label(score, 0.003)  # find general or adaptive threshold for other diseases

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

        rankings = run_adagio(steiner, disease_genes, disease_clusters) # change to steiner for quick testing
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