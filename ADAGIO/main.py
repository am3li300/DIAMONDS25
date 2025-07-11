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

def merge_rankings(rankings):
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
        return merge_rankings(rankings)


def gse_helper(graph, gene):
        return greedy_source_expansion(graph, source=gene)
        
def supervised_clustering(network_path, genelist_path):
        def jaccard(set1, set2):
                union = set1.union(set2)
                intersection = set1.intersection(set2)
                return float(len(intersection))/len(union) if union else 0

        full_graph = nx.read_weighted_edgelist(network_path)
        disease_genes = get_disease_genes(genelist_path)
   
        gse_args = [(full_graph, gene) for gene in disease_genes]

        with multiprocessing.Pool() as pool:
                clusters = pool.starmap(gse_helper, gse_args)

        # clusters = list(starmap(lambda g, s: greedy_source_expansion(g, source=s), gse_args))
        # for i, gene in enumerate(disease_genes):
        #         # can explore changing the cutoff parameter for this function later
        #         print(i)
        #         clusters.append(greedy_source_expansion(full_graph, source=gene))

        # get jaccard matrix of clusters

        n = len(clusters)
        jmat = [[0]*n for _ in range(n)]
        dick = {}
        for i in range(n - 1):
                jmat[i][i] = 1
                for j in range(i + 1, n):
                        jmat[i][j] = jmat[j][i] = jaccard(clusters[i], clusters[j])
                        num = round(jmat[i][j], 2)
                        dick[num] = dick.get(num, 0) + 1

        for i in range(n):
                jmat[n - 1][i] = jmat[i][n - 1]
        

        print(dick)
        
        rankings = get_cluster_rankings(clusters, full_graph, disease_genes)
        ranking_dict = {}

        for ranking in rankings:
                for gene, score in ranking:
                        ranking_dict[gene.name] = max(ranking_dict.get(gene.name, 0), score)
        
        return sorted([[Gene(x), ranking_dict[x]] for x in ranking_dict], key=lambda x: -x[1])


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
        # predictions = sorted(list(model.david_prioritize(1000, graph.graph)), key=lambda x: x[1], reverse=True)


        """
        unsupervised clustering - louvain, markov, walktrap
        """
        predictions = clustering(network_path, genelist_path, algorithm="walktrap")

        """ 
        supervised clustering
        """
        # predictions = supervised_clustering(network_path, genelist_path)


        with open(out_path, "w") as f:
                for gene, score in predictions:
                        f.write(f"{gene.name}\t{score}\n")

        end = time()
        print("Total time:", end-start)
        

if __name__ == "__main__":
        args = parser.parse_args()
        main(args.network, args.genelist, args.out)
        