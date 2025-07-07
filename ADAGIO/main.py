import argparse
from t_map.garbanzo.edgelist import EdgeListGarbanzo
from t_map.feta.glide import ADAGIO
from time import time
from heapq import heapify, heappush, heappop
import networkx as nx
from t_map.gene.gene import Gene

parser = argparse.ArgumentParser()
parser.add_argument('--network', '-n', type=str, required=True, help="Path to edgelist, assumes tab separated.")
parser.add_argument('--genelist', '-g', type=str, required=True, help="Path to genelist, assumes genes are newline separated and in the network.")
parser.add_argument('--out', '-o', type=str, default="adagio.out",help="Path to output results")


def clustering(network_path, genelist_path, out_path):
        def create_cluster_graph(graph, cluster):
                new_graph = graph.copy()
                for node in graph.nodes:
                        if node not in cluster:
                                new_graph.remove_node(node)
                
                return new_graph

        # make graph from network_path and get clusters
        full_graph = nx.read_weighted_edgelist(network_path)
        clusters = nx.community.louvain_communities(full_graph)

        n = len(clusters)
        gene_mapping = {} # store index for gene_groups
        gene_groups = [[] for _ in range(n)] # store the actual groups

        # store seed nodes as gene objects in a list from genelist_path
        with open(genelist_path) as f:
            for line in f.readlines():
                gene = line.strip()
                for i in range(n):
                        if gene in clusters[i]:
                                gene_mapping[gene] = i
                                gene_groups[i].append(Gene(name=gene))
                                
        # generate rankings for each cluster
        rankings = []
        for i, cluster in enumerate(clusters):
                if (len(gene_groups[i]) == 0):
                        continue
                sub_graph = create_cluster_graph(full_graph, cluster)
                model = ADAGIO()
                model.setup(sub_graph)
                model.set_add_edges_amount(20)
                rankings.append(sorted(list(model.prioritize(gene_groups[i], sub_graph)), key=lambda x: -x[1]))
        
        # merge rankings (merge k sorted lists)
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
        
        
def main(network_path: str, genelist_path: str, out_path: str="adagio.out"):
        start = time()

        """
        ADAGIO with no cluster handling
        """
        # graph = EdgeListGarbanzo(network_path, genelist_path)
        # print(len(graph.graph.edges))
        # model = ADAGIO()
        # model.setup(graph.graph)
        # model.set_add_edges_amount(20) # this will add edges to the graph

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
        constant k for disease nodes only (k = 20)
        """
        predictions = clustering(network_path, genelist_path, out_path) # sorted(list(model.prioritize(graph.genes, graph.graph)), key=lambda x: x[1], reverse=True)

        with open(out_path, "w") as f:
                for gene, score in predictions:
                        f.write(f"{gene.name}\t{score}\n")

        end = time()
        print("Total time:", end-start)
        

if __name__ == "__main__":
        args = parser.parse_args()
        main(args.network, args.genelist, args.out)
        