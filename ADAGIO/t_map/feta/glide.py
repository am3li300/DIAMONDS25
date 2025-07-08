from __future__ import annotations
import networkx as nx
import numpy as np
from contextlib import contextmanager
from copy import deepcopy
from math import floor
from collections import defaultdict
import pandas as pd

from t_map.gene.gene import Gene
from typing import Any, List, Set, Tuple, Union, Optional
from t_map.feta.description import Description
from t_map.feta.feta import PreComputeFeta
from t_map.feta.dada import Dada
from t_map.feta.randomwalk import RandomWalkWithRestart
from t_map.garbanzo.transforms.tissue_reweight import reweight_graph_by_tissue
from gfunc.command import glide_mat
from networkx.algorithms import tree

# import sys
# import os
# script_dir = os.path.abspath(os.path.join(__file__, '../../DScript'))
# sys.path.append(script_dir)

from DScript.predict import dscript_predict

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class ADAGIO(PreComputeFeta):

    def __init__(self, is_annotated=True, lamb: int = 1,
                 is_normalized: bool = False, glide_alph: float = 0.1,
                 glide_delta: int = 1,
                 with_dada: bool = False, dada_alpha: float = 0.85,
                 glide_loc="cw_normalized", extend_network: bool = False,
                 per_node_new_edges_count: int = 0,
                 global_new_edges_percentage: float = 0, **kwargs) -> None:
        self.__desc = Description(requires_training=False, training_opts=False,
                                  hyper_params={
                                      "lamb": lamb,
                                      "is_normalized": is_normalized,
                                      "glide_alph": glide_alph,
                                      "glide_delta": glide_delta,
                                      "glide_loc": glide_loc,
                                      "per_node_new_edges_count": per_node_new_edges_count,
                                      "global_new_edges_percentage": global_new_edges_percentage,
                                  })
        if with_dada:
            self.__dada = Dada(alpha=dada_alpha)

    def description(self) -> Description:
        return self.__desc

    def setup(self, graph: nx.Graph, *args, **kwargs) -> Any:
        elist = list(graph.edges(data=True))
        try:
            elist = list(
                map(lambda x_y_z: (x_y_z[0], x_y_z[1], x_y_z[2]['weight']), elist))
        except KeyError as _:
            print(bcolors.WARNING +
                  "Warning: Could not detect edgeweights, defaulting to 1." + bcolors.ENDC)
            elist = list(
                map(lambda x_y_z: (x_y_z[0], x_y_z[1], 1), elist))

        self.graph = deepcopy(graph)
        self.gmat, self.gmap = glide_mat(
            elist, self.description().hyper_params)
        self.rgmap = {self.gmap[key]: key for key in self.gmap}
        self._original_graph: nx.Graph = deepcopy(graph)
        # self._get_sorted_similarity_indexes()
        # self._get_sorted_similarity_indexes(descending=True)
        self.reweight_graph()
        self.candidate_pairs = pd.DataFrame()
        self.seq_dict = self.setup_fasta_dict()
        print("Finished graph set up...")

    def setup_fasta_dict(self) -> dict[str, str]:
        # data/UP000005640_9606.fasta
        file = open("../data/UP000005640_9606.fasta") #input("Enter file path for gene sequences: "))
        mapping = defaultdict(str)
        gene = ""
        sequence = []
        for line in file:
            if line[0] == '>':
                if gene:
                    mapping[gene] = "".join(sequence)
                    sequence = []

                x = "".join(line.split('|')[2])
                humanIndx = x.index("_HUMAN")
                gene = x[:humanIndx]

            else:
                sequence.append(line.strip())

        mapping[gene] = "".join(sequence)
        return mapping

    """
    change to handle dictionary argument
    """
    def set_add_edges_amount(self, amount: int | dict[str, int]) -> None:
        self.k_mat = amount if isinstance(amount, dict) else defaultdict(lambda: amount)

    def set_remove_edges_percent(self, percentage: float) -> None:
        self.to_remove = percentage

    def reweight_graph(self) -> None:
        self.graph = self._update_edge_weights(
            self.graph, self.gmat, self.gmap)

    def reset_graph(self) -> None:
        self.graph = deepcopy(self._original_graph)

    """
    adaptive k matrix
    """
    def construct_k_mat(self, graph, genes: List[Gene]) -> dict[str, int]:
        nodes = set(graph.nodes)
        total_sum = sum(graph.degree[node] for node in nodes)
        avg_degree = total_sum//len(nodes)
        max_edges_to_add = defaultdict(int)

        print("Number of queried genes:", len(genes))
        subset = [node.name for node in genes if node.name in nodes]
        print("Number of queried genes in graph:", len(subset))

        clustering_coefficients = nx.clustering(graph, nodes=subset)
        for name in subset:
            # adaptive k function
            max_edges_to_add[name] = floor(avg_degree*(1-clustering_coefficients[name]))

        return max_edges_to_add

    """
    slow
    """
    def get_k_value_for_node(self, graph, node):
        nodes = list(graph.nodes)
        total_sum = sum(graph.degree[gene] for gene in nodes)
        avg_degree = total_sum//len(nodes)
        return floor(avg_degree*(1-nx.clustering(graph, node.name)))

    def add_new_edges(self, global_new_edges_percentage: float,
                      in_place: bool = False) -> Optional[nx.Graph]:
        graph = deepcopy(self.graph)
        if global_new_edges_percentage > 0:
            new_edges_count = int(
                global_new_edges_percentage * len(graph.edges()))
            indexes = self._get_sorted_similarity_indexes()
            edges = graph.edges()
            print(bcolors.OKGREEN +
                  "Filtering {} new edges.".format(len(edges)) + bcolors.ENDC)
            avail_indexes = [(i, j) for i, j in indexes if (
                self.rgmap[i], self.rgmap[j]) not in edges][:new_edges_count]

            print(bcolors.OKGREEN +
                  "Adding {} new edges.".format(len(avail_indexes)) + bcolors.ENDC)
            for i, j in avail_indexes:
                graph.add_edge(
                    self.rgmap[i], self.rgmap[j], weight=self.gmat[i][j])

        if in_place:
            self.graph = graph
        else:
            return graph

    def remove_old_edges(self, edge_percentage_removal: float,
                         in_place: bool = False) -> Optional[nx.Graph]:
        graph = deepcopy(self.graph)
        if edge_percentage_removal > 0:
            removal_edges_count = int(
                edge_percentage_removal * len(graph.edges()))
            edges = graph.edges()
            indexes = self._get_sorted_similarity_indexes(descending=True)
            avail_indexes = []
            i = 0
            print(len(indexes))
            for i, j in indexes:
                if (self.rgmap[i], self.rgmap[j]) in edges:
                    avail_indexes.append((i, j))
                    i += 1
                    if i == removal_edges_count:
                        break

            print(bcolors.OKGREEN +
                  "Removing {} edges.".format(len(edges)) + bcolors.ENDC)

            for u, v in avail_indexes:
                try:
                    graph.remove_edge(self.rgmap[u], self.rgmap[v])
                except:
                    pass
        if in_place:
            self.graph = graph
        else:
            return graph

    """
    send dscript scores as a parameter
    """
    def add_edges_around_node(self, node: str, new_edges_count: int, variant: str = "none") -> List[Tuple[int, int]]:
        if node not in self.gmap:
            print("{0} does not exist in self.gmap".format(node))
            return []

        node_idx = self.gmap[node]
        indexes = self._get_sorted_similarity_indexes()
        """
        call dscript
        """
        # self.candidate_pairs = pd.DataFrame(indexes, columns=["Gene1", "Gene2"])
        graph_edges = self.graph.edges()
        add_cnt = 0
        pairs_to_add = []
        max_edges = self.graph.edges(node)
        for i, j in indexes:
            if node_idx == i or node_idx == j:
                """
                add additional condition to check dscript score
                """
                if (self.rgmap[i], self.rgmap[j]) not in graph_edges:
                    pairs_to_add.append((i, j))
                    add_cnt += 1
                else:
                    if variant == "min":
                        break
                    elif variant == "max":
                        if add_cnt >= len(max_edges):
                            break
                    else:
                        pass
            if add_cnt == new_edges_count:
                break
        return pairs_to_add

    def remove_edges_around_node(self, node: str, edge_percentage_removal: float, mst: nx.Graph) -> List[Tuple[int, int]]:
        graph_edges = self.graph.edges(node)
        pairs_to_remove = []

        edges_idx = [(self.gmap[u], self.gmap[v], self.gmat[self.gmap[u]]
                      [self.gmap[v]]) for u, v in graph_edges]
        sorted_edges_idx = sorted(edges_idx, key=lambda x: x[2])
        num_to_remove = int(edge_percentage_removal * len(sorted_edges_idx))

        cnt = 0
        for i, j, w in sorted_edges_idx:
            if cnt >= num_to_remove:
                break
            edge = (self.rgmap[i], self.rgmap[j])
            if edge in mst:
                pass
            else:
                pairs_to_remove.append((i, j))
                cnt += 1

        return pairs_to_remove

    def _get_sorted_similarity_indexes(self, descending=False) -> List[Tuple[int, int]]:
        if hasattr(self, "_sorted_similarity_indexes") and not descending:
            return self._sorted_similarity_indexes
        elif hasattr(self, "_sorted_similarity_indexes_desc") and descending:
            return self._sorted_similarity_indexes_desc
        else:
            shape = self.gmat.shape
            indexes = np.unravel_index(np.argsort(
                self.gmat.ravel(), axis=None), shape)
            if descending:
                indexes = [(indexes[0][i], indexes[1][i])
                           for i in range(len(indexes[0]))]
                self._sorted_similarity_indexes_desc = indexes
            else:
                indexes = list(reversed([(indexes[0][i], indexes[1][i])
                                        for i in range(len(indexes[0]))]))
                self._sorted_similarity_indexes = indexes
            return indexes

    def _update_edge_weights(self, graph: nx.Graph,
                             gmat: np.ndarray,
                             gmap: dict) -> nx.Graph:
        for u, v in graph.edges():
            graph[u][v]['weight'] = gmat[gmap[u]][gmap[v]]
        return graph

    def prioritize(self, disease_genes: List[Gene],
                   graph: Union[nx.Graph, None],
                   tissue_file: Optional[str] = None,
                   variant: str = "none",
                   **kwargs) -> Set[Tuple[Gene, float]]:

        if tissue_file:
            graph = reweight_graph_by_tissue(graph, tissue_file)

        graph = deepcopy(self.graph)
        """
        testing dscript - candidate_pairs is currently invalid
        """
        # dscript_predict(self.candidate_pairs, "DScript/model.safetensors", "a.out", self.seq_dict, 0.5)

        print("The graph originally has " + str(len(graph.edges)) + " edges")

        if hasattr(self, "k_mat"): # originally to_add
            k = self.k_mat.default_factory()
            for disease_gene in disease_genes:
                to_add_pairs = self.add_edges_around_node(
                    disease_gene.name, k, variant)
                for (i, j) in to_add_pairs:
                    # print(self.rgmap[i], self.rgmap[j], i, j)
                    graph.add_edge(
                        self.rgmap[i], self.rgmap[j], weight=self.gmat[i][j])
        if hasattr(self, "to_remove"):
            mst = tree.maximum_spanning_edges(
                graph, algorithm="prim", data=False)
            for disease_gene in disease_genes:
                to_remove_pairs = self.remove_edges_around_node(
                    disease_gene.name, self.to_remove, mst)
                for (i, j) in to_remove_pairs:
                    graph.add_edge(
                        self.rgmap[i], self.rgmap[j], weight=self.gmat[i][j])

        print("The graph now has " + str(len(graph.edges)) + " edges")
        if hasattr(self, '__dada'):
            return self.__dada.prioritize(disease_genes, graph)
        else:
            rwr = RandomWalkWithRestart(alpha=0.85)
            return rwr.prioritize(disease_genes, graph)


    def david_prioritize_2(self, disease_genes: List[Gene],
                   graph: Union[nx.Graph, None],
                   tissue_file: Optional[str] = None,
                   variant: str = "none",
                   **kwargs) -> Set[Tuple[Gene, float]]:

        print("Inside prioritize function")
        if tissue_file:
            graph = reweight_graph_by_tissue(graph, tissue_file)
 
        graph = deepcopy(self.graph)
        if hasattr(self, "k_mat"):
            k_mat = self.construct_k_mat(graph, disease_genes)
            for disease_gene in disease_genes:
                k_i = k_mat[disease_gene.name]
                if k_i > 0:
                    pairs = self.add_edges_around_node(disease_gene.name,
                                                    k_i,
                                                    variant)
                    for (i, j) in pairs:
                        graph.add_edge(self.rgmap[i],
                                    self.rgmap[j],
                                    weight=self.gmat[i][j])

        print(len(graph.edges))
        if hasattr(self, '__dada'):
            return self.__dada.prioritize(disease_genes, graph)
        else:
            rwr = RandomWalkWithRestart(alpha=0.85)
            return rwr.prioritize(disease_genes, graph)

    """
    David function
    """
    def david_prioritize(self, max_threshold: int, graph: Union[nx.Graph, None], adaptive_k: bool = False, tissue_file: Optional[str] = None) -> Set[Tuple[Gene, float]]:

        # add most confident edges with adaptive k and maximum edge addition threshold
        if tissue_file:
            graph = reweight_graph_by_tissue(graph, tissue_file)
 
        graph = deepcopy(self.graph)
        genes = [Gene(name=node) for node in graph.nodes]
        if max_threshold > 0:
            k_mat = self.construct_k_mat(graph, genes) if adaptive_k else defaultdict(lambda: float("inf"))
            indexes = self._get_sorted_similarity_indexes()
            for i, j in indexes:
                node_indx1 = self.rgmap[i] # string names
                node_indx2 = self.rgmap[j]

                # treat graph as undirected
                val = min(k_mat[node_indx1], k_mat[node_indx2])
                if val > 0:
                    graph.add_edge(self.rgmap[i],
                        self.rgmap[j],
                        weight=self.gmat[i][j])

                    k_mat[node_indx1] -= 1
                    k_mat[node_indx2] -= 1
                    max_threshold -= 1
                    if not max_threshold:
                        break
                    
        print(len(graph.edges))
        if hasattr(self, '__dada'):
            return self.__dada.prioritize(genes, graph)

        else:
            rwr = RandomWalkWithRestart(alpha=0.85)
            return rwr.prioritize(genes, graph)

    @classmethod
    def ADAGIO_with_pickle(cls, file_name: str, reset: bool = True) -> ADAGIO:
        glide: ADAGIO = cls.load(file_name)
        return glide._context(reset=reset)

    @classmethod
    def with_reset(cls, model) -> ADAGIO:
        new_model = deepcopy(model)
        return new_model._context(reset=True)

    @contextmanager
    def _context(self, reset: bool = True) -> ADAGIO:
        if reset:
            self.reset_graph()
        yield self
