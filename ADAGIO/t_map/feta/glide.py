from __future__ import annotations
import networkx as nx
import numpy as np
from contextlib import contextmanager
from copy import deepcopy

from t_map.gene.gene import Gene
from typing import Any, List, Set, Tuple, Union, Optional
from t_map.feta.description import Description
from t_map.feta.feta import PreComputeFeta
from t_map.feta.dada import Dada
from t_map.feta.randomwalk import RandomWalkWithRestart
from t_map.garbanzo.transforms.tissue_reweight import reweight_graph_by_tissue
from gfunc.command import glide_mat
from networkx.algorithms import tree


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

    def set_add_edges_amount(self, amount: int) -> None:
        self.to_add = amount

    def set_remove_edges_percent(self, percentage: float) -> None:
        self.to_remove = percentage

    def reweight_graph(self) -> None:
        self.graph = self._update_edge_weights(
            self.graph, self.gmat, self.gmap)

    def reset_graph(self) -> None:
        self.graph = deepcopy(self._original_graph)

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

    def add_edges_around_node(self, node: str, new_edges_count: int, variant: str = "none") -> List[Tuple[int, int]]:
        indexes = self._get_sorted_similarity_indexes()
        node_idx = self.gmap[node]
        graph_edges = self.graph.edges()
        add_cnt = 0
        pairs_to_add = []
        max_edges = self.graph.edges(node)
        for i, j in indexes:
            if node_idx == i or node_idx == j:
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
        if hasattr(self, "to_add"):
            for disease_gene in disease_genes:
                to_add_pairs = self.add_edges_around_node(
                    disease_gene.name, self.to_add, variant)
                for (i, j) in to_add_pairs:
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

        if hasattr(self, '__dada'):
            return self.__dada.prioritize(disease_genes, graph)
        else:
            rwr = RandomWalkWithRestart(alpha=0.85)
            return rwr.prioritize(disease_genes, graph)

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
