import networkx as nx
from typing import List
from t_map.garbanzo.garbanzo import Garbanzo
from t_map.gene.gene import Gene


class EdgeListGarbanzo(Garbanzo):
    def __init__(self, graph_path: str, gene_path: str):
        self._graph_path = graph_path
        self._gene_path = gene_path
        self._graph = self._read_graph_from_path(graph_path)
        self._genes = self._read_genes_from_path(gene_path)

    @property
    def graph(self) -> nx.Graph:
        return self._graph

    @property
    def graph_path(self) -> str:
        return self._graph_path

    @property
    def gene_path(self) -> str:
        return self._gene_path

    @property
    def genes(self) -> List[Gene]:
        return self._genes

    def __len__(self) -> int:
        return len(self._genes)

    def get(self, i: int) -> Gene:
        return self._genes[i]

    def _read_graph_from_path(self, path: str, weighted=True) -> nx.Graph:
        if weighted:
            return nx.read_weighted_edgelist(path)
        return nx.read_edgelist(path)

    def _read_genes_from_path(self, path: str) -> List[Gene]:
        genes = []
        with open(path) as f:
            for line in f.readlines():
                gene = line.strip()
                genes.append(Gene(name=gene))
        return genes
