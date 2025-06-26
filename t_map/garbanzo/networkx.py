import networkx as nx
from typing import List, Optional
from t_map.gene.gene import Gene
from t_map.garbanzo.garbanzo import Garbanzo


class Networkx(Garbanzo):
    def __init__(self, graph: nx.Graph, genes: Optional[List[Gene]] = None):
        self._graph = graph
        self._genes = genes if genes is not None else []

    @property
    def graph(self) -> nx.Graph:
        return self._graph

    @property
    def genes(self) -> List[Gene]:
        return self._genes

    @property
    def gene_path(self) -> str:
        raise Exception("Genes were passed as a parameter")

    @property
    def graph_path(self) -> str:
        raise Exception("Graph were passed as a parameter")

    def get(self, i: int) -> Gene:
        return self._genes[i]

    def __len__(self) -> int:
        return len(self._genes)
