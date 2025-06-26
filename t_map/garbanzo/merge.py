import networkx as nx
from itertools import chain
from t_map.garbanzo.garbanzo import Garbanzo
from typing import Iterable, Union, List
from t_map.gene.gene import Gene


class MergedGarbanzo(Garbanzo):

    def __init__(self, graph: nx.Graph, genes: List[Gene]):
        self._graph = graph
        self._genes = genes

    @property
    def graph(self) -> nx.Graph:
        return self._graph

    @property
    def graph_path(self) -> str:
        raise Exception("Merged Graphs Do Not Have a graph path")

    @property
    def gene_path(self) -> str:
        raise Exception("Merged Graphs Do Not Have a gene path")

    @property
    def genes(self) -> List[Gene]:
        return self._genes

    def get(self, i: int) -> Gene:
        return self._genes[i]

    def __len__(self) -> int:
        return len(self._genes)


def add_weights(data: Garbanzo) -> Garbanzo:
    """
    Adds weight 1 to all edges in a graph.
    :param graph:
    :return: A graph with weights 1.
    """
    for u, v in data.graph.edges():
        data.graph[u][v]['weight'] = 1
    return data


def normalize(data: Garbanzo) -> Garbanzo:
    """
    Given a graph, normalizes the weights of the edges 
    by the maximum edge weight.
    """
    max_weight = max(data.graph.edges(data='weight'), key=lambda x: x[2])[2]
    print(max_weight)
    for u, v in data.graph.edges():
        data.graph[u][v]['weight'] = data.graph[u][v]['weight'] / max_weight
    return data


def merge(items: Iterable[Union[Garbanzo, nx.Graph]]) -> Garbanzo:

    if not any([isinstance(item, Garbanzo) for item in items]):
        raise Exception("At least one item needs to be of instance Garbanzo")

    if any([nx.is_weighted(item.graph) for item in items]):
        weighted_networks = [
            item for item in items if nx.is_weighted(item.graph)]
        unweighted_networks = [
            item for item in items if not nx.is_weighted(item.graph)]
        reweighted_networks = [add_weights(item)
                               for item in unweighted_networks]
        items = weighted_networks + reweighted_networks
        items = [normalize(item) for item in items]

    garbanzos = [item for item in items if isinstance(item, Garbanzo)]
    garbanzo_graphs = [g.graph for g in garbanzos]
    list_of_genes = [g.genes for g in garbanzos]
    garbanzo_genes = list(set(chain(*list_of_genes)))
    graphs = [item for item in items if isinstance(item, nx.Graph)]
    merged_graph = nx.compose_all(list(chain(garbanzo_graphs, graphs)))
    return MergedGarbanzo(merged_graph, garbanzo_genes)
