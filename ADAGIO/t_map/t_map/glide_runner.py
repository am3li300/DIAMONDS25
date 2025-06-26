import networkx as nx
from t_map.garbanzo.garbanzo import Garbanzo
from t_map.gene.gene import Gene
from t_map.feta.glide import Glider


class Genedise(Garbanzo):
    def __init__(self, graph_path, disease_path):
        self._graph_path = graph_path
        self._graph = nx.read_graphml(self._graph_path)
        self.disease_path = disease_path
        self._genes = []

        with open(disease_path) as f:
            for line in f:
                self._genes.append(Gene(line.strip()))

    @property
    def graph(self) -> nx.Graph:
        return self._graph

    @property
    def graph_path(self) -> str:
        return self._graph_path

    @property
    def gene_path(self) -> str:
        return self.disease_path

    @property
    def genes(self) -> list:
        return self._genes

    def get(self, i: int) -> Gene:
        return self._genes[i]

    def __len__(self):
        return len(self._genes)


def run_glide(graph_path, disease_path):
    model = Glider()
    data = Genedise(graph_path, disease_path)
    model.setup(data.graph)
    return model.prioritize(data.genes)
