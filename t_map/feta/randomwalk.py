import networkx as nx
from t_map.gene.gene import Gene
from typing import Union, Tuple, Set, Any, List, Optional
import numpy as np
from networkx import to_dict_of_lists
from t_map.feta.feta import Feta, PreComputeFeta
from t_map.feta.description import Description

class RandomWalkWithRestart(Feta):

    def __init__(self, alpha=0.85):
        self._desc = Description(
            requires_training=False,
            training_opts=None,
            hyper_params={"alpha": alpha})

    def description(self) -> Description:
        return self.__desc

    def prioritize(self, disease_gene: List[Gene],
                   graph: Union[nx.Graph, None], **kwargs) -> Set[Tuple[Gene, float]]:
        personalization = dict()
        personalization = {gene.name: 1 for gene in disease_gene}
        pr = nx.pagerank(
            graph,
            alpha=self._desc.hyper_params["alpha"],
            personalization=personalization)
        data_set = set()

        for gene_id, score in pr.items():
            gene = Gene(name=gene_id)
            gene_pair = (gene, score)
            data_set.add(gene_pair)
        return data_set


class PreComputeRWR(PreComputeFeta):
    def __init__(self, alpha=0.85, **kwargs):
        self.__desc = Description(
            requires_training=False,
            training_opts=None,
            hyper_params={"alpha": alpha})

    def description(self) -> Description:
        return self.__desc

    def prioritize(self, disease_gene: List[Gene],
                   graph: Union[nx.Graph, None], **kwargs) -> Set[Tuple[Gene, float]]:     
        gene_list = to_dict_of_lists(graph)
        indicator_vector = []
        disease_gene_names = [gene.name for gene in disease_gene]
        for genes in gene_list:
            if (genes in disease_gene_names):
                indicator_vector.append(1)
            else:
                indicator_vector.append(0)
        indicator_vector = np.asmatrix(indicator_vector).transpose()
        final_vector = self.setup(graph) * indicator_vector
        gene_list = to_dict_of_lists(graph)
        final_rankings = set()
        for i, genes in enumerate(gene_list):
            pair = (genes, final_vector[i, 0])
            final_rankings.add(pair)

        return final_rankings

    def setup(self, graph: nx.Graph, *args, **kwargs) -> Any:
        """
        adjacency marix * all-one vector(matrix) = sum of each row
        flatten(sum of each row vector)
        convert vector into a diagonal matrix np.diag
        gives you D
        compute inv(D) * adjacency matrix = markov
        """

        temp = nx.adjacency_matrix(graph)
        n, _ = temp.shape
        d = (temp @ np.ones((n, 1))).flatten()
        d_i = np.array([1/dd if dd != 0 else 0 for dd in d])
        markov = np.diag(d_i) @ temp
        identity = np.identity(len(markov))

        prob = self.__desc.hyper_params["alpha"]

        # s = c(I - (1-c)markov)^-1r
        inverse_matrix = np.linalg.inv(
            identity - (1 - prob) * np.asmatrix(markov))
        s = prob * (inverse_matrix)

        return s
