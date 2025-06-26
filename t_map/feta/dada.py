import networkx as nx
from t_map.gene.gene import Gene
from typing import Union, Tuple, Set, List, Optional
import numpy as np
from t_map.feta.feta import Feta
from t_map.feta.description import Description


class Dada(Feta):

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

        pr_no_restart = nx.pagerank(
            graph, alpha=1,
            personalization=personalization)

        for gene_id, score in pr.items():
            gene = Gene(name=gene_id)
            if pr_no_restart[gene_id] == 0:
                score = 0
            else:
                score = np.log10(score / pr_no_restart[gene_id])
            gene_pair = (gene, score)
            data_set.add(gene_pair)
        return data_set
