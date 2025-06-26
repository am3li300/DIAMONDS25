import networkx as nx
from abc import ABC, abstractmethod, abstractproperty
from t_map.gene.gene import Gene
from typing import Union, Tuple, Set, List, Optional
from t_map.feta.description import Description
import dill as pickle


class Feta(ABC):
    """
        Feta is the model which holds different methods of prioritization.
        For any method there are three functions:

            description: Returns a description of the model.

            prioritize: Prioritizes a set of genes in a graph based
                            on a prioritization algorithm
                Parameters: Takes in a list of Genes
                            Optional: a networkx Graph
                Returns a set of tuples with the gene and float prioritization


            call: Calls prioritization on either a singular
                    gene or a list of genes
                Parameters: a list of disease gene or a singular disease gene
                Returns a set of tuples with the gene and float prioritization
                Exceptions:throws and exception if not given a gene
                        or a list of genes

        Within the Feta model there is also a subclass
        for any method that requires precomputing which
        requires an additional setup function.

            setup: Does precomputing that prioritize needs for a given method
                Parameters: a networkx Graph

    """
    @abstractproperty
    def description(self) -> Description:
        ...

    @abstractmethod
    def prioritize(self, disease_gene: List[Gene],
                   graph: Union[nx.Graph, None], tissue_file: Optional[str] = None, **kwargs) -> Set[Tuple[Gene, float]]:
        ...

    def __call__(self, disease_gene: Union[Gene, List[Gene]],
                 graph: Union[nx.Graph, None], **kwargs) -> Set[Tuple[Gene, float]]:
        if isinstance(disease_gene, Gene):
            return self.prioritize([disease_gene], graph, **kwargs)
        elif isinstance(disease_gene, list) and \
                all(isinstance(x, Gene) for x in disease_gene):
            sanitized_disease_gene = list(filter(
                lambda x: x.name in graph.nodes, disease_gene))
            if len(sanitized_disease_gene) != len(disease_gene):
                print(
                    f"{len(disease_gene) - len(sanitized_disease_gene)} genes were not found in the graph")
            return self.prioritize(sanitized_disease_gene, graph, **kwargs)
        else:
            Exception(
                "Need to pass a Gene or a list of Genes, instead passed:",
                type(disease_gene))


class PreComputeFeta(Feta):
    """
        Within the Feta model there is also a subclass
        for any method that requires precomputing which
        requires an additional setup function.

            setup: Does precomputing that prioritize needs for a given method
                Parameters: a networkx Graph

    """
    @abstractmethod
    def setup(self, graph: nx.Graph, *args, **kwargs):
        ...

    def save(self, file_name: str):
        with open(file_name, 'wb') as f:
            pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)

    @classmethod
    def load(cls, file_name: str):
        with open(file_name, 'rb') as f:
            return pickle.load(f)
