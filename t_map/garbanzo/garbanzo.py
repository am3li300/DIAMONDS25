import networkx as nx
from abc import ABC, abstractmethod, abstractproperty
from t_map.gene.gene import Gene
from typing import List


class Garbanzo(ABC):
    @abstractproperty
    def graph(self) -> nx.Graph:
        ...

    @abstractproperty
    def graph_path(self) -> str:
        ...

    @abstractproperty
    def gene_path(self) -> str:
        ...

    @abstractproperty
    def genes(self) -> List[Gene]:
        ...

    @abstractmethod
    def get(self, i: int) -> Gene:
        ...

    @abstractmethod
    def __len__(self) -> int:
        ...
