from __future__ import annotations
import networkx as nx
import random
from functools import partial
from t_map.garbanzo.garbanzo import Garbanzo
from typing import List, Tuple, Union, Generator, Literal, Optional
from typing import ContextManager
from t_map.gene.gene import Gene
from t_map.hummus.hummus_score import HummusScore, ScoreFn
from math import floor
from sklearn.model_selection import KFold

KnownGenes = List[Gene]
HeldoutGenes = List[Gene]


class Hummus:

    def __init__(self, data: Garbanzo,
                 with_scoring: Optional[HummusScore] = None,
                 train_test_split_ratio: Optional[float] = None):
        self.data = data
        self._genes = [self.data.get(i) for i in range(len(data))]
        self._score_module = with_scoring
        self._with_cv = False
        self._train_test_split_ratio = train_test_split_ratio

    def with_cv(self, k_fold: Union[int, Literal["LOO"]] = 5) -> Hummus:
        if k_fold == "LOO":
            k_fold = len(self.data)
        self._with_cv = True
        self._k = k_fold
        self._kf = KFold(n_splits=self._k)
        return self

    def train(self) -> HummusTrainer:
        '''
            checks if a _k has been set, if not
            will initialize a default train split
        '''
        if self._with_cv is False:
            train_size = floor(len(self.data) * self._train_test_split_ratio)
            items = {i for i in range(len(self.data))}
            self._train = random.sample(items, k=train_size)
            self._test = list(items - set(self._train))

        return HummusTrainer(self.data, self._train, self._score_module)

    def test(self) -> HummusTester:
        # If the user hasn't selected cv mode and hasn't run the
        # training loop assume they want to test on whole dataset
        if self._with_cv is False and not hasattr(self, "_test"):
            self._test = [i for i in range(len(self.data))]

        return HummusTester(self.data, self._test, self._score_module)

    def __iter__(self) -> Generator[Union[Tuple[HummusTrainer, HummusTester],
                                          HummusTester],
                                    None,
                                    None]:
        assert self._with_cv, f"{self} has not been set to cv mode"

        if self._train_test_split_ratio:
            raise NotImplementedError

        else:
            for (train, test) in self._kf.split(self._genes):
                self._train = train
                self._test = test
                yield self.test()


class HummusTester(ContextManager):

    def __init__(self, data: Garbanzo, test_idx: List[int],
                 with_scoring: Union[None,
                                     HummusScore] = None):

        self._data: Garbanzo = data
        self._test: List[int] = test_idx
        self._i = 0
        self._heldout_genes = self._collect_heldout_genes()
        self._with_scoring = with_scoring
        if self._with_scoring is not None:
            self._collect_known_genes()

    def _collect_known_genes(self):
        '''
            Will set self._known_targets to all the
            known disease/target genes. (perhaps we need
            to change this in the future)
        '''
        genes: List[Gene] = [
            self._data.get(i) for i in range(len(self._data))]
        self._known_genes = list(set(genes) - set([self._data.get(i)
                                                   for i in self._test]))

    def _collect_heldout_genes(self) -> HeldoutGenes:
        return [self._data.get(i) for i in self._test]

    def __enter__(self) -> Union[Tuple[KnownGenes, nx.Graph, HeldoutGenes],
                                 Tuple[KnownGenes, nx.Graph, ScoreFn]]:
        """
            Iterator over heldout genes. If a scoring module was provided
            it will use the scoring modules score function to score the heldout
            gene, otherwise it is the responsibility of the implementor to
            appropriately score the heldout gene.

            Returns:
                Tuple[KnownGenes, nx.Graph, HeldoutGenes]: List of known genes,
                    the graph and the heldout genes.
                Tuple[KnownGenes, nx.Graph, ScoreFn]: List of known genes, the
                    graph and the score function.
        """

        if isinstance(self._with_scoring, HummusScore):
            self._with_scoring.test()
            self._with_scoring.reset_ground_truth()

            score_fn = partial(self._with_scoring.score,
                               self._heldout_genes)

            return (self._known_genes, self._data.graph, score_fn)

        else:
            return (self._known_genes,
                    self._data.graph,
                    self._heldout_genes)

    def __exit__(self, exc_type, exc_value, traceback):
        ...


class HummusTrainer(ContextManager):
    """
        Ideally I expect this to be a training class, but
        after much thought I am a bit puzzled as to how to
        best structure this.
    """

    def __init__(self, data: Garbanzo, train_idx: List[int],
                 with_scoring: Union[None,
                                     HummusScore] = None):

        self._data: Garbanzo = data
        self._train: List[int] = train_idx
        self._i = 0
        self._with_scoring = with_scoring

        if self._with_scoring is not None:
            self._collect_known_genes()

    def _collect_known_genes(self):
        self._known_targets: List[Gene] = [
            self._data.get(i) for i in self._train]
        self._known_targets_set = set(self._known_targets)

    def __enter__(self) -> Union[Tuple[KnownGenes, nx.Graph],
                                 Tuple[KnownGenes, nx.Graph]]:

        if isinstance(self._with_scoring, HummusScore):
            self._with_scoring.train()
            return (list(self._known_targets), self._data.graph)

        else:
            return (list(self._known_targets), self._data.graph)

    def __exit__(self, exc_type, exc_value, traceback):
        ...
