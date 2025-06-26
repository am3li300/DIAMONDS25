

from collections import defaultdict
from typing import Tuple, List
from t_map.gene.gene import Gene


class Witness:

    def __init__(self, l_1: List[Tuple[Gene, float]],
                 *args: List[Tuple[Gene, float]]):
        self.l_1 = l_1
        self.l_2s = args

        self.l_1_map = dict()
        self.l_2_map = defaultdict(list)

        for i, (g, _) in enumerate(l_1):
            self.l_1_map[g.name] = i
        for l in self.l_2s:
            for i, (g, _) in enumerate(l):
                self.l_2_map[g.name].append(i)

    def find_range_of_witness(self, gene: Gene,
                              buffer=lambda x: x/2) -> Tuple[int, int]:
        start = list(filter(lambda x: x[1][0].name ==
                            gene.name, enumerate(self.l_1)))[0][0]
        end = start
        while all(map(lambda x: buffer(x) >= end, self.l_2_map[gene.name])):
            end += 1
        return start, end

    def find_all_witnesses(self, buffer=lambda x: x/2) -> List[Tuple[Gene, int, int]]:
        wittness_intervals = []
        for gene, _ in self.l_1:
            start, end = self.find_range_of_witness(gene, buffer)
            wittness_intervals.append((gene, start, end))
        return wittness_intervals

    @classmethod
    def compute_witness(cls, l_1: List[Tuple[Gene, int, int]]) -> List[int]:

        # Ensure data is formatted correctly.
        l_1 = list(filter(lambda x: x[1] != x[2], l_1))
        max_range = max(l_1, key=lambda x: x[2])[2]
        l = [0] * max_range
        for (_, (_, start, end)) in enumerate(l_1):
            for i in range(start, end):
                l[i] += 1

        return l
