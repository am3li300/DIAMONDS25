#!/usr/bin/env python
from t_map.garbanzo.edgelist import EdgeListGarbanzo
from t_map.garbanzo.huri import Huri
from t_map.garbanzo.merge import merge
from t_map.garbanzo.stringdb import StringDB


def test_edgelist():
    # testing for side effect
    EdgeListGarbanzo(
        graph_path="./tests/data/unweighted_edge.list",
        gene_path="./tests/data/gene.list")


def test_huri():
    data = Huri("./tests/data/gene.list")
    # maybe not the best test as only works for posix
    assert(data.graph_path == "/tmp/huri/huri.tsv")
    assert data.get(0)
    assert data.graph


def test_merge():
    huri = Huri("./tests/data/gene.list")
    stringdb = StringDB("./tests/data/gene.list")
    merged = merge([huri, stringdb])

    assert merged.graph
    assert list(huri.graph.nodes)[0] in merged.graph
