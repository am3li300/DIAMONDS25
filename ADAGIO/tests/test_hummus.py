from t_map.garbanzo.edgelist import EdgeListGarbanzo
from t_map.hummus.hummus import Hummus
from t_map.hummus.hummus_score import HummusScore


def test_cv():
    data = EdgeListGarbanzo(
        graph_path="./tests/data/unweighted_edge.list",
        gene_path="./tests/data/gene.list")
    score_module = HummusScore()
    runner = Hummus(data, with_scoring=score_module)

    # Check if we assert gracefully in cv is not set
    try:
        for tester in runner:
            ...
    except AssertionError:
        ...

    for tester in runner.with_cv("LOO"):
        with tester as (known_genes, graph, score_fn):
            ...
