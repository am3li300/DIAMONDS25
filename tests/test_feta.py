from t_map.feta.randomwalk import PreComputeRWR, RandomWalkWithRestart
from t_map.feta.glide import Glider 
from t_map.garbanzo.edgelist import EdgeListGarbanzo


def test_rwr():
    data = EdgeListGarbanzo(
        graph_path="./tests/data/unweighted_edge.list",
        gene_path="./tests/data/gene.list")

    model = RandomWalkWithRestart()
    model(data.get(0), data.graph)
    model([data.get(0), data.get(1)], data.graph)

    precompute = PreComputeRWR()
    precompute(data.get(0), data.graph)
    precompute([data.get(0), data.get(3)], data.graph)


def test_glider():
    TISSUE_PATH = "/cluster/tufts/cowenlab/mgelem01/E-MTAB-5214/tissue_files" + "/*.tsv"
    data = EdgeListGarbanzo(
        graph_path="./tests/data/unweighted_edge.list",
        gene_path="./tests/data/gene.list")
    import glob
    tissue_file = glob.glob(TISSUE_PATH)[0]
    print(tissue_file)
    model = Glider()
    model.setup(data.graph)
    model(data.get(0), data.graph, tissue_file=tissue_file)
