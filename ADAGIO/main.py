import argparse
from t_map.garbanzo.edgelist import EdgeListGarbanzo
from t_map.feta.glide import ADAGIO
from time import time

parser = argparse.ArgumentParser()
parser.add_argument('--network', '-n', type=str, required=True, help="Path to edgelist, assumes tab separated.")
parser.add_argument('--genelist', '-g', type=str, required=True, help="Path to genelist, assumes genes are newline separated and in the network.")
parser.add_argument('--out', '-o', type=str, default="adagio.out",help="Path to output results")

def main(network_path: str, genelist_path: str, out_path: str="adagio.out"):
        start = time()
        graph = EdgeListGarbanzo(network_path, genelist_path)
        print(len(graph.graph.edges))
        model = ADAGIO()
        model.setup(graph.graph)
        model.set_add_edges_amount(20)
        predictions = sorted(list(model.prioritize(graph.genes, graph.graph)), key=lambda x: x[1], reverse=True)
        with open(out_path, "w") as f:
                for gene, score in predictions:
                        f.write(f"{gene.name}\t{score}\n")

        end = time()
        print("Total time:", end-start)
        

if __name__ == "__main__":
        args = parser.parse_args()
        main(args.network, args.genelist, args.out)
        