#!/cluster/tufts/cowenlab/.envs/tmap/bin/python
import glob
import argparse
import matplotlib.pyplot as plt
from t_map.garbanzo.huri import Huri
from t_map.hummus.hummus import Hummus
from t_map.hummus.hummus_score import HummusScore, ScoreTypes
from t_map.feta.randomwalk import RandomWalkWithRestart
from t_map.garbanzo.stringdb import StringDB
from t_map.garbanzo.transforms.tissue_reweight import reweight_graph_by_tissue
from t_map.garbanzo.merge import merge
from t_map.feta.glide import Glider
from t_map.feta.dada import Dada
from t_map.feta.feta import PreComputeFeta

parser = argparse.ArgumentParser(description='Tissue reweight by algs')
parser.add_argument('-t', '--tissue-path', type=str, required=True,
                    help='Path to folder that contains tissue presence data')
parser.add_argument('--lamb', type=float, default=1)
parser.add_argument('--beta', type=float, default=1000)
parser.add_argument('-g', '--gene-set', type=str,
                    required=True, help="path to genesets")
parser.add_argument('-th', '--thresholds', type=list,
                    default=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

def main(args: argparse.Namespace):
    # Using glob get all of the tissue files.
    tissue_files = glob.glob(args.tissue_path + '/*.tsv')
    model = Glider(lamb=args.lamb, glide_beta = args.beta)
    for tissue_file in tissue_files:
        huri = Huri(args.gene_set, with_hugo=True)
        stringdb = StringDB(args.gene_set)
        merged = merge([huri, stringdb])
        if isinstance(model, PreComputeFeta):
	        model.setup(merged.graph)

        scoring = HummusScore(score_type=ScoreTypes.TOP_K, k=len(merged.graph.nodes()))
        runner = Hummus(merged, with_scoring=scoring)
        for i, test_runner in enumerate(runner.with_cv(k_fold="LOO")):
            with test_runner as (genes, graph, fn):
                predictions = model(genes, graph, tissue_file=tissue_file)
                fn(predictions)

        # This is how we get the gene set predictions without LOO cross
        # validation.
        gene_list = model(merged.genes, merged.graph)
        # TODO: Figure out how to score this gdamn thing.
        scoring.compute_roc()
        for threshold in args.thresholds:
            scoring.auc_roc(threshold)

        true_pos = [x[0] for x in scoring._roc]
        false_pos = [x[1] for x in scoring._roc]
        result_file = tissue_file.strip(".tsv")
        result_file = f"{result_file}_{args.beta}_{args.lamb}_glider.result"
        save_results(true_pos, false_pos, result_file)



def save_results(true_pos: list[float], false_pos: list[float], file_name: str) -> None:
    with open(file_name, "w") as f:
        true_pos_str = ",".join(map(str, true_pos))
        false_pos_str = ",".join(map(str, false_pos))
        f.write(f"TP:\t({true_pos_str})\n")
        f.write(f"FP:\t({false_pos_str})\n")

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
