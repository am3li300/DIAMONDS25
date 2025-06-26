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
parser.add_argument('-a', '--algorithms',
                    choices=['all', 'rwr', 'glider', 'dada', 'glider+dada'])
parser.add_argument('-g', '--gene-set', type=str,
                    required=True, help="path to genesets")
parser.add_argument('-th', '--thresholds', type=list,
                    default=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

models = {
    'rwr': RandomWalkWithRestart,
    'glider': Glider,
    'dada': Dada,
    'glider+dada': lambda: Glider(with_dada=True)
}


def main(args: argparse.Namespace):
    # Using glob get all of the tissue files.
    tissue_files = glob.glob(args.tissue_path + '/*.tsv')
    if args.algorithms == 'all':
        algs = ['rwr', 'glider', 'dada', 'glider+dada']
    else:
        algs = [args.algorithms]
    for tissue_file in tissue_files:
        huri = Huri(args.gene_set, with_hugo=True)
        stringdb = StringDB(args.gene_set)
        merged = merge([huri, stringdb])
        merged._graph = reweight_graph_by_tissue(merged.graph, tissue_file)
        for alg in algs:
            model = models[alg]()
            if isinstance(model, PreComputeFeta):
                model.setup(merged.graph)

            scoring = HummusScore(score_type=ScoreTypes.TOP_K, k=len(merged.graph.nodes()))
            runner = Hummus(merged, with_scoring=scoring)
            for i, test_runner in enumerate(runner.with_cv(k_fold="LOO")):
                with test_runner as (genes, graph, fn):
                    predictions = model(genes, graph)
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
            plt.plot(false_pos, true_pos, label=alg)
            plt.legend()
            plt.xlabel("False Positive")
            plt.ylabel("True Positive")
    
            print(scoring._roc)
   
    plt.savefig(tissue_file.strip(".tsv") + '.jpg')        
    plt.show()           


if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
