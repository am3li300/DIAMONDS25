"""
need to get network OR model pickle
need name of disease
need partitions subfolder name, default can be STRING
need method subfolder name
need to get seed genes and make partitions OR path to partitions folder
call and run adagio the specified number of times, write output to rankings/method subfolder/disease
label the rankings based on partitions, write to labels/method subfolder/disease
call positives on the labels, will have to manually save the graphs
maybe write positives output to a designated output file with the disease and method name, use write method that lets you add to file instead of overwriting
all done :)

"""

"""
python3 cross_validate.py \
  --network '../data/networks/STRING_protein_links_parsed.tsv' \
  --model 'adagio_model' \
  --disease 'hypertension' \
  --partition 'STRING' \
  --source 'genetic' \
  --method 0 \
  --jobs 2 \
  --folds 3
"""

from t_map.garbanzo.edgelist import EdgeListGarbanzo
from t_map.feta.glide import ADAGIO
from t_map.gene.gene import Gene

import argparse
from time import time
# from heapq import heapify, heappush, heappop

import networkx as nx
# from scipy.sparse import csr_matrix
# import markov_clustering as mc

# import igraph as ig

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# import pcst_fast
import numpy as np
from joblib import Parallel, delayed

from collections import defaultdict
import dill as pickle
import re


_MODEL = None

def _init_model(model_path):
    global _MODEL
    if _MODEL is None:
        with open(model_path, "rb") as f:
            _MODEL = pickle.load(f)

def _rank_from_paths(method_id, network_path, genelist_path, i):
    # Build the graph in the worker, not the parent
    graph = EdgeListGarbanzo(network_path, genelist_path)
    if method_id == 0:
        return i, sorted(list(_MODEL.prioritize(graph.genes, graph.graph)), key=lambda x: -x[1])

    elif method_id == 1:
        return i, sorted(list(_MODEL.david_prioritize_2(graph.genes, graph.graph)), key=lambda x: -x[1])

    elif method_id == 2:
        return i, sorted(list(_MODEL.david_prioritize(1000, graph.graph, True)), key=lambda x: -x[1])

    elif method_id == 3:
        pass

    else:
        # method_id == 4
        pass


parser = argparse.ArgumentParser()
parser.add_argument('--network', '-n', type=str, required=True, help="Path to edgelist, assumes tab separated.")
parser.add_argument('--model', '-m', type=str, required=True, help="Path to pickled model file")
parser.add_argument('--disease', '-d', type=str, required=True, help="Disease to cross validate")
parser.add_argument('--partition', '-p', type=str, required=True, help="Partition folder name (e.g. STRING)")
parser.add_argument('--source', '-s', type=str, required=True, help="Drug or genetic data")
parser.add_argument('--method', '-t', type=int, required=True,
                    choices=[0, 1, 2, 3, 4],
                    help="Ranking method: 0=baseline, 1=adaptive_k (seeds), 2=adaptive_k (all), 3=unsupervised, 4=disease-gene clustering")

parser.add_argument('--jobs', '-j', type=int, required=True, help="Number of jobs to run in parallel")
parser.add_argument('--folds', '-f', type=int, required=True, help="Cross validation folds (2-fold, 3-fold, etc.)")

def main(network_path, model_path, disease, partition_name, source, choice, jobs, folds):
    start_time = time()

    _init_model(model_path)

    source = source.lower()
    partition_folder = "../cross_validation/{0}/partitions/{1}/{2}".format(source, partition_name, disease)

    def extract_index(path):
        # grab the last number in the filename
        match = re.search(r'(\d+)(?=\D*$)', os.path.basename(path))
        return int(match.group(1))

    # get a list of all disease-gene files 
    gene_files = sorted(
        (os.path.join(partition_folder, e.name) for e in os.scandir(partition_folder) if "new_seeds" in e.name),
        key=extract_index
    )
    non_gene_files = sorted(
        (os.path.join(partition_folder, e.name) for e in os.scandir(partition_folder) if "non_seeds" in e.name),
        key=extract_index
    )
    n = len(gene_files)

    method_mapping = {0: "STRING_baseline",
                      1: "adaptive_k_cc",
                      2: "adaptive_k_cc_all_genes",
                      3: "network_clustering",
                      4: "disease_clustering_wip"}

    method = method_mapping[choice]

    # Make a lightweight job list
    jobspecs = [(network_path, genelist_path, i) for i, genelist_path in enumerate(gene_files)]

    with Parallel(
        n_jobs=jobs,
        backend="loky", # threading
        batch_size=1,
        pre_dispatch="2*n_jobs",

    ) as parallel:
        rankings = parallel(delayed(_rank_from_paths)(choice, npath, gpath, indx) for (npath, gpath, indx) in jobspecs)

    ranking_out_folder = "../cross_validation/{0}/rankings/{1}/{2}/".format(source, method, disease)
    label_out_folder = "../cross_validation/{0}/labels/{1}/{2}/".format(source, method, disease)

    for indx, ranking in rankings:
        # save each ranking to cross_validation/rankings
        ranking_file_name = "{0}_{1}_cross_validation_{2}.out".format(folds, disease, indx)
        ranking_out_file = open(ranking_out_folder + ranking_file_name, 'w')

        # save each ranking labels to cross_validation/labels
        label_file_name = "{0}_validation_labels_{1}".format(disease, indx)
        label_out_file = open(label_out_folder + label_file_name, 'w')

        # get seeds and nonseeds for current ranking
        with open(gene_files[indx]) as f:
            seeds = {line.strip() for line in f}

        with open(non_gene_files[indx]) as f:
            nonseeds = {line.strip() for line in f}

        for gene, score in ranking:
            name = gene.name
            ranking_out_file.write("{0}\t{1}\n".format(name, score))
            if name in seeds:
                continue

            label_out_file.write("{0} {1} {2}\n".format(name, score, 1 if name in nonseeds else 0))
            
        ranking_out_file.close()
        label_out_file.close()

    end_time = time()
    print("Total time:", end_time-start_time)


if __name__ == "__main__":
    args = parser.parse_args()
    main(args.network, args.model, args.disease, args.partition, args.source, args.method, args.jobs, args.folds)
