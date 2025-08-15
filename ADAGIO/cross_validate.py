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


_MODEL = None

def _init_model(model_path):
    global _MODEL
    if _MODEL is None:
        with open(model_path, "rb") as f:
            _MODEL = pickle.load(f)

def _rank_baseline(graph):
    # uses the global _MODEL already loaded in the worker
    ranking = sorted(list(_MODEL.prioritize(graph.genes, graph.graph)), key=lambda x: -x[1])
    return ranking

def _rank_adaptive_disease(graph):
    ranking = sorted(list(_MODEL.david_prioritize_2(graph.genes, graph.graph)), key=lambda x: -x[1])
    return ranking

def _rank_adaptive_all(graph):
    ranking = sorted(list(_MODEL.david_prioritize(1000, graph.graph, True)), key=lambda x: -x[1])
    return ranking

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

def main(network_path, model_path):
    start_time = time()

    _init_model(model_path)
    disease = input("Enter name of disease: ")

    # default partition folder for cross validation is cross_validation/partitions/STRING/disease
    change_flag = input("STRING is default partition. Use different partition (Y/N)? ")
    sub = input("Enter partition folder name: ") if change_flag[0].lower() == 'y' else "STRING"
    partition_folder = "../cross_validation/partitions/{0}/{1}".format(sub, disease)

    # two-fold/three-fold, etc.
    folds = int(input("Enter number of folds for validation: "))

    # get a list of all disease-gene files 
    gene_files = sorted(os.path.join(partition_folder, e.name) for e in os.scandir(partition_folder) if "new_seeds" in e.name)
    n = len(gene_files)

    # use joblib to cross-validate in parallel
    jobs = min(n, max(1, int(input("Enter number of jobs to run in parallel: "))))

    print("-------------------------------------------")
    print("              Select Method                ")
    print("-------------------------------------------")
    print(" 0) Baseline (constant k=20)")
    print(" 1) Adaptive k (disease-genes only)")
    print(" 2) Adaptive k (all genes)")
    print(" 3) Unsupervised network clustering")
    print(" 4) Disease-gene clustering")
    print("-------------------------------------------")

    choice = int(input())
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
        backend="loky",
        batch_size=1,
        pre_dispatch="2*n_jobs",

    ) as parallel:
        rankings = parallel(delayed(_rank_from_paths)(choice, npath, gpath, indx) for (npath, gpath, indx) in jobspecs)

    # save each ranking to cross_validation/rankings
    out_folder = "../cross_validation/rankings/{0}/{1}/".format(method, disease)
    for indx, ranking in rankings:
        file_name = "{0}_{1}_cross_validation_{2}.out".format(folds, disease, indx)
        out_file = open(out_folder + file_name, 'w')
        for gene, score in ranking:
            out_file.write("{0} {1}\n".format(gene.name, score))

        out_file.close()

    """
    add label functionality?
    """

    end_time = time()
    print("Total time:", end_time-start_time)


if __name__ == "__main__":
    args = parser.parse_args()
    main(args.network, args.model)






    """
    if choice == 0:
        method = "STRING_baseline"
        def baseline_helper(graph, i):
            ranking = sorted(list(model.prioritize(graph.genes, graph.graph)), key=lambda x: -x[1])
            return (i, ranking)

        rankings = Parallel(n_jobs=jobs)(delayed(baseline_helper)(graphs[i], i) for i in range(n))

    elif choice == 1:
        method = "adaptive_k_cc"
        def adaptive_k_disease_helper(graph, i):
            ranking = sorted(list(model.david_prioritize_2(graph.genes, graph.graph)), key=lambda x: -x[1])
            return (i, ranking)

        rankings = Parallel(n_jobs=jobs)(delayed(adaptive_k_disease_helper)(graphs[i], i) for i in range(n))

    elif choice == 2:
        method = "adaptive_k_cc_all_genes"
        def adaptive_k_all_helper(graph, i):
            ranking = sorted(list(model.david_prioritize(1000, graph.graph, True)), key=lambda x: -x[1])
            return (i, ranking)

        rankings = Parallel(n_jobs=jobs)(delayed(adaptive_k_all_helper)(graphs[i], i) for i in range(n))

    elif choice == 3:
        method = "network_clustering"
        pass
    
    elif choice == 4:
        method = "disease_clustering_wip"
        pass

    else:
        print("Invalid Method")
        return 
    """