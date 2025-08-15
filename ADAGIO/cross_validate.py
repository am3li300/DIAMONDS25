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
from heapq import heapify, heappush, heappop

import networkx as nx
from scipy.sparse import csr_matrix
import markov_clustering as mc

import igraph as ig

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pcst_fast
import numpy as np
from joblib import Parallel, delayed

from scipy.ndimage import gaussian_filter1d
from scipy.signal import argrelextrema

from collections import defaultdict
import dill as pickle

parser = argparse.ArgumentParser()
parser.add_argument('--network', '-n', type=str, required=True, help="Path to edgelist, assumes tab separated.")

def main(network_path):
    # load model
    model_file = open(input("Enter path to pickled model file: "), "rb")
    model = pickle.load(model_file)

    # default partition folder for cross validation is cross_validation/partitions/STRING
    disease = input("Enter name of disease: ")
    change_flag = input("STRING is default partition. Use different partition (Y/N)? ")
    if change_flag[0].lower() == 'y':
        partition_folder = "cross_validation/partitions/{0}/{1}".format(input("Enter partition folder name: "), disease)

    else:
        partition_folder = "cross_validation/partitions/STRING/{0}".format(disease)

    # get a list of all disease-gene files 
    gene_files = [entry.name for entry in os.scandir(partition_folder) if entry.is_file()]
    graphs = [EdgeListGarbanzo(network_path, genelist_path) for genelist_path in gene_files]

    # use joblib to cross-validate in parallel
    jobs = max(1, int(input("Enter number of jobs to run in parallel: ")))

    print("-------------------------------------------")
    print("              Select Method                ")
    print("-------------------------------------------")
    print(" 0) Baseline (constant k=20)")
    print(" 1) Adaptive k (disease-genes only)")
    print(" 2) Adaptive k (all genes)")
    print(" 3) Unsupervised network clustering")
    print(" 4) Disease-gene clustering")
    print("-------------------------------------------")

    method = int(input())
    if method == 0:
        def baseline_helper(graph):
            return sorted(list(model.prioritize(graph.genes, graph.graph)), key=lambda x: -x[1])

        rankings = Parallel(n_jobs=jobs)(delayed(baseline_helper)(graph) for graph in graphs)

    elif method == 1:
        def adaptive_k_disease_helper(graph):
            return sorted(list(model.david_prioritize_2(graph.genes, graph.graph)), key=lambda x: -x[1])

        rankings = Parallel(n_jobs=jobs)(delayed(adaptive_k_disease_helper)(graph) for graph in graphs)

    elif method == 2:
        def adaptive_k_all_helper(graph):
            return sorted(list(model.david_prioritize(1000, graph.graph, True)), key=lambda x: -x[1])

        rankings = Parallel(n_jobs=jobs)(delayed(adaptive_k_all_helper)(graph) for graph in graphs)

    elif method == 3:
        pass
    
    elif method == 4:
        pass

    else:
        print("Invalid Method")
        return 

    # save each ranking to cross_validation/rankings
    