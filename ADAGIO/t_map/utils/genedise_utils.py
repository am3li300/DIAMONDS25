import networkx as nx
import pandas as pd
from pathlib import Path
from t_map.gene.gene import Gene
from typing import Dict, List, Literal, Optional, Tuple, Union
from collections import defaultdict
from dataclasses import dataclass
from t_map.garbanzo.networkx import Networkx
from glob import glob
from t_map.feta.feta import Feta, PreComputeFeta

CLS_TYPE = Union[Literal["genetic"], Literal["drugs"]]


@dataclass
class Split:
    rep: int
    disease_name: str
    data_type: CLS_TYPE
    training_folds: List[str]
    validation_folds: List[str]


def run_genedise(graph_path: str,
                 splits_path: str,
                 output_path: str,
                 model: Union[type[Union[Feta, PreComputeFeta]],
                              Feta,
                              PreComputeFeta],
                 pickle_path: Optional[str] = None,
                 variant: str="none") -> None:

    graph = load_graph(graph_path)

#     graph = reweight_graph_by_tissue(
#         graph, "../data/is_brain.tsv", getter=getter)

    splits = glob(splits_path + "/*.csv")
    split_dict = restructure_splits(splits, splits_path)

    for rep, splits in split_dict.items():
        print(f"Running rep {rep} {len(splits)}")
        for split in splits:
            run_split(split, graph, model, output_path, pickle_path, variant)


def getter(string: str):
    name, _, is_expressed = string.split("\t")
    return name, int(is_expressed)


def run_split(split: Split, graph: nx.Graph,
              model: Union[type[Union[Feta, PreComputeFeta]],
                           Feta,
                           PreComputeFeta],
              output_path: str,
              pickle_path: Optional[str] = None,
              variant:str="none") -> None:
    needs_to_load_pickle = pickle_path is not None
    training_folds = split.training_folds
    validation_folds = split.validation_folds
    if isinstance(model, Feta) or isinstance(model, PreComputeFeta):
        needs_to_load_pickle = False
    else:
        model = model()

    validation_score_per_fold = []
    for i, (training_fold, validation_fold) in enumerate(zip(training_folds,
                                                             validation_folds)):
        training_genes = read_genelist(training_fold)
        data = Networkx(graph, training_genes)

        if needs_to_load_pickle:
            try:
                model = model.load(pickle_path)
            except OSError:
                model.setup(data.graph)
            predictions = model(data.genes, model.graph, variant=variant)
        else:
            predictions = model(data.genes, data.graph, variant=variant)

        # Extract validation genes from predictions
        validation_scores = [(g, s) for g, s in predictions]
        validation_score_per_fold.append(validation_scores)
        output_results(validation_score_per_fold,
                       f"{output_path}/{split.disease_name}_{split.rep}_{split.data_type}_Fold{i}.csv")
        validation_score_per_fold = []
        # If PreComputeFeta, save the model
        #if isinstance(model, PreComputeFeta) and needs_to_load_pickle:
        #    needs_to_load_pickle = False
        #    model.save(pickle_path)


def output_results(results: List[List[Tuple[Gene, float]]],
                   out_path: str) -> None:
    import itertools
    results = itertools.chain.from_iterable(results)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for (gene, score) in results:
            f.write(f"{gene.name},{score}\n")


def read_genelist(path: str) -> List[Gene]:
    df = pd.read_csv(path, sep=",", skiprows=1, names=["gene", "is_disease"])
    df = df[df["is_disease"] == 1]
    gene_list = [Gene(g) for g in df["gene"].tolist()]
    return gene_list


def restructure_splits(splits: List[str], splits_path: str) -> Dict[int, List[Split]]:
    rep_to_files = defaultdict(list)
    rep_to_split = dict()
    prefix = splits_path + "/"
    for split in splits:
        rep = split.removeprefix(prefix).split(".")[0]
        rep_to_files[rep].append(split)

    for rep in rep_to_files.keys():
        files = rep_to_files[rep]
        validation_fold = [f for f in files if "_val" in f]
        training_fold = [f for f in files if "_val" not in f]
        int_rep = int(rep[3:])
        rep_to_split[int_rep] = generate_split_for_fold(
            int_rep, validation_fold, training_fold, prefix)

    return rep_to_split


def generate_split_for_fold(fold: int, validation_fold: List[str],
                            training_fold: List[str], prefix: str) -> List[Split]:
    # A dictionary that takes the type of the split
    # and returns a tuple of training and validation folds
    ty = Dict[CLS_TYPE, Tuple[List[str], List[str]]]

    disease_to_data_type_to_fold: ty = dict()

    for tf, vf in zip(sorted(training_fold), sorted(validation_fold)):
        _, disease_name, disease_type, *_ = tf.removeprefix(prefix).split("_")
        if disease_type not in disease_to_data_type_to_fold:
            disease_to_data_type_to_fold[disease_type] = dict()
        if disease_name not in disease_to_data_type_to_fold[disease_type]:
            disease_to_data_type_to_fold[disease_type][disease_name] = ([], [])

        # Add the training and validation folds to correct split list
        disease_to_data_type_to_fold[disease_type][disease_name][0].append(tf)
        disease_to_data_type_to_fold[disease_type][disease_name][1].append(vf)

    splits = []
    for disease_type, disease_to_fold in disease_to_data_type_to_fold.items():
        for disease_name, (training_folds,
                           validation_folds) in disease_to_fold.items():
            splits.append(
                Split(fold, disease_name, disease_type.removesuffix(".csv"), training_folds,
                      validation_folds))
    return splits


def load_and_relabel_graph(path: str) -> nx.Graph:
    graph = nx.read_graphml(path)
    relabeling_dict = dict()

    for node, data in graph.nodes(data=True):
        relabeling_dict[node] = data["name"]
    graph = nx.relabel_nodes(graph, relabeling_dict)

    return graph


def load_graph(path: str) -> nx.Graph:
    # load tab seperated edge list from path
    edge_list = pd.read_csv(path, sep="\t")
    # convert pandas dataframe to networkx graph
    graph = nx.from_pandas_edgelist(
        edge_list, source="src", target="dst", edge_attr=True)
    return graph
