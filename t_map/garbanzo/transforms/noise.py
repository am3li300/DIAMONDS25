from copy import deepcopy
import networkx as nx
from networkx.algorithms import tree
import numpy as np


def gen_noise(graph: nx.Graph, p: float = 0.05) -> nx.Graph:
    """
    Generate noise on the graph.
    :param graph:
    :return:

    # Given a percentage p of the edges to be added/removed in the network
    # 1) If p is less than 0 (removal):-> preserve the maximum spanning tree. For the remaining edges not in the tree,
    # sample from the uniform distribution $U$, for each edge in the network. If the u_{ij}
    # sampled for edge (i, j) is less than p, then remove.
    # If p is greater than 0 (addition) -> randomly generate [(p, q, w)] 
    """
    graph = deepcopy(graph)
    if p < 0:
        """
        Edge removal operation:
        """
        mst = tree.maximum_spanning_edges(graph, algorithm="prim", data=False)
        filtered_edges = graph.copy()
        filtered_edges.remove_edges_from(list(mst))
        for e in filtered_edges.edges():
            pval = np.random.uniform()
            if pval < np.abs(p):
                graph.remove_edge(*e)
    else:
        edges_to_add = int(p * graph.number_of_edges())
        graph_nodes = list(graph.nodes())
        rand_tuples = np.random.randint(
            len(graph_nodes), size=(2 * edges_to_add, 2))
        count = 0
        for p_, q_ in rand_tuples:
            p, q = [graph_nodes[p_], graph_nodes[q_]]
            if not graph.has_edge(p, q):
                count += 1
                graph.add_edge(p, q, weight=np.random.uniform())
                if count >= edges_to_add:
                    break
    return graph
