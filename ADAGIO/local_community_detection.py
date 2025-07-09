"""
Local Community Detection Algorithms

Local Community Detection (LCD) aims to detect one or a few communities
starting from certain source nodes in the network. This differs from Global
Community Detection (GCD), which aims to partition an entire network into
communities.

LCD is often useful when only a portion of the graph is known or the
graph is large enough that GCD is infeasible.

References
----------
.. [1] Baltsou, Georgia, Konstantinos Christopoulos, and Konstantinos Tsichlas.
   Local community detection: A survey. IEEE Access 10 (2022): 110701-110726.
   https://doi.org/10.1109/ACCESS.2022.3213980
"""

__all__ = ["greedy_source_expansion"]

import networkx as nx


def _clauset_greedy_source_expansion(G, *, source, cutoff=None):
    if cutoff is None:
        cutoff = float("inf")
    C = {source}
    B = {source}
    U = set(G[source]) - C
    T = {frozenset([node, nbr]) for node in B for nbr in G.neighbors(node)}
    I = {edge for edge in T if all(node in C for node in edge)}

    R_value = 0
    while len(C) < cutoff:
        if not U:
            break

        max_R = 0
        best_node = None
        best_node_B = best_node_T = best_node_I = set()

        for v in U:
            R_tmp, B_tmp, T_tmp, I_tmp = _calculate_local_modularity_for_candidate(
                G, v, C, B, T, I
            )
            if R_tmp > max_R:
                max_R = R_tmp
                best_node = v
                best_node_B = B_tmp
                best_node_T = T_tmp
                best_node_I = I_tmp

        if best_node is None:
            break

        C.add(best_node)
        U.update(set(G[best_node]) - C)
        U.discard(best_node)
        B = best_node_B
        T = best_node_T
        I = best_node_I

        if max_R < R_value:
            break
        R_value = max_R

    return C


def _calculate_local_modularity_for_candidate(G, v, C, B, T, I):
    C_tmp = C | {v}
    B_tmp = B.copy()
    T_tmp = T.copy()
    I_tmp = I.copy()
    removed_B_nodes = set()

    for nbr in G[v]:
        if nbr not in C_tmp:
            B_tmp.add(v)
            T_tmp.add(frozenset([v, nbr]))

        if nbr in B:
            if all(nbr_nbr in C_tmp for nbr_nbr in G[nbr]):
                B_tmp.remove(nbr)
                removed_B_nodes.add(nbr)

        if nbr in C_tmp:
            I_tmp.add(frozenset([v, nbr]))

    for removed_node in removed_B_nodes:
        for removed_node_nbr in G[removed_node]:
            if removed_node_nbr not in B_tmp:
                T_tmp.discard(frozenset([removed_node_nbr, removed_node]))
                I_tmp.discard(frozenset([removed_node_nbr, removed_node]))

    R_tmp = len(I_tmp) / len(T_tmp) if T_tmp else 1
    return R_tmp, B_tmp, T_tmp, I_tmp


ALGORITHMS = {
    "clauset": _clauset_greedy_source_expansion,
}


def greedy_source_expansion(G, *, source, cutoff=None, method="clauset"):
    """Find the local community around a source node using Greedy Source Expansion.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    source : node
        The source node from which the community expansion begins.
    cutoff : int, optional
        The maximum number of nodes to include in the community.
    method : str, optional
        The algorithm to use (default 'clauset').

    Returns
    -------
    set
        A set of nodes representing the local community.
    """
    try:
        algo = ALGORITHMS[method]
    except KeyError as e:
        raise ValueError(f"{method} is not a valid choice for an algorithm.") from e

    return algo(G, source=source, cutoff=cutoff)
