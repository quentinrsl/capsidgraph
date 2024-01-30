import networkx as nx
import random
from typing import Dict, List


def probability_fragment(G: nx.Graph, settings: Dict) -> nx.Graph:
    """
    Fragment the graph G by removing each node or edge with a gibe probability

    Parameters
    ----------
    G : nx.Graph
        The graph to fragment
    settings : Dict
        The settings of the fragmentation.

        The `fragmentation` entry is the probability of removal. Its value is a float between 0 and 1.

        The `fragmentation_type` determines whether to remove nodes or edges. Its value can be `"edges"` or `"nodes"`

    Returns
    -------
    nx.Graph
        The fragmented graph
    """
    G_ = G.copy()
    p = settings["fragmentation"]
    fragmentation_type = settings["fragmentation_type"]
    if fragmentation_type == "nodes":
        # Remove each nodes with probability p
        for node in G.nodes:
            if random.random() < p:
                G_.remove_node(node)
    elif fragmentation_type == "edges":
        # Remove each edge with probability p
        for a, b in G.edges:
            if random.random() < p:
                G_.remove_edge(a, b)
    return G_


def strength_edges_fragment(G: nx.Graph, settings: dict) -> nx.Graph:
    """
    Fragment the graph G by randomly removing edges until the strength of the strength removed from the graph is equal to a given value.
    A probability weight is assigned to each edge, inversly proportional to its `strength` attribute.
    The strength of the graph is the sum of the strength of the edges.

    Parameters
    ----------
    G : nx.Graph
        The graph to fragment
    settings : Dict
        The settings of the fragmentation.

        The `fragmentation` entry is the strength to remove from the graph. Its value is a float.

    Returns
    -------
    nx.Graph
        The fragmented graph
    """
    # Get the attributes of the edges of the graph
    bond_strength = nx.get_edge_attributes(G, "strength")
    edges = list(G.edges)
    # Compute edge probability weights
    weights = []
    for e in edges:
        weights.append(1 / bond_strength[e])
    min_bond_strength = 1 / max(weights)  # strength of the weakest bond
    G_ = G.copy()
    remaining_nodes = list(G_.edges)
    removed_edges = []
    remainig_edges_weights = weights.copy()
    strength = settings["fragmentation"]
    while strength > min_bond_strength:
        i = None
        # randomly pick a bond if there are any left, if the strength of the bond picked is higher than the percolationStrength left, pick another one.
        # This loop has to stop because the strength is bigger than the minimum bond strenth
        while (i == None or strength <= bond_strength[remaining_nodes[i]]) and len(
            remaining_nodes
        ) > 0:
            [i] = random.choices(
                list(range(len(remaining_nodes))), weights=remainig_edges_weights, k=1
            )

        # If a bond has been picked, remove it and make subsequent updates
        if i != None:
            strength -= bond_strength[remaining_nodes[i]]
            removed_edges.append(remaining_nodes[i])
            del remaining_nodes[i]
            del remainig_edges_weights[i]
            # update weakest bond strength
            if len(remainig_edges_weights) > 0:
                min_bond_strength = 1 / max(remainig_edges_weights)
        else:
            # All edges have been removed
            break
    # Remove the edges from the graoh
    for a, b in removed_edges:
        G_.remove_edge(a, b)
    return G_


def _remove_node(
    G: nx.Graph,
    remaining_nodes: List[int],
    removed_nodes: List[int],
    weights: List[float],
    i: int,
):
    """
    Remove a node from the graph and update the strength of the neighbours
    """
    removed_neighbours = []
    for n1, n2, edgeAttributes in G.edges(remaining_nodes[i], True):
        neighbour = n1 if n1 != remaining_nodes[i] else n2
        if neighbour not in removed_nodes:
            neighbourIndex = remaining_nodes.index(neighbour)
            G.nodes[neighbour]["strength"] -= edgeAttributes["strength"]

            # Put the neighbours to be removed in a separate list to keep the current node to remove at the index i in the remainingNode list
            if G.nodes[neighbour]["strength"] > 1e-15:
                weights[neighbourIndex] = abs(1 / G.nodes[neighbour]["strength"])
            else:
                removed_neighbours.append(remaining_nodes[neighbourIndex])

    # remove node
    del remaining_nodes[i]
    del weights[i]
    # Remove 0 strength neighbours
    for neighbour in removed_neighbours:
        neighbourIndex = remaining_nodes.index(neighbour)
        del remaining_nodes[neighbourIndex]
        del weights[neighbourIndex]
        removed_nodes.append(neighbour)


def strength_nodes_fragment(G: nx.Graph, settings: Dict) -> nx.Graph:
    """
    Fragment the graph G by randomly removing nodes until the strength of the graph is less that a given value.
    A probability weight is assigned to each node, inversly proportional to the sum of its `strength` attribute.
    The strength of the graph is the sum of the strength of the edges.

    Parameters
    ----------
    G : nx.Graph
        The graph to fragment
    settings : Dict
        The settings of the fragmentation.

        The `fragmentation` entry is the strength to remove from the graph. Its value is a float.

    Returns
    -------
    nx.Graph
        The fragmented graph
    """
    weights = []
    for e in G.nodes:
        weights.append(1 / G.nodes[e]["strength"])
    G_ = G.copy()
    remaining_nodes = list(G_.nodes)
    removed_nodes = []
    strength = settings["fragmentation"]

    min_node_strength = 1 / max(weights)

    while strength + 1e-15 > min_node_strength:
        i = None
        while (
            i == None or strength + 1e-15 < G_.nodes[remaining_nodes[i]]["strength"]
        ) and len(remaining_nodes) > 0:
            [i] = random.choices(
                list(range(len(remaining_nodes))), weights=weights, k=1
            )
        if i != None:
            strength -= G_.nodes[remaining_nodes[i]]["strength"]
            removed_nodes.append(remaining_nodes[i])
            # update the strength / probability weights of the neighbouring nodes
            _remove_node(G_, remaining_nodes, removed_nodes, weights, i)

            # update weakest bond strength
            if len(weights) > 0:
                min_node_strength = 1 / max(weights)
        else:
            # All edges have been removed
            break

    for node in removed_nodes:
        G_.remove_node(node)

    return G_
