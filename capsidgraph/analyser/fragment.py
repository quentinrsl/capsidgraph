import networkx as nx
import random
from typing import Dict, List


def _probability_fragment(G: nx.Graph, settings: Dict) -> nx.Graph:
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


def _energy_edges_fragment(G: nx.Graph, settings: dict) -> nx.Graph:
    # Get the attributes of the edges of the graph
    bond_energy = nx.get_edge_attributes(G, "energy")
    edges = list(G.edges)
    # Compute edge probability weights
    weights = []
    for e in edges:
        weights.append(1 / bond_energy[e])
    min_bond_energy = 1 / max(weights)  # Energy of the weakest bond
    G_ = G.copy()
    remaining_nodes = list(G_.edges)
    removed_edges = []
    remainig_edges_weights = weights.copy()
    energy = settings["fragmentation"]
    while energy > min_bond_energy:
        i = None
        # randomly pick a bond if there are any left, if the energy of the bond picked is higher than the percolationEnergy left, pick another one.
        # This loop has to stop because the energy is bigger than the minimum bond energy
        while (i == None or energy <= bond_energy[remaining_nodes[i]]) and len(
            remaining_nodes
        ) > 0:
            [i] = random.choices(
                list(range(len(remaining_nodes))), weights=remainig_edges_weights, k=1
            )

        # If a bond has been picked, remove it and make subsequent updates
        if i != None:
            energy -= bond_energy[remaining_nodes[i]]
            removed_edges.append(remaining_nodes[i])
            del remaining_nodes[i]
            del remainig_edges_weights[i]
            # update weakest bond energy
            if len(remainig_edges_weights) > 0:
                min_bond_energy = 1 / max(remainig_edges_weights)
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
    removed_neighbours = []
    for n1, n2, edgeAttributes in G.edges(remaining_nodes[i], True):
        neighbour = n1 if n1 != remaining_nodes[i] else n2
        if neighbour not in removed_nodes:
            neighbourIndex = remaining_nodes.index(neighbour)
            G.nodes[neighbour]["energy"] -= edgeAttributes["energy"]

            # Put the neighbours to be removed in a separate list to keep the current node to remove at the index i in the remainingNode list
            if G.nodes[neighbour]["energy"] > 1e-15:
                weights[neighbourIndex] = abs(1 / G.nodes[neighbour]["energy"])
            else:
                removed_neighbours.append(remaining_nodes[neighbourIndex])

    # remove node
    del remaining_nodes[i]
    del weights[i]
    # Remove 0 energy neighbours
    for neighbour in removed_neighbours:
        neighbourIndex = remaining_nodes.index(neighbour)
        del remaining_nodes[neighbourIndex]
        del weights[neighbourIndex]
        removed_nodes.append(neighbour)


def _energy_nodes_fragment(G: nx.Graph, settings: Dict) -> nx.Graph:
    weights = []
    for e in G.nodes:
        weights.append(1 / G.nodes[e]["energy"])
    G_ = G.copy()
    remaining_nodes = list(G_.nodes)
    removed_nodes = []
    energy = settings["fragmentation"]

    min_node_energy = 1 / max(weights)

    while energy + 1e-15 > min_node_energy:
        i = None
        while (
            i == None or energy + 1e-15 < G_.nodes[remaining_nodes[i]]["energy"]
        ) and len(remaining_nodes) > 0:
            [i] = random.choices(
                list(range(len(remaining_nodes))), weights=weights, k=1
            )
        if i != None:
            energy -= G_.nodes[remaining_nodes[i]]["energy"]
            removed_nodes.append(remaining_nodes[i])
            # update the energy / probability weights of the neighbouring nodes
            _remove_node(G_, remaining_nodes, removed_nodes, weights, i)

            # update weakest bond energy
            if len(weights) > 0:
                min_node_energy = 1 / max(weights)
        else:
            # All edges have been removed
            break

    for node in removed_nodes:
        G_.remove_node(node)

    return G_
