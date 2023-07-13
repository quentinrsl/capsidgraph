import networkx as nx
from typing import List


# intialize weights on nodes
# Create weight on nodes based on the sum of the energy of the bonds attached to it
# Ignore nodes that are in the blacklist (act if they were removed)hist
def init_nodes_energy(G: nx.Graph) -> List[float]:
    # Compute node weights
    for node in G.nodes():
        energy = 0
        for edge in G.edges(node):
            energy += G[edge[0]][edge[1]]["energy"]
        G.nodes[node]["energy"] = energy
    # Compute edge probability weights
    nodes_probability_weights = []
    for e in G.nodes:
        nodes_probability_weights.append(1 / G.nodes[e]["energy"])
        energy = 0
        for edge in G.edges(node):
            energy += G[edge[0]][edge[1]]["energy"]
    return nodes_probability_weights
