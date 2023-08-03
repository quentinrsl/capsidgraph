import numpy as np
import matplotlib.pyplot as plt
from capsidgraph.generator import (
    icosahedral_patterns,
    create_icosahedral_face_edges,
    create_icosahedral_capsid_graph,
)
from capsidgraph.analyser import get_fragmentation_probability_random_edge_removal
from capsidgraph.analyser import get_fragmentation_probability_random_node_removal

"""
This example demonstrate how one can use the get_fragmentation_probability_random_node_removal and get_fragmentation_probability_random_edge_removal functions to study a weighted interaction network under random edge and node removal.
This example computes the percolation threshold for different removal probability and plots the result.
"""

fragmentation_type = "nodes"  # "nodes" or "edges", type of removal
processes = 4  # Number of processes
iterations = 10000  # Number of iterations
pointNumber = 20  # Number of points in the plot

X = np.linspace(0, 1, pointNumber)

h = 1
k = 1
[edges, Tx, Ty, Tscale] = icosahedral_patterns.PATTERN_666
face_edges, axis = create_icosahedral_face_edges(edges, Tx, Ty, h, k)
G = create_icosahedral_capsid_graph(face_edges, axis)
print("T=", Tscale * (h * h + k * k + h * k))


if __name__ == "__main__":
    Y = []
    for p in X:
        if fragmentation_type == "nodes":
            Y.append(get_fragmentation_probability_random_node_removal(G, p, iterations, debug=True, process_number=processes))
        elif fragmentation_type == "edges":
            Y.append(get_fragmentation_probability_random_edge_removal(G, p, iterations, debug=True, process_number=processes))
    print(Y)
    plt.plot(X, Y)
    plt.show()
