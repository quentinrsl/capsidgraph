import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from capsidgraph.generator import (
    icosahedral_patterns,
    create_icosahedral_face_edges,
    create_icosahedral_capsid_graph,
)
from capsidgraph.analyser import get_fragmentation_probability_random_edge_removal
from capsidgraph.analyser import get_fragmentation_probability_random_node_removal

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


def worker(p):
    print("Computing for p=", p)
    if fragmentation_type == "nodes":
        return get_fragmentation_probability_random_node_removal(G, p, iterations)
    elif fragmentation_type == "edges":
        return get_fragmentation_probability_random_edge_removal(G, p, iterations)


if __name__ == "__main__":
    with Pool(processes) as pool:
        Y = pool.map(worker, X)
    print(Y)
    plt.plot(X, Y)
    plt.show()
