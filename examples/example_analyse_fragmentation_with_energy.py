import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from capsidgraph.generator import (
    icosahedral_patterns,
    create_icosahedral_face_edges,
    create_icosahedral_capsid_graph,
)
from capsidgraph.analyser import get_fragmentation_probability_energy_edge_removal
from capsidgraph.analyser import get_fragmentation_probability_energy_node_removal


"""
This example demonstrate how one can use the get_fragmentation_probability_energy_node_removal and get_fragmentation_probability_energy_edge_removal functions to study a weighted interaction network under random edge and node removal.
This example computes the energy percolation threshold for different fraction of energy removed and plots the result.
"""

fragmentation_type = "nodes"  # "nodes" or "edges", type of removal
processes = 12  # Number of processes
iterations = 10000  # Number of iterations
pointNumber = 20  # Number of points in the plot
fb = 0.5  #fraction of b edges' energy compared to a edges
fc = 1 #fraction of c edges' energy compared to a edges

X = np.linspace(0, 1, pointNumber)

face_edges = [((1,0),(1,1)), #c
    ((1,0),(1,-1)),#b
    ((1,0),(0,1)),#b
    ((1,0),(2,-1)),#c
    ((1,0),(0,0)),#a
    ((1,0),(2,0)),#c
    ((1,1),(2,0)),#b
    ((1,-1),(2,-1)),#c
    ((2,0),(2,1)),#a
    ((2,0),(2,-1)),#c
    ((2,0),(3,-1)),#c
    ((2,0),(3,0)),#b
    ((2,-1),(2,-2)),#b
    ((2,-1),(3,-2)),#a
    ((2,-1),(3,-1))#b
]
axis = ((0, 0), (2, 1), (3, -2))
Ea = 1 / (60 + 60 * fb + 90 * fc)
Eb = Ea * fb
Ec = Ea * fc
energy =   [Ec,Eb,Eb,Ec,Ea,Ec,Eb,Ec,Ea,Ec,Ec,Eb,Eb,Ea,Eb]
G = create_icosahedral_capsid_graph(face_edges, axis, energy)

def worker(p):
    print("Computing for p=", p)
    if fragmentation_type == "nodes":
        return get_fragmentation_probability_energy_node_removal(G, p, iterations,debug=True)
    elif fragmentation_type == "edges":
        return get_fragmentation_probability_energy_edge_removal(G, p, iterations,debug=True)


if __name__ == "__main__":
    with Pool(processes) as pool:
        Y = pool.map(worker, X)
    print(Y)
    plt.plot(X, Y)
    plt.show()
