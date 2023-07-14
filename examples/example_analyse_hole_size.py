import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

from capsidgraph.generator import (
    icosahedral_patterns,
    create_icosahedral_face_edges,
    create_icosahedral_capsid_graph,
)
from capsidgraph.analyser import (
    energy_nodes_fragment,
    init_nodes_energy,
    get_hole_size_distribution
)

iterations = 10000
E = 1

h = 1
k = 1
[edges, Tx, Ty, Tscale] = icosahedral_patterns.PATTERN_666
face_edges, axis = create_icosahedral_face_edges(edges, Tx, Ty, h, k)
G = create_icosahedral_capsid_graph(face_edges, axis)
nx.set_edge_attributes(G, 1 / len(G.edges), "energy")
init_nodes_energy(G)
res = get_hole_size_distribution(G, iterations, energy_nodes_fragment, {"fragmentation":E, "fragmentation_type":"nodes"})

plt.bar(list(range(len(res))), res)
plt.show()